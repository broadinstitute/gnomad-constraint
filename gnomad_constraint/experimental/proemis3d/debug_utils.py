"""Debug utility functions for the Proemis3D pipeline."""

import logging
from typing import Dict, Optional

import hail as hl
import numpy as np

# Note: getAIC is imported inside functions that use it to avoid circular imports
# OE upper functions are imported inside functions that use them to avoid
# circular imports

logger = logging.getLogger("proemis3d_debug")
logger.setLevel(logging.INFO)

# ANSI color codes
RESET = "\033[0m"
RED = "\033[91m"  # High PAE (above cutoff), moderate constraint
BOLD = "\033[1m"
ORANGE = "\033[38;5;214m"  # Filtered residues (bright orange, 256-color mode)
# Bold, underline, green for minimum OE upper, very constrained
HIGHLIGHT = "\033[1m\033[4m\033[92m"
UNDERLINE = "\033[4m"
GREEN = "\033[92m"  # Good values (low AIC, very negative NLL)
YELLOW = "\033[93m"  # Moderate values
# Additional color codes for gradients
GREEN_IDX = "\033[92m"  # First in sorted array (closest)
YELLOW_IDX = "\033[93m"  # Middle
MAGENTA_IDX = "\033[95m"  # Last in sorted array (furthest)
BLUE = "\033[94m"  # Closest (bright blue)
CYAN = "\033[96m"  # Medium-close (cyan)
WHITE = "\033[97m"  # Furthest (white)


def _calculate_oe_upper(
    obs: int, exp: float, oe_upper_method: str = "gamma", alpha: float = 0.05
) -> Optional[float]:
    """
    Calculate OE upper confidence interval for debug output.

    :param obs: Observed count
    :param exp: Expected count
    :param oe_upper_method: Method to use ("gamma" or "chisq")
    :param alpha: Significance level (default 0.05)
    :return: OE upper value or None if calculation fails
    """
    if obs is None or exp is None or exp <= 0:
        return None

    try:
        # Lazy import to avoid circular dependency
        from gnomad_constraint.experimental.proemis3d.utils import (
            chisq_upper_ci,
            gamma_upper_ci,
        )

        oe_upper_func = gamma_upper_ci if oe_upper_method == "gamma" else chisq_upper_ci
        oe_upper_val = hl.eval(
            oe_upper_func(
                hl.literal(obs),
                hl.literal(exp),
                alpha,
            )
        )
        return float(oe_upper_val)
    except Exception:
        return None


def debug_print_pae_matrix_for_region(
    dist_mat_expr: hl.expr.ArrayExpression,
    center_residue_index_expr: hl.expr.Int32Expression,
    max_pae: float,
    title: str = "PAE matrix for region filtering",
    min_plddt: Optional[float] = None,
    plddt_cutoff_method: Optional[str] = None,
) -> Dict[int, str]:
    """
    Generate debug output showing PAE values between each residue and all subsequent residues
    in the sorted-by-distance array. This helps visualize how filter_on_pairwise_pae_in_region works.
    Returns a dictionary where keys are center residue indices and values are output strings.

    :param dist_mat_expr: Array expression with distance matrix entries (sorted by distance)
    :param center_residue_index_expr: Center residue index expression
    :param max_pae: Maximum allowed PAE value
    :param title: Title for the debug output
    :return: Dictionary mapping center residue indices to their debug output strings
    """
    # Dictionary to store output per center residue
    output_dict = {}

    # Get first key for debugging
    ht = dist_mat_expr._indices.source
    uniprot_id = None
    transcript_id = None
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.enst.collect()[0]

    # Filter and collect data
    _ht_debug = ht.annotate(
        dist_mat=dist_mat_expr, center_idx=center_residue_index_expr
    )
    _ht_debug = _ht_debug.filter(
        (_ht_debug.uniprot_id == uniprot_id) & (_ht_debug.enst == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return {}

    # Collect data for all center residues
    debug_data = _ht_debug.select("dist_mat", "center_idx").collect()

    if not debug_data:
        return {}

    # Group by center residue to show matrix for each
    center_groups = {}
    for row in debug_data:
        center_idx = row.center_idx
        if center_idx not in center_groups:
            center_groups[center_idx] = []
        center_groups[center_idx].append(row)

    # Show PAE matrix for each center residue
    for center_idx in sorted(center_groups.keys()):
        rows = center_groups[center_idx]
        # Use first row for this center (they should all have the same dist_mat)
        row = rows[0]
        dist_mat = row.dist_mat

        # Build output for this center residue
        center_output = f"    {BOLD}=== {title} ==={RESET}\n\n"
        center_output += f"        {BOLD}PAE matrix (row = from residue, column = to residue):{RESET}\n\n"
        center_output += f"        {BOLD}Color legend:{RESET}\n"
        center_output += f"            {BOLD}PAE colors:{RESET} Uncolored (≤{max_pae})  {RED}Red (>{max_pae}){RESET}\n"

        if min_plddt is not None:
            center_output += f"            {BOLD}pLDDT colors:{RESET} Uncolored (≥{min_plddt})  {RED}Red (<{min_plddt}){RESET}\n"

        center_output += f"            {BOLD}Filtered residues:{RESET} {ORANGE}Orange{RESET} (excluded from region)\n"

        center_output += f"\n"

        # Check if PAE data is available in dist_mat
        # Note: residue_i.pae might be a float (PAE from center) or an array (full PAE row)
        # Also check for pae_array field which contains the full PAE array for
        # each residue
        has_pae_data = (
            dist_mat
            and len(dist_mat) > 0
            and (hasattr(dist_mat[0], "pae") or hasattr(dist_mat[0], "pae_array"))
        )
        pae_is_array = False
        has_pae_array_field = False
        if has_pae_data:
            # Check if pae_array field exists (preferred for pairwise PAE)
            if hasattr(dist_mat[0], "pae_array") and dist_mat[0].pae_array is not None:
                try:
                    len(dist_mat[0].pae_array)
                    has_pae_array_field = True
                    pae_is_array = True
                except (TypeError, AttributeError):
                    pass
            # Fall back to checking pae field
            if (
                not has_pae_array_field
                and hasattr(dist_mat[0], "pae")
                and dist_mat[0].pae is not None
            ):
                try:
                    len(dist_mat[0].pae)
                    pae_is_array = True
                except (TypeError, AttributeError):
                    pae_is_array = False

        # Check if pLDDT data is available
        has_plddt_data = (
            dist_mat and len(dist_mat) > 0 and hasattr(dist_mat[0], "plddt")
        )

        # Get all residue indices for column headers
        all_residue_indices = [entry.residue_index for entry in dist_mat]

        # Determine which residues would be filtered by pLDDT (if pLDDT filtering
        # is used)
        plddt_filtered_residues = set()
        # Also track low pLDDT residues for exclude_low_plddt_from_stats (they
        # remain but are excluded from PAE checks)
        low_plddt_residues = set()
        if min_plddt is not None and has_plddt_data and plddt_cutoff_method:
            if plddt_cutoff_method == "truncate_at_first_low_plddt":
                # Find first residue with low pLDDT and mark all after it
                found_first_low = False
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if found_first_low:
                        plddt_filtered_residues.add(residue_idx)
                    elif (
                        hasattr(residue, "plddt")
                        and residue.plddt is not None
                        and residue.plddt < min_plddt
                    ):
                        found_first_low = True
                        # Don't add the first low one, but mark all subsequent ones
            elif plddt_cutoff_method == "remove_low_plddt_residues":
                # Mark residues with low pLDDT
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if (
                        hasattr(residue, "plddt")
                        and residue.plddt is not None
                        and residue.plddt < min_plddt
                    ):
                        plddt_filtered_residues.add(residue_idx)
            elif plddt_cutoff_method == "exclude_low_plddt_from_stats":
                # Mark residues with low pLDDT (they remain in region but are excluded
                # from PAE checks)
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if (
                        hasattr(residue, "plddt")
                        and residue.plddt is not None
                        and residue.plddt < min_plddt
                    ):
                        low_plddt_residues.add(residue_idx)
            # For "mask_low_confidence_plddt", residues are not filtered from the region

        # Simulate the filtering logic to determine which residues would be filtered by PAE
        # This mimics the scan operation in filter_on_pairwise_pae_in_region
        # The logic: add residue if (pae_array is None) OR (NOT any(PAE > max_pae to residues in region))
        # Note: pLDDT filtering happens FIRST in the actual code, so we simulate pLDDT filtering first,
        # then PAE filtering on the remaining residues
        pae_filtered_residues = set()
        region = []

        if has_pae_array_field:
            # First, apply pLDDT filtering to get the residues that would remain after pLDDT filtering
            # (if pLDDT filtering removes residues)
            remaining_after_plddt = []
            if plddt_cutoff_method in [
                "truncate_at_first_low_plddt",
                "remove_low_plddt_residues",
            ]:
                # Only consider residues that pass pLDDT filtering
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if residue_idx not in plddt_filtered_residues:
                        remaining_after_plddt.append(residue)
            else:
                # For exclude_low_plddt_from_stats, all residues remain (they're just
                # marked)
                remaining_after_plddt = dist_mat

            # Simulate the scan: build region incrementally from residues that passed pLDDT filtering
            # A residue is added if its PAE to ALL residues in the current region is <= max_pae
            # A residue is filtered if its PAE to ANY residue in the current region is > max_pae
            # Note: For exclude_low_plddt_from_stats, low pLDDT residues are still checked for PAE,
            # but when checking PAE from a new residue to the region, we skip low
            # pLDDT residues in the region
            for residue in remaining_after_plddt:
                residue_idx = residue.residue_index

                # Check if this residue should be added to the region based on PAE
                # Low pLDDT residues are still checked (they're not automatically added)
                should_add = True

                if hasattr(residue, "pae_array") and residue.pae_array is not None:
                    # Check PAE from this residue to each residue already in the region
                    # Skip low pLDDT residues in the region (they don't count for PAE
                    # filtering)
                    for region_residue in region:
                        region_residue_idx = region_residue.residue_index
                        # Skip low pLDDT residues when checking PAE (they don't block
                        # other residues)
                        if region_residue_idx in low_plddt_residues:
                            continue
                        try:
                            if region_residue_idx < len(residue.pae_array):
                                pae_to_region = residue.pae_array[region_residue_idx]
                                if pae_to_region > max_pae:
                                    # PAE exceeds threshold, so this residue should be
                                    # filtered
                                    should_add = False
                                    break
                        except (TypeError, AttributeError, IndexError):
                            # If we can't get PAE value, treat as should not add
                            # (conservative)
                            should_add = False
                            break
                # else: pae_array is None or doesn't exist
                # In Hail code: (residue.pae_array is None) | ... means if None, add it
                # So should_add = True is already set, which is correct

                if should_add:
                    region.append(residue)
                else:
                    pae_filtered_residues.add(residue_idx)
        else:
            # If no PAE array field, all residues that passed pLDDT filtering are in
            # the region
            if plddt_cutoff_method in [
                "truncate_at_first_low_plddt",
                "remove_low_plddt_residues",
            ]:
                for residue in dist_mat:
                    residue_idx = residue.residue_index
                    if residue_idx not in plddt_filtered_residues:
                        region.append(residue)
            else:
                # For exclude_low_plddt_from_stats or no pLDDT filtering, all residues
                # are in the region
                region = dist_mat

        # Determine which residues are in the final region (not filtered)
        region_residue_indices = {res.residue_index for res in region}

        # Create column header (use variable for backslash since f-strings can't
        # have backslashes)
        from_to_label = "From\\To"
        header = f"        {from_to_label:<8}"
        for res_idx in all_residue_indices:
            # Color orange if residue is not in the final region (filtered out)
            is_filtered = res_idx not in region_residue_indices
            if is_filtered:
                header += f"{ORANGE}{res_idx:>8}{RESET}"
            else:
                header += f"{res_idx:>8}"
        # Add pLDDT column header if pLDDT filtering is used
        if min_plddt is not None and has_plddt_data:
            header += f"{'pLDDT':>8}"
        center_output += header + "\n"
        header_sep_len = 8 + 8 * len(all_residue_indices)
        if min_plddt is not None and has_plddt_data:
            header_sep_len += 8
        center_output += "        " + "-" * header_sep_len + "\n"

        # For each residue in the sorted array, show its PAE to all residues
        for i, residue_i in enumerate(dist_mat):
            residue_i_idx = residue_i.residue_index
            row_str = f"        {residue_i_idx:>8}"

            # Check if this source residue is filtered (not in final region)
            # If filtered, all PAE values in this row should be "-" since they're not
            # used for filtering
            is_source_filtered = residue_i_idx not in region_residue_indices

            # For each column (target residue)
            for j, residue_j in enumerate(dist_mat):
                residue_j_idx = residue_j.residue_index

                if j <= i:
                    # Lower triangle or diagonal: show "-" or "0.0" for self
                    if j == i:
                        # row_str += f"{'0.0':>8}"
                        row_str += f"{'-':>8}"
                    else:
                        row_str += f"{'-':>8}"
                else:
                    # Upper triangle: show actual PAE value
                    # If source residue is filtered, show "-" (its PAE values don't
                    # count for filtering)
                    if is_source_filtered:
                        pae_str = f"{'-':>8}"
                    # For exclude_low_plddt_from_stats, show "-" if source residue
                    # (residue_i) has low pLDDT
                    elif (
                        plddt_cutoff_method == "exclude_low_plddt_from_stats"
                        and residue_i_idx in low_plddt_residues
                    ):
                        pae_str = f"{'-':>8}"
                    else:
                        pae_val = None
                        if has_pae_data:
                            # First try pae_array field (preferred for pairwise PAE)
                            if (
                                has_pae_array_field
                                and hasattr(residue_i, "pae_array")
                                and residue_i.pae_array is not None
                            ):
                                try:
                                    if residue_j_idx < len(residue_i.pae_array):
                                        pae_val = residue_i.pae_array[residue_j_idx]
                                except (TypeError, AttributeError, IndexError):
                                    pass
                            # Fall back to pae field if it's an array
                            elif (
                                pae_is_array
                                and hasattr(residue_i, "pae")
                                and residue_i.pae is not None
                            ):
                                try:
                                    if residue_j_idx < len(residue_i.pae):
                                        pae_val = residue_i.pae[residue_j_idx]
                                except (TypeError, AttributeError, IndexError):
                                    pass

                        # Format PAE value (ensure consistent width accounting for ANSI
                        # codes)
                        if pae_val is not None:
                            if pae_val > max_pae:
                                # Format with color, ensuring the visible part is 5
                                # chars wide
                                pae_str = f"{RED}{pae_val:5.1f}{RESET}"
                            else:
                                pae_str = f"{pae_val:5.1f}"
                        else:
                            # No PAE data available
                            pae_str = "  ?  "

                    # Calculate visible width (without ANSI codes) and pad accordingly
                    # ANSI codes don't count toward visible width, so we need to pad the
                    # visible content
                    visible_pae_str = pae_str.replace(RED, "").replace(RESET, "")
                    padding_needed = max(0, 8 - len(visible_pae_str))
                    row_str += " " * padding_needed + pae_str

            # Add pLDDT column if pLDDT filtering is used
            if min_plddt is not None and has_plddt_data:
                plddt_val = None
                if hasattr(residue_i, "plddt"):
                    plddt_val = residue_i.plddt

                if plddt_val is not None:
                    if plddt_val < min_plddt:
                        plddt_str = f"{RED}{plddt_val:>8.1f}{RESET}"
                    else:
                        plddt_str = f"{plddt_val:>8.1f}"
                else:
                    plddt_str = "     NA"
                row_str += plddt_str

            center_output += row_str + "\n"

        # Add pLDDT row if pLDDT filtering is used
        if min_plddt is not None and has_plddt_data:
            plddt_row = f"        {'pLDDT':<8}"
            for residue in dist_mat:
                residue_idx = residue.residue_index
                plddt_val = None
                if hasattr(residue, "plddt"):
                    plddt_val = residue.plddt

                if plddt_val is not None:
                    if plddt_val < min_plddt:
                        plddt_str = f"{RED}{plddt_val:>8.1f}{RESET}"
                    else:
                        plddt_str = f"{plddt_val:>8.1f}"
                else:
                    plddt_str = "     NA"
                plddt_row += plddt_str
            # Add pLDDT column value (same as row label, or could show average/NA)
            plddt_row += f"{'pLDDT':>8}"
            center_output += "        " + "-" * header_sep_len + "\n"
            center_output += "\n" + plddt_row + "\n"

        # Store output for this center residue
        output_dict[center_idx] = center_output

    return output_dict


def get_debug_oe_table_string(
    oe_list, plddt_lookup=None, min_plddt=None, pae_lookup=None, max_pae=None
):
    # Print OE array (per-residue observed/expected)
    output_string = ""

    # Determine which columns to include
    has_plddt = (
        min_plddt is not None and plddt_lookup is not None and len(plddt_lookup) > 0
    )
    has_pae = max_pae is not None and pae_lookup is not None and len(pae_lookup) > 0

    # Add color legend
    legend_parts = []
    if has_pae:
        legend_parts.append(
            f"        {BOLD}PAE colors:{RESET} Uncolored (≤{max_pae})  {RED}Red (>{max_pae}){RESET}"
        )
    if has_plddt:
        legend_parts.append(
            f"        {BOLD}pLDDT colors:{RESET} Uncolored (≥{min_plddt})  {RED}Red (<{min_plddt}){RESET}"
        )
    if legend_parts:
        legend = "\n".join(legend_parts)
        output_string += f"\n    {BOLD}Color legend:{RESET}\n{legend}\n"

    # Get all unique residue indices for PAE columns
    # First, determine the residue field name from oe_list
    residue_field_name = None
    if len(oe_list) > 0:
        if "residue_index" in oe_list[0]:
            residue_field_name = "residue_index"
        elif "aa_index" in oe_list[0]:
            residue_field_name = "aa_index"

    # Get all residue indices from pae_lookup (all residues that have PAE data)
    all_residue_indices = []
    if has_pae:
        # Get residue indices from pae_lookup keys (all residues with PAE data)
        all_residue_indices = sorted(pae_lookup.keys())
        # Also include any residue indices from oe_list that might not be in pae_lookup
        if residue_field_name:
            oe_residue_indices = set([entry[residue_field_name] for entry in oe_list])
            all_residue_indices = sorted(set(all_residue_indices) | oe_residue_indices)

    # Build header
    header = f"    {'Residue Index':<15} {'Observed':<12} {'Expected':<12} {'O/E':<10}"
    if has_plddt:
        header += f" {'pLDDT':<10}"
    if has_pae:
        header += " "  # Empty column before PAE
        for res_idx in all_residue_indices:
            header += f"{res_idx:>8}"

    # Add PAE header row if PAE columns are present
    if has_pae:
        # Calculate the width of columns before PAE section
        pre_pae_width = 15 + 12 + 12 + 10  # Residue Index + Observed + Expected + O/E
        if has_plddt:
            pre_pae_width += 10  # pLDDT
        pre_pae_width += 4  # "    " prefix
        pre_pae_width += 1  # Empty column before PAE

        # Create PAE header row with "PAE to residue:" label
        pae_header = (
            "    " + " " * (15 + 1) + " " * (12 + 1) + " " * (12 + 1) + " " * (10 + 1)
        )
        if has_plddt:
            pae_header += " " * (10 + 1)
        pae_header += " "  # Empty column before PAE
        # Add "PAE to residue:" label, centered over the PAE columns
        pae_col_width = 8 * len(all_residue_indices)
        pae_label = "PAE to residue:"
        padding = max(0, (pae_col_width - len(pae_label)) // 2)
        pae_header += " " * padding + pae_label
        output_string += pae_header + "\n"

    output_string += header + "\n"

    # Calculate separator length
    sep_len = 50
    if has_plddt:
        sep_len += 11
    if has_pae:
        sep_len += 1 + 8 * len(all_residue_indices)
    output_string += "    " + "-" * sep_len + "\n"

    residue_field_name = None
    if len(oe_list) > 0:
        if "residue_index" in oe_list[0]:
            residue_field_name = "residue_index"
        elif "aa_index" in oe_list[0]:
            residue_field_name = "aa_index"
        else:
            oe_list = [x.annotate(residue_index=i) for i, x in enumerate(oe_list)]
            residue_field_name = "residue_index"

    for oe_entry in oe_list:
        residue_idx = oe_entry[residue_field_name]
        obs = oe_entry.obs
        exp = oe_entry.exp
        oe_ratio = obs / exp if exp > 0 else 0.0

        row_str = f"    {residue_idx:<15} {obs:<12} {exp:<12.2f} {oe_ratio:<10.3f}"

        if has_plddt:
            plddt_val = plddt_lookup.get(residue_idx, None)
            if plddt_val is not None:
                # Color low pLDDT (below min_plddt) red
                if plddt_val < min_plddt:
                    plddt_str = f"{RED}{plddt_val:<10.1f}{RESET}"
                else:
                    plddt_str = f"{plddt_val:<10.1f}"
            else:
                plddt_str = "NA        "
            row_str += f" {plddt_str}"

        if has_pae:
            row_str += " "  # Empty column before PAE
            # Get PAE array for this residue
            pae_array = pae_lookup.get(residue_idx, None)
            for target_res_idx in all_residue_indices:
                pae_val = None
                if pae_array is not None:
                    try:
                        if target_res_idx < len(pae_array):
                            pae_val = pae_array[target_res_idx]
                    except (TypeError, IndexError):
                        pass

                if pae_val is not None:
                    if pae_val > max_pae:
                        pae_str = f"{RED}{pae_val:>8.2f}{RESET}"
                    else:
                        pae_str = f"{pae_val:>8.2f}"
                else:
                    pae_str = "     NA"
                row_str += pae_str

        output_string += row_str + "\n"

    return output_string


def debug_print_oe_table(
    oe_expr: hl.expr.ArrayExpression,
    min_plddt: Optional[float] = None,
    max_pae: Optional[float] = None,
    dist_mat_expr: Optional[hl.expr.ArrayExpression] = None,
) -> None:
    output_string = "\n\n\n"

    # Get first key for debugging.
    uniprot_id = None
    transcript_id = None
    ht = oe_expr._indices.source
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.enst.collect()[0]

    # Filter to the specific uniprot/transcript
    _ht_debug = ht.filter((ht.uniprot_id == uniprot_id) & (ht.enst == transcript_id))

    output_string += (
        f"{BOLD}=== Per-residue observed and expected values ==={RESET}\n\n"
    )

    # Collect data - also get plddt_lookup, pae_lookup, and dist_mat if available
    select_fields = ["oe"]
    # Check if fields exist in the row schema
    row_schema = _ht_debug.row.dtype
    has_plddt_lookup = "plddt_lookup" in row_schema.fields
    has_pae_lookup = "pae_lookup" in row_schema.fields
    has_dist_mat = "dist_mat" in row_schema.fields
    if has_plddt_lookup:
        select_fields.append("plddt_lookup")
    if has_pae_lookup:
        select_fields.append("pae_lookup")
    if has_dist_mat and max_pae is not None:
        # Also get dist_mat to extract PAE arrays if pae_lookup is not available
        select_fields.append("dist_mat")

    debug_data = _ht_debug.select(*select_fields).collect()

    if not debug_data:
        return

    row = debug_data[0]

    # Build pLDDT lookup if available
    plddt_lookup = None
    if (
        has_plddt_lookup
        and hasattr(row, "plddt_lookup")
        and row.plddt_lookup is not None
    ):
        # Create dictionary mapping residue_index to pLDDT value
        plddt_lookup = {
            entry.residue_index: entry.plddt
            for entry in row.plddt_lookup
            if entry is not None
        }

    # Build PAE lookup if available
    # PAE lookup structure: array of structs with {residue_index: int,
    # pae_array: array<float>}
    pae_lookup = None
    if max_pae is not None:
        # First try to get it from dist_mat_expr parameter (most reliable since
        # it's passed directly)
        if dist_mat_expr is not None:
            # Collect dist_mat_expr to extract PAE arrays
            _ht_dist = dist_mat_expr._indices.source
            _ht_dist = _ht_dist.annotate(dist_mat=dist_mat_expr)
            _ht_dist = _ht_dist.filter(
                (_ht_dist.uniprot_id == uniprot_id) & (_ht_dist.enst == transcript_id)
            )
            dist_data = _ht_dist.select("dist_mat").head(1).collect()
            if dist_data and len(dist_data) > 0:
                dist_row = dist_data[0]
                if (
                    hasattr(dist_row, "dist_mat")
                    and dist_row.dist_mat is not None
                    and len(dist_row.dist_mat) > 0
                ):
                    pae_lookup = {}
                    for entry in dist_row.dist_mat:
                        if (
                            entry is not None
                            and hasattr(entry, "residue_index")
                            and hasattr(entry, "pae_array")
                            and entry.pae_array is not None
                        ):
                            pae_lookup[entry.residue_index] = entry.pae_array
        # If not available from dist_mat_expr, try pae_lookup field in table
        if (
            (pae_lookup is None or len(pae_lookup) == 0)
            and has_pae_lookup
            and hasattr(row, "pae_lookup")
            and row.pae_lookup is not None
        ):
            # Create dictionary mapping residue_index to PAE array
            pae_lookup = {}
            for entry in row.pae_lookup:
                if (
                    entry is not None
                    and hasattr(entry, "residue_index")
                    and hasattr(entry, "pae_array")
                ):
                    pae_lookup[entry.residue_index] = entry.pae_array
        # If still not available, try to get it from dist_mat if available
        if (
            (pae_lookup is None or len(pae_lookup) == 0)
            and has_dist_mat
            and hasattr(row, "dist_mat")
            and row.dist_mat is not None
        ):
            # Extract PAE arrays from dist_mat entries
            pae_lookup = {}
            for entry in row.dist_mat:
                if (
                    entry is not None
                    and hasattr(entry, "residue_index")
                    and hasattr(entry, "pae_array")
                    and entry.pae_array is not None
                ):
                    pae_lookup[entry.residue_index] = entry.pae_array

    # Print OE array (per-residue observed/expected)
    output_string += get_debug_oe_table_string(
        row.oe,
        plddt_lookup=plddt_lookup,
        min_plddt=min_plddt,
        pae_lookup=pae_lookup,
        max_pae=max_pae,
    )

    logger.info(output_string)


def debug_print_oe_and_regions(
    ht: hl.Table,
    title: str = "OE and Regions",
) -> None:
    """
    Print debug output showing OE values and min_oe_upper regions for a given UniProt/transcript.

    :param ht: Hail Table with schema containing 'oe' and 'min_oe_upper' arrays
    :param uniprot_id: UniProt ID to filter to
    :param transcript_id: Transcript ID to filter to
    :param title: Title for the debug output
    :param pae_cutoff: PAE threshold - values above this are colored red, below are uncolored. Default is 15.0
    """
    output_string = "\n\n\n"

    # Get first key for debugging.
    uniprot_id = None
    transcript_id = None
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.transcript_id.collect()[0]

    # Filter to the specific uniprot/transcript
    _ht_debug = ht.filter(
        (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return

    output_string += f"{BOLD}=== {title} ==={RESET}\nUniProt ID: {uniprot_id}, Transcript ID: {transcript_id}\n"

    # Collect data
    debug_data = _ht_debug.select("oe", "min_oe_upper").collect()

    if not debug_data:
        return

    row = debug_data[0]

    # Print min_oe_upper array (best region for each center residue)
    output_string += f"\n    {BOLD}=== Min OE Upper Array (best region for each center residue) ==={RESET}\n\n"

    if not row.min_oe_upper or len(row.min_oe_upper) == 0:
        output_string += "No min_oe_upper entries found.\n"
        logger.info(output_string)
        return

    # Find the entry with the minimum OE value (not oe_upper)
    min_oe_idx = None
    min_oe_val = None
    for idx, entry in enumerate(row.min_oe_upper):
        if hasattr(entry, "oe") and entry.oe is not None:
            if min_oe_val is None or entry.oe < min_oe_val:
                min_oe_val = entry.oe
                min_oe_idx = idx

    for idx, entry in enumerate(row.min_oe_upper):
        is_min = min_oe_idx is not None and idx == min_oe_idx

        output_string += f"\n    {BOLD}Center Residue: {entry.residue_index}{RESET}\n"
        output_string += f"      Observed: {entry.obs:7d}\n"
        output_string += f"      Expected: {entry.exp:7.2f}\n"

        # Highlight minimum OE
        if is_min:
            output_string += (
                f"      {HIGHLIGHT}O/E:      {entry.oe:7.2f} (MINIMUM){RESET}\n"
            )
        else:
            output_string += f"      O/E:      {entry.oe:7.2f}\n"

        output_string += f"      OE Upper: {entry.oe_upper:7.2f}\n"

        # Region (array of residue indices)
        if hasattr(entry, "region") and entry.region:
            region_str = ", ".join([f"{r:7d}" for r in entry.region])
            output_string += (
                f"      Region residues ({len(entry.region)}): {region_str}\n"
            )
        else:
            output_string += f"      Region: None or empty\n"

        # Distances (array of distances)
        if hasattr(entry, "dists") and entry.dists:
            dists_str = ", ".join([f"{d:7.2f}" for d in entry.dists])
            output_string += (
                f"      Distances ({len(entry.dists)}):       {dists_str}\n"
            )
        else:
            output_string += f"      Distances: None or empty\n"

    logger.info(output_string)


def _debug_calculate_min_max_for_coloring(
    dist_mat_expr: hl.expr.ArrayExpression,
) -> Optional[Dict[str, int]]:
    """
    Calculate min/max values from dist_mat_expr for consistent coloring across debug statements.

    This ensures that residue indices and distances use the same color scale across
    all debug output, making it easier to compare values visually.

    :param dist_mat_expr: Hail expression for distance matrix array
    :return: Dictionary with min_res_idx, max_res_idx, min_dist, max_dist, or None if calculation fails
    """
    try:
        # Get min/max values from the original dist_mat_expr for consistent coloring
        # We'll calculate these from the first row's dist_mat
        _ht_debug = dist_mat_expr._indices.source
        _ht_debug = _ht_debug.annotate(dist_mat=dist_mat_expr).head(1)
        if _ht_debug.count() > 0:
            debug_data = _ht_debug.select("dist_mat").collect()
            if debug_data and len(debug_data[0].dist_mat) > 0:
                all_residue_indices = [
                    entry.residue_index for entry in debug_data[0].dist_mat
                ]
                all_distances = [entry.dist for entry in debug_data[0].dist_mat]
                if all_residue_indices and all_distances:
                    return {
                        "min_res_idx": min(all_residue_indices),
                        "max_res_idx": max(all_residue_indices),
                        "min_dist": min(all_distances),
                        "max_dist": max(all_distances),
                    }
    except Exception:
        # If calculation fails, return None (debug functions will handle None
        # gracefully)
        pass
    return None


def debug_print_dist_mat_with_colors(
    dist_mat_expr: hl.expr.ArrayExpression,
    title: str = "Distance matrix with PAE",
    max_pae: Optional[float] = None,
    min_plddt: Optional[float] = None,
    oe_expr: Optional[hl.expr.ArrayExpression] = None,
    min_res_idx: Optional[int] = None,
    max_res_idx: Optional[int] = None,
    min_dist: Optional[float] = None,
    max_dist: Optional[float] = None,
) -> Dict[int, str]:
    """
    Generate colored debug output showing distance matrix with PAE values color-coded.
    Returns a dictionary where keys are center residue indices and values are output strings.

    :param dist_mat_expr: Array expression with distance matrix entries (must have dist, residue_index, and optionally pae, plddt)
    :param title: Title for the debug output
    :param max_pae: Maximum allowed PAE value - values above this are colored red, below are uncolored. If None, PAE won't be shown.
    :param min_plddt: Minimum allowed pLDDT value. If None, pLDDT won't be shown.
    :param oe_expr: Optional array expression with OE data (obs, exp, oe, oe_upper). If provided, will print
        observed, expected, cumulative observed, cumulative expected, cumulative O/E, and cumulative OE upper.
    :return: Dictionary mapping center residue indices to their debug output strings.
        For the "Before sorting" full matrix case, returns a dict with key -1 and the full matrix output.
    """

    # Get first key for debugging.
    ht = dist_mat_expr._indices.source
    uniprot_id = None
    transcript_id = None
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.enst.collect()[0]

    output_string = ""

    # For "Before sorting by nearest neighbor", if PAE and pLDDT are not shown,
    # show the full distance matrix (all residues as rows, distances to all
    # residues as columns)
    show_full_matrix = (
        "Before sorting by nearest neighbor" in title
        and max_pae is None
        and min_plddt is None
        and oe_expr is None
    )

    # Filter and collect data
    select_fields = {"dist_mat": dist_mat_expr}
    if oe_expr is not None:
        select_fields["oe_expr"] = oe_expr

    _ht_debug = ht.annotate(**select_fields)
    _ht_debug = _ht_debug.filter(
        (_ht_debug.uniprot_id == uniprot_id) & (_ht_debug.enst == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return

    header_string = ""
    if not show_full_matrix:
        header_string += "    "

    header_string += f"{BOLD}=== {title} ==={RESET}\n"

    # Collect data for colored output
    select_fields = ["dist_mat"]
    if oe_expr is not None:
        select_fields.append("oe_expr")

    debug_data = _ht_debug.select(*select_fields).collect()

    if not debug_data:
        return {}

    # Dictionary to store output per center residue
    output_dict = {}

    # Check if PAE and pLDDT fields exist in the struct type (only if requested)
    has_pae_field = False
    has_plddt_field = False
    if debug_data and len(debug_data[0].dist_mat) > 0:
        first_entry = debug_data[0].dist_mat[0]
        if max_pae is not None:
            has_pae_field = hasattr(first_entry, "pae")
        if min_plddt is not None:
            has_plddt_field = hasattr(first_entry, "plddt")

    # Check if OE fields exist
    has_oe_data = oe_expr is not None and debug_data and len(debug_data[0].oe_expr) > 0

    # ANSI color codes for distances (blue to white gradient)
    BLUE = "\033[94m"  # Closest (bright blue)
    CYAN = "\033[96m"  # Medium-close (cyan)
    WHITE = "\033[97m"  # Furthest (white)

    # ANSI color codes for residue indices (gradient from green to magenta)
    GREEN_IDX = "\033[92m"  # First in sorted array (closest)
    YELLOW_IDX = "\033[93m"  # Middle
    MAGENTA_IDX = "\033[95m"  # Last in sorted array (furthest)

    # ANSI color code for highlighting minimum OE upper
    HIGHLIGHT = "\033[1m\033[4m\033[92m"  # Bold and underline for minimum OE upper

    # Add color legends
    legend_parts = [
        f"            {BOLD}Residue index colors:{RESET} {GREEN_IDX}Lowest index{RESET} → {YELLOW_IDX}Middle{RESET} → {MAGENTA_IDX}Highest index{RESET}",
        f"            {BOLD}Distance colors:{RESET} {BLUE}Closest{RESET} → {CYAN}Medium{RESET} → {WHITE}Furthest{RESET}",
    ]
    if max_pae is not None and has_pae_field:
        legend_parts.append(
            f"            {BOLD}PAE colors:{RESET} Uncolored (≤{max_pae})  {RED}Red (>{max_pae}){RESET}"
        )
    if min_plddt is not None and has_plddt_field:
        legend_parts.append(
            f"            {BOLD}pLDDT colors:{RESET} Uncolored (≥{min_plddt})  {RED}Red (<{min_plddt}){RESET}"
        )
    legend = "\n".join(legend_parts)
    legend_string = f"\n        {BOLD}Color legend:{RESET}\n{legend}\n\n"

    # If showing full matrix, we need all rows to build it, but only show it once
    if show_full_matrix:
        # Build full distance matrix from all center residues
        # First, get all unique residue indices from the first row
        if not debug_data or len(debug_data) == 0:
            return {}

        first_row = debug_data[0]
        all_unique_residues = sorted(
            [entry.residue_index for entry in first_row.dist_mat]
        )

        # Build distance lookup: center_residue -> {target_residue: distance}
        dist_lookup = {}
        for center_row in debug_data:
            center_res = center_row.aa_index
            if center_res not in dist_lookup:
                dist_lookup[center_res] = {}
            for entry in center_row.dist_mat:
                target_res = entry.residue_index
                dist_lookup[center_res][target_res] = entry.dist

        # Calculate min/max distances for coloring
        all_distances = []
        for center_res in dist_lookup:
            for target_res in dist_lookup[center_res]:
                all_distances.append(dist_lookup[center_res][target_res])

        if all_distances:
            row_min_dist = min(all_distances) if min_dist is None else min_dist
            row_max_dist = max(all_distances) if max_dist is None else max_dist
            dist_range = (
                row_max_dist - row_min_dist if row_max_dist > row_min_dist else 1.0
            )
        else:
            row_min_dist = 0.0
            row_max_dist = 1.0
            dist_range = 1.0

        # Calculate residue index range for coloring
        if all_unique_residues:
            row_min_res_idx = (
                min(all_unique_residues) if min_res_idx is None else min_res_idx
            )
            row_max_res_idx = (
                max(all_unique_residues) if max_res_idx is None else max_res_idx
            )
            res_idx_range = (
                row_max_res_idx - row_min_res_idx
                if row_max_res_idx > row_min_res_idx
                else 1.0
            )
        else:
            row_min_res_idx = 0
            row_max_res_idx = 1
            res_idx_range = 1.0

        # Helper function to get color for residue index
        def get_residue_color(res_idx):
            if res_idx_range > 0:
                normalized_res_idx = (res_idx - row_min_res_idx) / res_idx_range
                if normalized_res_idx < 0.33:
                    return GREEN_IDX
                elif normalized_res_idx < 0.67:
                    return YELLOW_IDX
                else:
                    return MAGENTA_IDX
            else:
                return GREEN_IDX

        # Create matrix header with colored residue indices
        header = f"\n        {'Residue':<10}"
        for res_idx in all_unique_residues:
            color = get_residue_color(res_idx)
            header += f"{color}{res_idx:>8}{RESET}"
        output_string += header + "\n"
        output_string += "        " + "-" * (10 + 8 * len(all_unique_residues)) + "\n"

        # Create matrix rows with colored row headers
        for from_res in all_unique_residues:
            if from_res in dist_lookup:
                row_color = get_residue_color(from_res)
                row_str = f"        {row_color}{from_res:<10}{RESET}"
                for to_res in all_unique_residues:
                    if to_res in dist_lookup[from_res]:
                        dist_val = dist_lookup[from_res][to_res]
                        # Color code distance
                        if dist_range > 0:
                            normalized = (dist_val - row_min_dist) / dist_range
                            if normalized < 0.5:
                                color = BLUE if normalized < 0.25 else CYAN
                            else:
                                color = CYAN if normalized < 0.75 else WHITE
                        else:
                            color = BLUE
                        row_str += f"{color}{dist_val:>8.2f}{RESET}"
                    else:
                        row_str += f"{'?':>8}"
                output_string += row_str + "\n"

        # For full matrix, store with key -1 (special key for full matrix)
        full_output = "\n\n" + header_string + legend_string + output_string
        output_dict[-1] = full_output
        return output_dict

    # Normal processing: show data for each center residue
    for row in debug_data:
        center_idx = row.aa_index

        # Use provided min/max values if available, otherwise calculate from current row
        residue_indices = [entry.residue_index for entry in row.dist_mat]
        if not residue_indices:
            # Empty dist_mat - skip this row
            continue

        if min_res_idx is None or max_res_idx is None:
            # Calculate from current row if not provided
            row_min_res_idx = min(residue_indices)
            row_max_res_idx = max(residue_indices)
        else:
            # Use provided values
            row_min_res_idx = min_res_idx
            row_max_res_idx = max_res_idx
        res_idx_range = (
            row_max_res_idx - row_min_res_idx
            if row_max_res_idx > row_min_res_idx
            else 1
        )  # Avoid division by zero

        # Find min and max distances for gradient
        distances = [entry.dist for entry in row.dist_mat]
        if not distances:
            # Empty dist_mat - skip this row
            continue

        if min_dist is None or max_dist is None:
            # Calculate from current row if not provided
            row_min_dist = min(distances)
            row_max_dist = max(distances)
        else:
            # Use provided values
            row_min_dist = min_dist
            row_max_dist = max_dist
        dist_range = (
            row_max_dist - row_min_dist if row_max_dist > row_min_dist else 1.0
        )  # Avoid division by zero

        # Get OE data if available (should be parallel to dist_mat)
        oe_entries = row.oe_expr if has_oe_data else None

        # Find minimum OE upper value and its index (if OE data is available)
        min_oe_upper_idx = None
        min_oe_upper_val = None
        if has_oe_data and oe_entries:
            # Filter to entries that have oe_upper defined and not None
            valid_oe_entries = [
                (idx, entry)
                for idx, entry in enumerate(oe_entries)
                if hasattr(entry, "oe_upper") and entry.oe_upper is not None
            ]
            if valid_oe_entries:
                # Find the entry with minimum OE upper
                min_entry = min(valid_oe_entries, key=lambda x: x[1].oe_upper)
                min_oe_upper_idx = min_entry[0]
                min_oe_upper_val = min_entry[1].oe_upper

        # Build colored output
        residue_strs = []
        dist_strs = []
        pae_strs = []
        plddt_strs = []
        obs_strs = []
        exp_strs = []
        cum_obs_strs = []
        cum_exp_strs = []
        cum_oe_strs = []
        cum_oe_upper_strs = []

        for idx, entry in enumerate(row.dist_mat):
            # Color code residue indices by their actual index value (not position in
            # array)
            res_idx_val = entry.residue_index
            if res_idx_range > 0:
                # Normalize residue index to 0-1 range
                normalized_res_idx = (res_idx_val - row_min_res_idx) / res_idx_range
                # Interpolate between green (lowest index) and magenta (highest index)
                if normalized_res_idx < 0.33:
                    color_idx = GREEN_IDX
                elif normalized_res_idx < 0.67:
                    color_idx = YELLOW_IDX
                else:
                    color_idx = MAGENTA_IDX
            else:
                # All residue indices are the same
                color_idx = GREEN_IDX
            residue_strs.append(f"{color_idx}{res_idx_val:7d}{RESET}")

            # Color code distances (blue for closest, white for furthest)
            dist_val = entry.dist
            if dist_range > 0:
                # Normalize distance to 0-1 range
                normalized = (dist_val - row_min_dist) / dist_range
                # Interpolate between blue (0.0) and white (1.0)
                if normalized < 0.5:
                    # Blue to cyan gradient
                    color = BLUE if normalized < 0.25 else CYAN
                else:
                    # Cyan to white gradient
                    color = CYAN if normalized < 0.75 else WHITE
            else:
                # All distances are the same
                color = BLUE
            dist_strs.append(f"{color}{dist_val:7.2f}{RESET}")

            # Color code PAE values if available and requested
            if max_pae is not None and has_pae_field:
                pae_val = entry.pae
                if pae_val > max_pae:
                    # Color red if above cutoff
                    pae_strs.append(f"{RED}{pae_val:7.2f}{RESET}")
                else:
                    # No color if below cutoff
                    pae_strs.append(f"{pae_val:7.2f}")
            elif max_pae is not None:
                pae_strs.append("     NA")

            # Add pLDDT values if available and requested
            if min_plddt is not None and has_plddt_field:
                plddt_val = entry.plddt
                if plddt_val is not None:
                    # Color low pLDDT (below min_plddt) red
                    if plddt_val < min_plddt:
                        plddt_strs.append(f"{RED}{plddt_val:7.1f}{RESET}")
                    else:
                        plddt_strs.append(f"{plddt_val:7.1f}")
                else:
                    plddt_strs.append("     NA")
            elif min_plddt is not None:
                plddt_strs.append("     NA")

            # Add OE data if available (no coloring)
            # Check if this residue is excluded from stats (for
            # exclude_low_plddt_from_stats)
            is_excluded = False
            if hasattr(entry, "exclude_from_stats"):
                is_excluded = entry.exclude_from_stats is True

            if has_oe_data and oe_entries and idx < len(oe_entries):
                oe_entry = oe_entries[idx]

                # If excluded from stats, show NA for all OE values
                if is_excluded:
                    obs_strs.append("     NA")
                    exp_strs.append("     NA")
                    cum_obs_strs.append("     NA")
                    cum_exp_strs.append("     NA")
                    cum_oe_strs.append("     NA")
                    cum_oe_upper_strs.append("     NA")
                else:
                    # After get_cumulative_oe, obs and exp are cumulative
                    # Calculate per-residue values from cumulative
                    if idx == 0:
                        # First entry: per-residue = cumulative
                        per_residue_obs = oe_entry.obs
                        per_residue_exp = oe_entry.exp
                    else:
                        # Subsequent entries: per-residue = current cumulative -
                        # previous cumulative
                        prev_entry = oe_entries[idx - 1]
                        per_residue_obs = (
                            oe_entry.obs - prev_entry.obs
                            if (oe_entry.obs is not None and prev_entry.obs is not None)
                            else None
                        )
                        per_residue_exp = (
                            oe_entry.exp - prev_entry.exp
                            if (oe_entry.exp is not None and prev_entry.exp is not None)
                            else None
                        )

                    # Handle None values in formatting
                    if per_residue_obs is not None:
                        obs_strs.append(f"{per_residue_obs:7d}")
                    else:
                        obs_strs.append("     NA")

                    if per_residue_exp is not None:
                        exp_strs.append(f"{per_residue_exp:7.2f}")
                    else:
                        exp_strs.append("     NA")

                    # Cumulative values
                    if oe_entry.obs is not None:
                        cum_obs_strs.append(f"{oe_entry.obs:7d}")
                    else:
                        cum_obs_strs.append("     NA")

                    if oe_entry.exp is not None:
                        cum_exp_strs.append(f"{oe_entry.exp:7.2f}")
                    else:
                        cum_exp_strs.append("     NA")

                    if oe_entry.oe is not None:
                        cum_oe_strs.append(f"{oe_entry.oe:7.2f}")
                    else:
                        cum_oe_strs.append("     NA")

                    # Highlight minimum OE upper value
                    if oe_entry.oe_upper is not None:
                        if min_oe_upper_idx is not None and idx == min_oe_upper_idx:
                            cum_oe_upper_strs.append(
                                f"{HIGHLIGHT}{oe_entry.oe_upper:7.2f}{RESET}"
                            )
                        else:
                            cum_oe_upper_strs.append(f"{oe_entry.oe_upper:7.2f}")
                    else:
                        cum_oe_upper_strs.append("     NA")
            elif has_oe_data:
                # OE data exists but this index is out of range
                obs_strs.append("     NA")
                exp_strs.append("     NA")
                cum_obs_strs.append("     NA")
                cum_exp_strs.append("     NA")
                cum_oe_strs.append("     NA")
                cum_oe_upper_strs.append("     NA")

        # Build output for this center residue (content only, no header/legend)
        output = (
            f"        Residue indices:     {', '.join(residue_strs)}\n"
            f"        Distances:           {', '.join(dist_strs)}\n"
        )

        if max_pae is not None and has_pae_field:
            output += f"        PAE:                 {', '.join(pae_strs)}\n"

        if min_plddt is not None and has_plddt_field:
            output += f"        pLDDT:               {', '.join(plddt_strs)}\n"

        if has_oe_data:
            output += (
                f"        Observed:            {', '.join(obs_strs)}\n"
                f"        Expected:            {', '.join(exp_strs)}\n"
                f"        Cumulative Observed: {', '.join(cum_obs_strs)}\n"
                f"        Cumulative Expected: {', '.join(cum_exp_strs)}\n"
                f"        Cumulative O/E:      {', '.join(cum_oe_strs)}\n"
                f"        Cumulative OE Upper: {', '.join(cum_oe_upper_strs)}\n"
            )
            if min_oe_upper_idx is not None:
                min_oe_entry = oe_entries[min_oe_upper_idx]
                output += f"\n        {BOLD}Minimum OE Upper:{RESET} {HIGHLIGHT}{min_oe_upper_val:.2f}{RESET} at index {min_oe_upper_idx} (residue {row.dist_mat[min_oe_upper_idx].residue_index})\n\n"

                # Add Observed, Expected, O/E, and OE Upper values
                if hasattr(min_oe_entry, "obs") and min_oe_entry.obs is not None:
                    output += f"            Region observed:    {min_oe_entry.obs:7d}\n"
                else:
                    output += f"            Region observed:    {'NA':>7}\n"

                if hasattr(min_oe_entry, "exp") and min_oe_entry.exp is not None:
                    output += (
                        f"            Region expected:    {min_oe_entry.exp:7.2f}\n"
                    )
                else:
                    output += f"            Region expected:    {'NA':>7}\n"

                if hasattr(min_oe_entry, "oe") and min_oe_entry.oe is not None:
                    output += (
                        f"            Region O/E:         {min_oe_entry.oe:7.2f}\n"
                    )
                else:
                    output += f"            Region O/E:         {'NA':>7}\n"

                if (
                    hasattr(min_oe_entry, "oe_upper")
                    and min_oe_entry.oe_upper is not None
                ):
                    output += f"            Region OE Upper:    {min_oe_entry.oe_upper:7.2f}\n"
                else:
                    output += f"            Region OE Upper:    {'NA':>7}\n"

                # Add region information if available
                if (
                    hasattr(min_oe_entry, "region")
                    and min_oe_entry.region is not None
                    and len(min_oe_entry.region) > 0
                ):
                    region_str = ", ".join([f"{r}" for r in min_oe_entry.region])
                    output += f"            Region residues ({len(min_oe_entry.region)}): {region_str}\n"
                else:
                    # If region is not available, build it from dist_mat up to
                    # min_oe_upper_idx
                    region_residues = [
                        entry.residue_index
                        for entry in row.dist_mat[: min_oe_upper_idx + 1]
                    ]
                    region_str = ", ".join([f"{r}" for r in region_residues])
                    output += f"            Region residues ({len(region_residues)}): {region_str}\n"

        # Store output for this center residue
        # Include header and legend in the output string
        full_output = "\n\n" + header_string + legend_string + output
        output_dict[center_idx] = full_output

    return output_dict


def _debug_get_3d_residue(
    dist_mat_expr: hl.expr.ArrayExpression,
    oe_expr: hl.expr.ArrayExpression,
    max_pae: Optional[float],
    min_plddt: Optional[float],
    pae_cutoff_method: str,
    plddt_cutoff_method: Optional[str],
    center_residue_index_expr: hl.expr.Int32Expression,
    debug_min_max: Optional[dict],
    debug_outputs_by_residue: dict,
    step_name: str,
    oe_expr_after_calc: Optional[hl.expr.ArrayExpression] = None,
):
    """
    Handle all debug output for get_3d_residue.

    :param dist_mat_expr: Distance matrix expression at current step.
    :param oe_expr: Original OE expression.
    :param max_pae: Maximum PAE cutoff.
    :param min_plddt: Minimum pLDDT cutoff.
    :param pae_cutoff_method: PAE cutoff method.
    :param plddt_cutoff_method: pLDDT cutoff method.
    :param center_residue_index_expr: Center residue index expression.
    :param debug_min_max: Dict with min/max values for consistent coloring.
    :param debug_outputs_by_residue: Dict to store debug outputs.
    :param step_name: Name of the current step.
    :param oe_expr_after_calc: OE expression after calculate_oe_upper (optional).
    """
    if step_name == "initial":
        # Extract UniProt ID and Transcript ID once for debug output
        _ht_debug = dist_mat_expr._indices.source
        first_row = _ht_debug.head(1)
        if first_row.count() > 0:
            debug_uniprot_id = first_row.uniprot_id.collect()[0]
            debug_transcript_id = first_row.enst.collect()[0]
            # Print UniProt ID and Transcript ID once at the beginning
            logger.info(
                f"\nUniProt ID: {debug_uniprot_id}, Transcript ID: {debug_transcript_id}\n"
            )
            debug_print_oe_table(
                oe_expr,
                min_plddt=min_plddt,
                max_pae=max_pae,
                dist_mat_expr=dist_mat_expr,
            )

    elif step_name == "before_sorting":
        debug_dict = debug_print_dist_mat_with_colors(
            dist_mat_expr,
            title="get_3d_residue: Before sorting by nearest neighbor",
            min_res_idx=debug_min_max["min_res_idx"] if debug_min_max else None,
            max_res_idx=debug_min_max["max_res_idx"] if debug_min_max else None,
            min_dist=debug_min_max["min_dist"] if debug_min_max else None,
            max_dist=debug_min_max["max_dist"] if debug_min_max else None,
        )
        # Store debug output (key -1 is for full matrix, otherwise it's center_res_idx)
        for center_idx, output_str in debug_dict.items():
            if center_idx not in debug_outputs_by_residue:
                debug_outputs_by_residue[center_idx] = {}
            debug_outputs_by_residue[center_idx]["Before sorting"] = output_str

    elif step_name == "after_sorting":
        debug_dict = debug_print_dist_mat_with_colors(
            dist_mat_expr,
            title="get_3d_residue: After sorting by nearest neighbor",
            max_pae=max_pae,
            min_plddt=min_plddt,
            min_res_idx=debug_min_max["min_res_idx"] if debug_min_max else None,
            max_res_idx=debug_min_max["max_res_idx"] if debug_min_max else None,
            min_dist=debug_min_max["min_dist"] if debug_min_max else None,
            max_dist=debug_min_max["max_dist"] if debug_min_max else None,
        )
        # Store debug output
        for center_idx, output_str in debug_dict.items():
            if center_idx not in debug_outputs_by_residue:
                debug_outputs_by_residue[center_idx] = {}
            debug_outputs_by_residue[center_idx]["After sorting"] = output_str

    elif step_name == "pae_matrix_before_filter":
        debug_dict = debug_print_pae_matrix_for_region(
            dist_mat_expr,
            center_residue_index_expr,
            max_pae=max_pae,
            title="get_3d_residue: PAE matrix before filter_on_pairwise_pae_in_region (after pLDDT filtering)",
            min_plddt=min_plddt,
            plddt_cutoff_method=plddt_cutoff_method,
        )
        # Store debug output
        for center_idx, output_str in debug_dict.items():
            if center_idx not in debug_outputs_by_residue:
                debug_outputs_by_residue[center_idx] = {}
            debug_outputs_by_residue[center_idx][
                "PAE matrix before filter"
            ] = output_str

    elif step_name == "after_calculate_oe_upper":
        debug_dict = debug_print_dist_mat_with_colors(
            dist_mat_expr,
            title="get_3d_residue: After calculate_oe_upper",
            max_pae=max_pae,
            min_plddt=min_plddt,
            oe_expr=oe_expr_after_calc,
            min_res_idx=debug_min_max["min_res_idx"] if debug_min_max else None,
            max_res_idx=debug_min_max["max_res_idx"] if debug_min_max else None,
            min_dist=debug_min_max["min_dist"] if debug_min_max else None,
            max_dist=debug_min_max["max_dist"] if debug_min_max else None,
        )
        # Store debug output
        for center_idx, output_str in debug_dict.items():
            if center_idx not in debug_outputs_by_residue:
                debug_outputs_by_residue[center_idx] = {}
            debug_outputs_by_residue[center_idx][
                "After calculate_oe_upper"
            ] = output_str

    elif step_name == "final":
        # Print all debug outputs grouped by residue
        if debug_outputs_by_residue:
            debug_string = ""
            # First, handle the full matrix case (key -1) if it exists
            if -1 in debug_outputs_by_residue:
                debug_string += debug_outputs_by_residue[-1]["Before sorting"]
                del debug_outputs_by_residue[-1]

            # Then print per-residue outputs
            for center_idx in sorted(debug_outputs_by_residue.keys()):
                residue_outputs = debug_outputs_by_residue[center_idx]

                # Print center residue index once before all steps.
                combined_output = f"\n{BOLD}Center residue index: {center_idx}{RESET}\n"

                # Combine all steps for this residue.
                for step_name_inner in [
                    "Before sorting",
                    "After sorting",
                    "PAE matrix before filter",
                    "After calculate_oe_upper",
                ]:
                    if step_name_inner in residue_outputs:
                        combined_output += residue_outputs[step_name_inner] + "\n"
                if combined_output:
                    debug_string += combined_output

            logger.info(debug_string)


def debug_print_forward_round(
    ht: hl.Table,
    uniprot_id: str,
    transcript_id: str,
    model_comparison_method: str,
    title: str = "Forward Algorithm Round",
    oe_upper_method: str = "gamma",
    min_exp_mis: int = 1,
) -> str:
    """
    Print debug output for a round of the forward algorithm.

    :param ht: Hail Table with forward algorithm state
    :param uniprot_id: UniProt ID
    :param transcript_id: Transcript ID
    :param model_comparison_method: Model comparison method being used
    :param title: Title for the debug output
    """
    output_string = "\n"

    # Filter to the specific uniprot/transcript
    _ht_debug = ht.filter(
        (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
    )

    if _ht_debug.count() == 0:
        logger.info(f"{title}: No data found for {uniprot_id} / {transcript_id}")
        return

    output_string += f"{BOLD}    === {title} ==={RESET}\n"

    # Collect data
    # Note: key fields (idx, uniprot_id, transcript_id) are automatically included in collect()
    # Check if excluded_residues is available
    has_excluded_res_debug = "excluded_residues" in _ht_debug.row.dtype.fields
    debug_data = _ht_debug.select(
        "region",
        "selected",
        "selected_nll",
        "null_model",
        "best_aic",
        "found_best",
        "_region",
        "_updated_null",
        "oe",
        "center_residue_index",
        *["excluded_residues"] if has_excluded_res_debug else [],
    ).collect()

    if not debug_data:
        return ""

    # Print current state (use first row for state info, should be same for all)
    row = debug_data[0] if debug_data else None
    if row is None:
        return ""

    output_string += f"        {BOLD}Current State:{RESET}\n"
    output_string += (
        f"            Selected regions: {len(row.selected) if row.selected else 0}\n"
    )
    if row.selected_nll is None:
        selected_nll_str = "NA"
    elif row.selected_nll == 0:
        selected_nll_str = "NA (set to 0.0 as initial state)"
    else:
        selected_nll_str = f"{row.selected_nll:.4f}"
    output_string += f"            Selected NLL: {selected_nll_str}\n"
    best_aic_str = f"{row.best_aic:.4f}" if row.best_aic is not None else "NA"
    output_string += f"            Best AIC: {best_aic_str}\n"
    output_string += f"            Found best: {row.found_best}\n"

    if row.null_model:
        output_string += f"\n        {BOLD}Null Model (catch-all):{RESET}\n"
        obs_str = (
            f"{row.null_model.obs:7d}" if row.null_model.obs is not None else "     NA"
        )
        output_string += f"            Observed: {obs_str}\n"
        exp_str = (
            f"{row.null_model.exp:7.2f}"
            if row.null_model.exp is not None
            else "     NA"
        )
        output_string += f"            Expected: {exp_str}\n"
        # Color aggregate O/E values
        if row.null_model.oe is not None:
            oe_str = f"{row.null_model.oe:7.2f}"
        else:
            oe_str = "     NA"
        output_string += f"            O/E:      {oe_str}\n"
        # Color NLL values (more negative = better fit)
        if row.null_model.nll is not None:
            nll_str = f"{BOLD}{UNDERLINE}{row.null_model.nll:.4f}{RESET}"
        else:
            nll_str = "NA"
        output_string += f"            NLL:      {nll_str}\n\n"
        if hasattr(row.null_model, "region") and row.null_model.region:
            region_str = ", ".join([f"{r:7d}" for r in row.null_model.region])
            output_string += f"            Region residues: {region_str}\n"

            # Print per-residue observed and expected values for null model
            # Get oe from first row (should be same for all)
            if debug_data and len(debug_data) > 0:
                first_row = debug_data[0]
                if hasattr(first_row, "oe") and first_row.oe:
                    obs_list = []
                    exp_list = []
                    # Get excluded_residues if available
                    excluded_residues_dict = None
                    if (
                        hasattr(first_row, "excluded_residues")
                        and first_row.excluded_residues is not None
                    ):
                        excluded_residues_dict = first_row.excluded_residues
                    for res_idx in row.null_model.region:
                        # Check if this residue is excluded from stats
                        is_excluded = False
                        if excluded_residues_dict is not None:
                            if hasattr(excluded_residues_dict, "get"):
                                exclude_flag = excluded_residues_dict.get(res_idx)
                                is_excluded = (
                                    exclude_flag is True
                                    if exclude_flag is not None
                                    else False
                                )
                            elif isinstance(excluded_residues_dict, dict):
                                is_excluded = (
                                    excluded_residues_dict.get(res_idx, False) is True
                                )

                        if is_excluded:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                        elif res_idx < len(first_row.oe):
                            oe_entry = first_row.oe[res_idx]
                            obs_val = (
                                oe_entry.obs
                                if hasattr(oe_entry, "obs") and oe_entry.obs is not None
                                else None
                            )
                            exp_val = (
                                oe_entry.exp
                                if hasattr(oe_entry, "exp") and oe_entry.exp is not None
                                else None
                            )

                            obs_list.append(
                                f"{obs_val:7d}" if obs_val is not None else "    NA"
                            )
                            exp_list.append(
                                f"{exp_val:7.2f}" if exp_val is not None else "    NA"
                            )
                        else:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                    output_string += f"            Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                    output_string += f"            Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"

    # Print all candidate regions being evaluated
    # First, filter to only valid candidates (those with valid _region data
    # and exp >= min_exp_mis)
    valid_candidates = []
    for row in debug_data:
        if (
            hasattr(row, "_region")
            and row._region
            and hasattr(row._region, "obs")
            and row._region.obs is not None
            and hasattr(row._region, "exp")
            and row._region.exp is not None
            and row._region.exp >= min_exp_mis
            and hasattr(row._region, "nll")
            and row._region.nll is not None
        ):
            valid_candidates.append(row)

    # Find the best candidate using the same logic as the aggregation:
    # 1. Minimum NLL
    # 2. If NLL is equal, use OE upper (lower OE upper is better)
    best_candidate_idx = None
    best_nll = None
    best_oe_upper = None

    # Import OE upper functions for tie-breaking
    from gnomad_constraint.experimental.proemis3d.utils import (
        chisq_upper_ci,
        gamma_upper_ci,
    )

    oe_upper_func = gamma_upper_ci if oe_upper_method == "gamma" else chisq_upper_ci

    for row_idx, row in enumerate(valid_candidates):
        if (
            hasattr(row, "_region")
            and row._region
            and hasattr(row._region, "nll")
            and row._region.nll is not None
        ):
            # Calculate OE upper for tie-breaking
            oe_upper_val = None
            if (
                hasattr(row._region, "obs")
                and row._region.obs is not None
                and hasattr(row._region, "exp")
                and row._region.exp is not None
                and row._region.exp > 0
            ):
                try:
                    oe_upper_val = hl.eval(
                        oe_upper_func(
                            hl.literal(row._region.obs),
                            hl.literal(row._region.exp),
                            0.05,
                        )
                    )
                except Exception:
                    pass

            # Check if this is a better candidate
            is_better = False
            if best_nll is None:
                is_better = True
            elif row._region.nll < best_nll:
                is_better = True
            elif row._region.nll == best_nll and oe_upper_val is not None:
                # Tie-breaking: lower OE upper is better
                if best_oe_upper is None or oe_upper_val < best_oe_upper:
                    is_better = True

            if is_better:
                best_candidate_idx = row_idx
                best_nll = row._region.nll
                best_oe_upper = oe_upper_val

    output_string += (
        f"\n        {BOLD}Candidate Regions ({len(valid_candidates)} total):{RESET}\n"
    )
    for row_idx, row in enumerate(valid_candidates):
        if hasattr(row, "_region") and row._region:
            idx_str = str(row.idx) if row.idx is not None else "NA"
            center_res_str = (
                f" (center residue: {row.center_residue_index})"
                if hasattr(row, "center_residue_index")
                and row.center_residue_index is not None
                else ""
            )
            output_string += (
                f"\n            {BOLD}Candidate {row_idx + 1}{center_res_str}:{RESET}\n"
            )
            obs_str = (
                f"{row._region.obs:7d}"
                if (hasattr(row._region, "obs") and row._region.obs is not None)
                else "     NA"
            )
            output_string += f"                Observed: {obs_str}\n"
            exp_str = (
                f"{row._region.exp:7.2f}"
                if (hasattr(row._region, "exp") and row._region.exp is not None)
                else "     NA"
            )
            output_string += f"                Expected: {exp_str}\n"
            # Keep O/E uncolored for candidates (selection is based on NLL, not O/E)
            oe_str = (
                f"{row._region.oe:7.2f}"
                if (hasattr(row._region, "oe") and row._region.oe is not None)
                else "     NA"
            )
            output_string += f"                O/E:      {oe_str}\n"
            # Calculate and display OE upper
            if (
                hasattr(row._region, "obs")
                and row._region.obs is not None
                and hasattr(row._region, "exp")
                and row._region.exp is not None
            ):
                oe_upper_val = _calculate_oe_upper(
                    row._region.obs, row._region.exp, oe_upper_method
                )
                oe_upper_str = (
                    f"{oe_upper_val:7.2f}" if oe_upper_val is not None else "     NA"
                )
            else:
                oe_upper_str = "     NA"
            output_string += f"                OE Upper: {oe_upper_str}\n"
            # Color NLL values (more negative = better fit)
            if (
                hasattr(row._region, "region_nll")
                and row._region.region_nll is not None
            ):
                region_nll_str = f"{row._region.region_nll:.4f}"
            else:
                region_nll_str = "NA"
            output_string += f"                Region NLL: {region_nll_str}\n"
            # Highlight the best candidate (lowest NLL, or lower OE upper if tied)
            if hasattr(row._region, "nll") and row._region.nll is not None:
                if best_candidate_idx is not None and row_idx == best_candidate_idx:
                    total_nll_str = f"{HIGHLIGHT}{row._region.nll:.4f} (BEST){RESET}"
                else:
                    total_nll_str = f"{BOLD}{UNDERLINE}{row._region.nll:.4f}{RESET}"
            else:
                total_nll_str = "NA"
            output_string += f"                {BOLD}{UNDERLINE}Total NLL: {total_nll_str}{RESET}\n\n"
            if hasattr(row._region, "region") and row._region.region:
                region_str = ", ".join([f"{r:7d}" for r in row._region.region])
                output_string += f"                Region residues: {region_str}\n"

                # Print per-residue observed and expected values
                if hasattr(row, "oe") and row.oe:
                    obs_list = []
                    exp_list = []
                    # Get excluded_residues if available
                    excluded_residues_dict = None
                    if (
                        hasattr(row, "excluded_residues")
                        and row.excluded_residues is not None
                    ):
                        excluded_residues_dict = row.excluded_residues
                    for res_idx in row._region.region:
                        # Check if this residue is excluded from stats
                        is_excluded = False
                        if excluded_residues_dict is not None:
                            # excluded_residues_dict is a dict-like structure
                            # Check if res_idx is in the dict and if its value is True
                            if hasattr(excluded_residues_dict, "get"):
                                exclude_flag = excluded_residues_dict.get(res_idx)
                                is_excluded = (
                                    exclude_flag is True
                                    if exclude_flag is not None
                                    else False
                                )
                            elif isinstance(excluded_residues_dict, dict):
                                is_excluded = (
                                    excluded_residues_dict.get(res_idx, False) is True
                                )

                        if is_excluded:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                        elif res_idx < len(row.oe):
                            oe_entry = row.oe[res_idx]
                            obs_val = (
                                oe_entry.obs
                                if hasattr(oe_entry, "obs") and oe_entry.obs is not None
                                else None
                            )
                            exp_val = (
                                oe_entry.exp
                                if hasattr(oe_entry, "exp") and oe_entry.exp is not None
                                else None
                            )

                            obs_list.append(
                                f"{obs_val:7d}" if obs_val is not None else "    NA"
                            )
                            exp_list.append(
                                f"{exp_val:7.2f}" if exp_val is not None else "    NA"
                            )
                        else:
                            obs_list.append("     NA")
                            exp_list.append("     NA")
                    output_string += f"                Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                    output_string += f"                Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"

            # Print updated null model after removing candidate region
            if hasattr(row, "_updated_null") and row._updated_null:
                output_string += f"\n                {BOLD}Updated Null Model:{RESET}\n"
                obs_str = (
                    f"{row._updated_null.obs:7d}"
                    if (
                        hasattr(row._updated_null, "obs")
                        and row._updated_null.obs is not None
                    )
                    else "     NA"
                )
                output_string += f"                    Observed: {obs_str}\n"
                exp_str = (
                    f"{row._updated_null.exp:7.2f}"
                    if (
                        hasattr(row._updated_null, "exp")
                        and row._updated_null.exp is not None
                    )
                    else "     NA"
                )
                output_string += f"                    Expected: {exp_str}\n"
                # Color aggregate O/E
                if (
                    hasattr(row._updated_null, "oe")
                    and row._updated_null.oe is not None
                ):
                    oe_str = f"{row._updated_null.oe:7.2f}"
                else:
                    oe_str = "     NA"
                output_string += f"                    O/E:      {oe_str}\n"
                # Color NLL values
                if (
                    hasattr(row._updated_null, "nll")
                    and row._updated_null.nll is not None
                ):
                    nll_str = f"{row._updated_null.nll:.4f}"
                else:
                    nll_str = "NA"
                output_string += f"                    NLL:      {nll_str}\n\n"
                if hasattr(row._updated_null, "region") and row._updated_null.region:
                    region_str = ", ".join(
                        [f"{r:7d}" for r in row._updated_null.region]
                    )
                    output_string += (
                        f"                    Region residues: {region_str}\n"
                    )

                    # Print per-residue observed and expected values for updated null
                    # model
                    if hasattr(row, "oe") and row.oe:
                        obs_list = []
                        exp_list = []
                        # Get excluded_residues if available
                        excluded_residues_dict = None
                        if (
                            hasattr(row, "excluded_residues")
                            and row.excluded_residues is not None
                        ):
                            excluded_residues_dict = row.excluded_residues
                        for res_idx in row._updated_null.region:
                            # Check if this residue is excluded from stats
                            is_excluded = False
                            if excluded_residues_dict is not None:
                                if hasattr(excluded_residues_dict, "get"):
                                    exclude_flag = excluded_residues_dict.get(res_idx)
                                    is_excluded = (
                                        exclude_flag is True
                                        if exclude_flag is not None
                                        else False
                                    )
                                elif isinstance(excluded_residues_dict, dict):
                                    is_excluded = (
                                        excluded_residues_dict.get(res_idx, False)
                                        is True
                                    )

                            if is_excluded:
                                obs_list.append("     NA")
                                exp_list.append("     NA")
                            elif res_idx < len(row.oe):
                                oe_entry = row.oe[res_idx]
                                obs_val = (
                                    oe_entry.obs
                                    if hasattr(oe_entry, "obs")
                                    and oe_entry.obs is not None
                                    else None
                                )
                                exp_val = (
                                    oe_entry.exp
                                    if hasattr(oe_entry, "exp")
                                    and oe_entry.exp is not None
                                    else None
                                )

                                obs_list.append(
                                    f"{obs_val:7d}" if obs_val is not None else "    NA"
                                )
                                exp_list.append(
                                    f"{exp_val:7.2f}"
                                    if exp_val is not None
                                    else "    NA"
                                )
                            else:
                                obs_list.append("     NA")
                                exp_list.append("     NA")
                        output_string += f"                    Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                        output_string += f"                    Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"

    # Print model comparison info if available (from first row)
    if hasattr(row, "aic_cand") and row.aic_cand is not None:
        output_string += (
            f"\n    {BOLD}Model Comparison ({model_comparison_method}):{RESET}\n"
        )
        output_string += f"      Candidate AIC: {row.aic_cand:.4f}\n"
        if hasattr(row, "p_lrt") and row.p_lrt is not None:
            output_string += f"      LRT p-value: {row.p_lrt:.6f}\n"
        if hasattr(row, "w_cand") and row.w_cand is not None:
            output_string += f"      AIC weight: {row.w_cand:.4f}\n"

    return output_string


def _debug_run_forward_round_start(
    ht,
    round_num,
    region_expr,
    model_comparison_method,
    debug_outputs,
    oe_upper_method: str = "gamma",
    min_exp_mis: int = 1,
):
    """Handle debug output for the start of a round."""
    debug_outputs.append(f"\n\n{BOLD}=== run_forward: Round {round_num} ==={RESET}\n")
    first_row = ht.head(1)
    if first_row.count() > 0:
        uniprot_id = first_row.uniprot_id.collect()[0]
        transcript_id = first_row.transcript_id.collect()[0]
        debug_outputs.append(
            f"\nUniProt ID: {uniprot_id}, Transcript ID: {transcript_id}\n"
        )
        debug_output = debug_print_forward_round(
            ht.annotate(_region=region_expr),
            uniprot_id,
            transcript_id,
            model_comparison_method,
            title=f"run_forward: Round {round_num} - After preparing candidate region",
            oe_upper_method=oe_upper_method,
            min_exp_mis=min_exp_mis,
        )
        if debug_output:
            debug_outputs.append(debug_output)


def _debug_run_forward_best_candidate(
    ht, ht2, round_num, debug_outputs, oe_upper_method: str = "gamma"
):
    """Handle debug output for the best candidate."""
    _ht_debug = ht.head(1)
    if _ht_debug.count() > 0:
        uniprot_id = _ht_debug.uniprot_id.collect()[0]
        transcript_id = _ht_debug.transcript_id.collect()[0]
        _ht_debug = ht.filter(
            (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
        )
        has_regions = _ht_debug.aggregate(hl.agg.any(hl.is_defined(_ht_debug.region)))
        if has_regions:
            best_region_debug = (
                ht2.filter(
                    (ht2.uniprot_id == uniprot_id)
                    & (ht2.transcript_id == transcript_id)
                )
                .select("best_region", "m_candidates", "min_idx")
                .collect()
            )
            _ht_debug_oe = _ht_debug.select("oe").collect()
            center_residue_idx = None
            if best_region_debug and best_region_debug[0].min_idx is not None:
                _ht_center = (
                    _ht_debug.filter(
                        (_ht_debug.uniprot_id == uniprot_id)
                        & (_ht_debug.transcript_id == transcript_id)
                        & (_ht_debug.idx == best_region_debug[0].min_idx)
                    )
                    .select("center_residue_index")
                    .collect()
                )
                if _ht_center and len(_ht_center) > 0:
                    center_residue_idx = _ht_center[0].center_residue_index
            if best_region_debug:
                br_row = best_region_debug[0]
                output_string = f"\n    {BOLD}=== run_forward: Round {round_num} - Best Candidate ==={RESET}\n\n"
                output_string += (
                    f"        Number of candidates evaluated: {br_row.m_candidates}\n"
                )
                if br_row.best_region:
                    center_res_str = (
                        f" (center residue: {center_residue_idx})"
                        if center_residue_idx is not None
                        else ""
                    )
                    output_string += f"\n        {BOLD}Best Candidate Region{center_res_str}:{RESET}\n"
                    output_string += (
                        f"            Observed: {br_row.best_region.obs:7d}\n"
                    )
                    output_string += (
                        f"            Expected: {br_row.best_region.exp:7.2f}\n"
                    )
                    if br_row.best_region.oe is not None:
                        oe_str = f"{br_row.best_region.oe:7.2f}"
                    else:
                        oe_str = "     NA"
                    output_string += f"            O/E:      {oe_str}\n"
                    # Calculate and display OE upper
                    if (
                        hasattr(br_row.best_region, "obs")
                        and br_row.best_region.obs is not None
                        and hasattr(br_row.best_region, "exp")
                        and br_row.best_region.exp is not None
                    ):
                        oe_upper_val = _calculate_oe_upper(
                            br_row.best_region.obs,
                            br_row.best_region.exp,
                            oe_upper_method,
                        )
                        oe_upper_str = (
                            f"{oe_upper_val:7.2f}"
                            if oe_upper_val is not None
                            else "     NA"
                        )
                    else:
                        oe_upper_str = "     NA"
                    output_string += f"            OE Upper: {oe_upper_str}\n"
                    if br_row.best_region.region_nll is not None:
                        region_nll_str = f"{br_row.best_region.region_nll:.4f}"
                    else:
                        region_nll_str = "NA"
                    output_string += f"            Region NLL: {region_nll_str}\n"
                    if br_row.best_region.nll is not None:
                        total_nll_str = f"{br_row.best_region.nll:.4f}"
                    else:
                        total_nll_str = "NA"
                    output_string += f"            {BOLD}{UNDERLINE}Total NLL: {total_nll_str}{RESET}\n\n"
                    if (
                        hasattr(br_row.best_region, "region")
                        and br_row.best_region.region
                    ):
                        region_str = ", ".join(
                            [f"{r:7d}" for r in br_row.best_region.region]
                        )
                        output_string += f"            Region residues: {region_str}\n"
                        if (
                            _ht_debug_oe
                            and len(_ht_debug_oe) > 0
                            and _ht_debug_oe[0].oe
                        ):
                            obs_list = []
                            exp_list = []
                            for res_idx in br_row.best_region.region:
                                if res_idx < len(_ht_debug_oe[0].oe):
                                    oe_entry = _ht_debug_oe[0].oe[res_idx]
                                    obs_val = (
                                        oe_entry.obs
                                        if hasattr(oe_entry, "obs")
                                        and oe_entry.obs is not None
                                        else None
                                    )
                                    exp_val = (
                                        oe_entry.exp
                                        if hasattr(oe_entry, "exp")
                                        and oe_entry.exp is not None
                                        else None
                                    )
                                    obs_list.append(
                                        f"{obs_val:7d}"
                                        if obs_val is not None
                                        else "    NA"
                                    )
                                    exp_list.append(
                                        f"{exp_val:7.2f}"
                                        if exp_val is not None
                                        else "    NA"
                                    )
                            output_string += f"            Observed ({len(obs_list):3d}):  {', '.join(obs_list)}\n"
                            output_string += f"            Expected ({len(exp_list):3d}):  {', '.join(exp_list)}\n"
                debug_outputs.append(output_string)


def _debug_run_forward_model_comparison(
    ht,
    ht2,
    round_num,
    model_comparison_method,
    lrt_alpha,
    lrt_df_added,
    bonferroni_per_round,
    aic_weight_thresh,
    debug_outputs,
):
    """Handle debug output for model comparison."""
    _ht_debug = ht.head(1)
    if _ht_debug.count() > 0:
        uniprot_id = _ht_debug.uniprot_id.collect()[0]
        transcript_id = _ht_debug.transcript_id.collect()[0]
        _ht_debug = ht.filter(
            (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
        )
        has_regions = _ht_debug.aggregate(hl.agg.any(hl.is_defined(_ht_debug.region)))
        if has_regions:
            _ht_debug2 = ht2.filter(
                (ht2.uniprot_id == uniprot_id) & (ht2.transcript_id == transcript_id)
            )
            _best_region_debug = _ht_debug2[
                (_ht_debug.uniprot_id, _ht_debug.transcript_id)
            ].best_region
            _m_candidates_debug = hl.or_else(
                _ht_debug2[
                    (_ht_debug.uniprot_id, _ht_debug.transcript_id)
                ].m_candidates,
                1,
            )
            _current_nll_debug = _ht_debug.selected_nll + _ht_debug.null_model.nll
            _candidate_nll_debug = hl.or_missing(
                hl.is_defined(_best_region_debug), _best_region_debug.nll
            )
            # Lazy import to avoid circular dependency
            from gnomad_constraint.experimental.proemis3d.utils import getAIC

            _aic_cand_debug = hl.or_missing(
                hl.is_defined(_best_region_debug),
                getAIC(_best_region_debug.null_model, 0)
                + getAIC(
                    _ht_debug.selected.append(_best_region_debug.drop("null_model")),
                    _best_region_debug.nll,
                ),
            )
            _found_best_debug = _ht_debug.found_best | hl.is_missing(_best_region_debug)
            if model_comparison_method == "lrt":
                _lrt_stat_debug = hl.or_missing(
                    hl.is_defined(_best_region_debug),
                    2 * (_current_nll_debug - _candidate_nll_debug),
                )
                _p_lrt_debug = hl.or_missing(
                    hl.is_defined(_best_region_debug),
                    hl.pchisqtail(_lrt_stat_debug, lrt_df_added),
                )
                _adj_alpha_debug = hl.if_else(
                    bonferroni_per_round,
                    lrt_alpha / hl.max(1, _m_candidates_debug),
                    lrt_alpha,
                )
                _found_best_debug |= _p_lrt_debug > _adj_alpha_debug
            elif model_comparison_method == "aic":
                _found_best_debug |= _aic_cand_debug >= _ht_debug.best_aic
            elif model_comparison_method == "aic_weight":
                _w_cand_debug = 1 / (
                    1 + hl.exp(0.5 * (_aic_cand_debug - _ht_debug.best_aic))
                )
                _found_best_debug |= _w_cand_debug < aic_weight_thresh

            debug_data = (
                _ht_debug.annotate(
                    current_nll=_current_nll_debug,
                    candidate_nll=_candidate_nll_debug,
                    aic_cand=_aic_cand_debug,
                    found_best=_found_best_debug,
                )
                .select(
                    "selected_nll",
                    "null_model",
                    "best_aic",
                    "current_nll",
                    "candidate_nll",
                    "aic_cand",
                    "found_best",
                )
                .collect()
            )
            debug_data2 = _ht_debug2.select(
                "best_region",
                "m_candidates",
            ).collect()
            if debug_data and debug_data2:
                row = debug_data[0]
                row2 = debug_data2[0]
                output_string = f"\n{BOLD}    === run_forward: Round {round_num} - Model Comparison ==={RESET}\n"
                output_string += f"        {BOLD}Current Model:{RESET}\n"
                if row.current_nll is not None:
                    current_nll_str = f"{BOLD}{UNDERLINE}{row.current_nll:.4f}{RESET}"
                else:
                    current_nll_str = "NA"
                output_string += f"            Current NLL: {current_nll_str}\n"
                if row.best_aic is not None:
                    best_aic_str = f"{row.best_aic:.4f}"
                else:
                    best_aic_str = "NA"
                output_string += f"            Best AIC: {best_aic_str}\n"
                if row2.best_region and row.candidate_nll is not None:
                    output_string += f"\n          {BOLD}Candidate Model:{RESET}\n"
                    if row.candidate_nll is not None:
                        candidate_nll_str = (
                            f"{BOLD}{UNDERLINE}{row.candidate_nll:.4f}{RESET}"
                        )
                    else:
                        candidate_nll_str = "NA"
                    output_string += f"            Candidate NLL: {candidate_nll_str}\n"
                    # For AIC method, Candidate AIC will be highlighted below
                    # For other methods, show it uncolored here
                    if row.aic_cand is not None and model_comparison_method != "aic":
                        if row.best_aic is not None and row.aic_cand < row.best_aic:
                            aic_cand_str = f"{GREEN}{row.aic_cand:.4f}{RESET}"
                        else:
                            aic_cand_str = f"{row.aic_cand:.4f}"
                        output_string += f"            Candidate AIC: {aic_cand_str}\n"
                    if model_comparison_method == "lrt":
                        lrt_stat = 2 * (row.current_nll - row.candidate_nll)
                        p_lrt_val = hl.eval(hl.pchisqtail(lrt_stat, lrt_df_added))
                        adj_alpha = (
                            lrt_alpha / max(1, row2.m_candidates)
                            if bonferroni_per_round
                            else lrt_alpha
                        )
                        output_string += f"            LRT statistic: {lrt_stat:.4f}\n"
                        # Highlight LRT p-value (the key comparison statistic)
                        accept = p_lrt_val <= adj_alpha
                        if accept:
                            p_lrt_str = (
                                f"{BOLD}{UNDERLINE}{GREEN}{p_lrt_val:.6f}{RESET}"
                            )
                        else:
                            p_lrt_str = f"{BOLD}{UNDERLINE}{RED}{p_lrt_val:.6f}{RESET}"
                        output_string += f"            LRT p-value: {p_lrt_str}\n"
                        output_string += (
                            f"            Adjusted alpha: {adj_alpha:.6f}\n"
                        )
                        if accept:
                            accept_str = f"{GREEN}{accept}{RESET}"
                        else:
                            accept_str = f"{RED}{accept}{RESET}"
                        output_string += f"            Accept candidate: {accept_str}\n"
                    elif model_comparison_method == "aic":
                        accept = (
                            row.aic_cand < row.best_aic
                            if row.aic_cand is not None
                            else False
                        )
                        # Highlight Candidate AIC (the key comparison statistic)
                        if accept:
                            aic_cand_str = (
                                f"{BOLD}{UNDERLINE}{GREEN}{row.aic_cand:.4f}{RESET}"
                            )
                        else:
                            aic_cand_str = (
                                f"{BOLD}{UNDERLINE}{RED}{row.aic_cand:.4f}{RESET}"
                            )
                        output_string += f"            Candidate AIC: {aic_cand_str}\n"
                        if accept:
                            accept_str = f"{GREEN}{accept}{RESET}"
                        else:
                            accept_str = f"{RED}{accept}{RESET}"
                        output_string += f"            Accept candidate: {accept_str}\n"
                    elif model_comparison_method == "aic_weight":
                        if row.aic_cand is not None:
                            w_cand_val = 1 / (
                                1 + np.exp(0.5 * (row.aic_cand - row.best_aic))
                            )
                            accept = w_cand_val >= aic_weight_thresh
                            # Highlight AIC weight (the key comparison statistic)
                            if accept:
                                w_cand_str = (
                                    f"{BOLD}{UNDERLINE}{GREEN}{w_cand_val:.4f}{RESET}"
                                )
                            else:
                                w_cand_str = (
                                    f"{BOLD}{UNDERLINE}{RED}{w_cand_val:.4f}{RESET}"
                                )
                            output_string += f"            AIC weight: {w_cand_str}\n"
                            if accept:
                                accept_str = f"{GREEN}{accept}{RESET}"
                            else:
                                accept_str = f"{RED}{accept}{RESET}"
                            output_string += (
                                f"            Accept candidate: {accept_str}\n"
                            )
                output_string += f"\n        Found best (stop): {row.found_best}\n"
                debug_outputs.append(output_string)


def _debug_collect_candidates_before_update(ht: hl.Table) -> list:
    """
    Collect candidate regions before they are updated/removed.

    This is used to show which candidates were removed when a region is selected.

    :param ht: Hail Table with candidate regions
    :return: List of candidate data (center_residue_index, region, oe) for the first uniprot/transcript
    """
    candidates_before_filter = []
    _ht_debug_before_update = ht.head(1)
    if _ht_debug_before_update.count() > 0:
        uniprot_id_before = _ht_debug_before_update.uniprot_id.collect()[0]
        transcript_id_before = _ht_debug_before_update.transcript_id.collect()[0]
        _ht_candidates_before = ht.filter(
            (ht.uniprot_id == uniprot_id_before)
            & (ht.transcript_id == transcript_id_before)
        )
        # Collect region (candidate regions) before they're removed, along with oe data
        candidates_before_filter = _ht_candidates_before.select(
            "center_residue_index", "region", "oe"
        ).collect()
    return candidates_before_filter


def _debug_run_forward_after_update(
    ht,
    round_num,
    candidates_before_filter,
    debug_outputs,
    oe_upper_method: str = "gamma",
):
    """Handle debug output after round update."""
    _ht_debug = ht.head(1)
    if _ht_debug.count() > 0:
        uniprot_id = _ht_debug.uniprot_id.collect()[0]
        transcript_id = _ht_debug.transcript_id.collect()[0]
        _ht_debug = ht.filter(
            (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
        )
        has_regions = _ht_debug.aggregate(hl.agg.any(hl.is_defined(_ht_debug.region)))
        if has_regions:
            debug_data = _ht_debug.select(
                "selected", "selected_nll", "best_aic", "found_best", "region"
            ).collect()
            if debug_data:
                row = debug_data[0]
                output_string = f"\n\n\n{BOLD}    === run_forward: Round {round_num} - After Update ==={RESET}\n"
                output_string += f"        Selected regions: {len(row.selected) if row.selected else 0}\n"
                if row.selected and len(row.selected) > 0:
                    latest_region = row.selected[-1]
                    output_string += f"\n        {BOLD}Newly Selected Region:{RESET}\n"
                    if hasattr(latest_region, "obs") and latest_region.obs is not None:
                        output_string += (
                            f"            Observed: {latest_region.obs:7d}\n"
                        )
                    if hasattr(latest_region, "exp") and latest_region.exp is not None:
                        output_string += (
                            f"            Expected: {latest_region.exp:7.2f}\n"
                        )
                    if hasattr(latest_region, "oe") and latest_region.oe is not None:
                        output_string += (
                            f"            O/E:      {latest_region.oe:7.2f}\n"
                        )
                    if hasattr(latest_region, "region") and latest_region.region:
                        region_str = ", ".join(
                            [f"{r:7d}" for r in latest_region.region]
                        )
                        output_string += f"            {BOLD}Residues ({len(latest_region.region)}): {region_str}{RESET}\n"
                if row.selected_nll is not None:
                    selected_nll_str = f"{BOLD}{UNDERLINE}{row.selected_nll:.4f}{RESET}"
                else:
                    selected_nll_str = "NA"
                output_string += f"\n            Selected NLL: {selected_nll_str}\n"
                if row.best_aic is not None:
                    best_aic_str = f"{row.best_aic:.4f}"
                else:
                    best_aic_str = "NA"
                output_string += f"            Best AIC: {best_aic_str}\n"
                output_string += f"            Found best: {row.found_best}\n"

                if row.selected and len(row.selected) > 0 and candidates_before_filter:
                    latest_region = row.selected[-1]
                    chosen_residues = (
                        set(latest_region.region)
                        if hasattr(latest_region, "region") and latest_region.region
                        else set()
                    )
                    output_string += f"\n        {BOLD}All Candidates from This Round (with chosen region residues colored red):{RESET}\n"
                    for cand_idx, cand_row in enumerate(candidates_before_filter):
                        if hasattr(cand_row, "region") and cand_row.region:
                            center_res_str = (
                                f" (center residue: {cand_row.center_residue_index})"
                                if hasattr(cand_row, "center_residue_index")
                                and cand_row.center_residue_index is not None
                                else ""
                            )
                            output_string += f"\n            {BOLD}Candidate {cand_idx + 1}{center_res_str}:{RESET}\n"
                            region_residues = cand_row.region
                            if hasattr(cand_row, "oe") and cand_row.oe:
                                total_obs = 0
                                total_exp = 0.0
                                for res_idx in region_residues:
                                    if res_idx < len(cand_row.oe):
                                        oe_entry = cand_row.oe[res_idx]
                                        if (
                                            hasattr(oe_entry, "obs")
                                            and oe_entry.obs is not None
                                        ):
                                            total_obs += oe_entry.obs
                                        if (
                                            hasattr(oe_entry, "exp")
                                            and oe_entry.exp is not None
                                        ):
                                            total_exp += oe_entry.exp
                                if total_exp > 0:
                                    total_oe = total_obs / total_exp
                                else:
                                    total_oe = None
                                output_string += (
                                    f"                Observed: {total_obs:7d}\n"
                                )
                                output_string += (
                                    f"                Expected: {total_exp:7.2f}\n"
                                )
                                if total_oe is not None:
                                    output_string += (
                                        f"                O/E:      {total_oe:7.2f}\n"
                                    )
                                # Calculate and display OE upper
                                if total_exp > 0:
                                    oe_upper_val = _calculate_oe_upper(
                                        total_obs, total_exp, oe_upper_method
                                    )
                                    oe_upper_str = (
                                        f"{oe_upper_val:7.2f}"
                                        if oe_upper_val is not None
                                        else "     NA"
                                    )
                                else:
                                    oe_upper_str = "     NA"
                                output_string += (
                                    f"                OE Upper: {oe_upper_str}\n"
                                )
                            residue_strs = []
                            for res in region_residues:
                                if res in chosen_residues:
                                    residue_strs.append(f"{RED}{res}{RESET}")
                                else:
                                    residue_strs.append(f"{res}")
                            region_str = ", ".join(residue_strs)
                            output_string += f"                {BOLD}Residues ({len(region_residues)}): {region_str}{RESET}\n"
                            remaining_residues = [
                                r for r in region_residues if r not in chosen_residues
                            ]
                            if remaining_residues:
                                remaining_str = ", ".join(
                                    [f"{r}" for r in remaining_residues]
                                )
                                output_string += f"                Remaining after removal: {remaining_str} ({len(remaining_residues)} residues)\n"
                            else:
                                output_string += f"                {RED}(All residues removed - candidate eliminated){RESET}\n"
                debug_outputs.append(output_string)


def _debug_run_forward_final(ht, debug_outputs):
    """Handle debug output for final state."""
    _ht_debug = ht.head(1)
    if _ht_debug.count() > 0:
        uniprot_id = _ht_debug.uniprot_id.collect()[0]
        transcript_id = _ht_debug.transcript_id.collect()[0]
        _ht_debug = ht.filter(
            (ht.uniprot_id == uniprot_id) & (ht.transcript_id == transcript_id)
        )
        debug_data = _ht_debug.select(
            "selected", "selected_nll", "best_aic", "null_model", "oe"
        ).collect()
        if debug_data:
            row = debug_data[0]
            output_string = f"\n{BOLD}=== run_forward: Final State ==={RESET}\n"
            output_string += f"    {BOLD}Total selected regions: {len(row.selected) if row.selected else 0}{RESET}\n"
            if row.selected_nll is not None:
                selected_nll_str = f"{BOLD}{UNDERLINE}{row.selected_nll:.4f}{RESET}"
            else:
                selected_nll_str = "NA"
            output_string += (
                f"    {BOLD}Final selected NLL: {selected_nll_str}{RESET}\n"
            )
            if row.best_aic is not None:
                best_aic_str = f"{row.best_aic:.4f}"
            else:
                best_aic_str = "NA"
            output_string += f"    {BOLD}Final best AIC: {best_aic_str}{RESET}\n"
            if row.selected:
                output_string += f"\n    {BOLD}Selected Regions:{RESET}\n"
                for idx, region in enumerate(row.selected):
                    output_string += f"        {BOLD}Region {idx}:{RESET}\n"
                    output_string += (
                        f"            {BOLD}Observed: {region.obs:7d}{RESET}\n"
                    )
                    output_string += (
                        f"            {BOLD}Expected: {region.exp:7.2f}{RESET}\n"
                    )
                    if region.oe is not None:
                        oe_str = f"{region.oe:7.2f}"
                    else:
                        oe_str = "     NA"
                    output_string += f"            {BOLD}O/E:      {oe_str}{RESET}\n"
                    if hasattr(region, "region") and region.region:
                        region_str = ", ".join([f"{r:7d}" for r in region.region])
                        output_string += f"            {BOLD}Residues ({len(region.region):3d}): {region_str}{RESET}\n"
                        if row.oe:
                            obs_list = []
                            exp_list = []
                            for res_idx in region.region:
                                if res_idx < len(row.oe):
                                    oe_entry = row.oe[res_idx]
                                    obs_val = (
                                        oe_entry.obs
                                        if hasattr(oe_entry, "obs")
                                        and oe_entry.obs is not None
                                        else None
                                    )
                                    exp_val = (
                                        oe_entry.exp
                                        if hasattr(oe_entry, "exp")
                                        and oe_entry.exp is not None
                                        else None
                                    )
                                    obs_list.append(
                                        f"{obs_val:7d}"
                                        if obs_val is not None
                                        else "    NA"
                                    )
                                    exp_list.append(
                                        f"{exp_val:7.2f}"
                                        if exp_val is not None
                                        else "    NA"
                                    )
                            output_string += f"            {BOLD}Observed ({len(obs_list):3d}): {', '.join(obs_list)}{RESET}\n"
                            output_string += f"            {BOLD}Expected ({len(exp_list):3d}): {', '.join(exp_list)}{RESET}\n"
            output_string += f"\n\n    {BOLD}Null Model (catch-all):{RESET}\n"
            output_string += f"        {BOLD}Observed: {row.null_model.obs:7d}{RESET}\n"
            output_string += (
                f"        {BOLD}Expected: {row.null_model.exp:7.2f}{RESET}\n"
            )
            if row.null_model.oe is not None:
                oe_str = f"{row.null_model.oe:7.2f}"
            else:
                oe_str = "     NA"
            output_string += f"        {BOLD}O/E:      {oe_str}{RESET}\n"
            # Add residues and per-residue observed/expected if available
            if hasattr(row.null_model, "region") and row.null_model.region and row.oe:
                region_str = ", ".join([f"{r:7d}" for r in row.null_model.region])
                output_string += f"        {BOLD}Residues ({len(row.null_model.region):3d}): {region_str}{RESET}\n"
                obs_list = []
                exp_list = []
                for res_idx in row.null_model.region:
                    if res_idx < len(row.oe):
                        oe_entry = row.oe[res_idx]
                        obs_val = (
                            oe_entry.obs
                            if hasattr(oe_entry, "obs") and oe_entry.obs is not None
                            else None
                        )
                        exp_val = (
                            oe_entry.exp
                            if hasattr(oe_entry, "exp") and oe_entry.exp is not None
                            else None
                        )
                        obs_list.append(
                            f"{obs_val:7d}" if obs_val is not None else "     NA"
                        )
                        exp_list.append(
                            f"{exp_val:7.2f}" if exp_val is not None else "     NA"
                        )
                    else:
                        obs_list.append("     NA")
                        exp_list.append("     NA")
                output_string += f"        {BOLD}Observed ({len(obs_list):3d}): {', '.join(obs_list)}{RESET}\n"
                output_string += f"        {BOLD}Expected ({len(exp_list):3d}): {', '.join(exp_list)}{RESET}\n"
            debug_outputs.append(output_string)


def debug_run_forward(stage: str, debug_outputs: list, **kwargs):
    """
    Unified debug function for run_forward algorithm.

    :param stage: Stage name - one of: "round_start", "best_candidate",
                  "model_comparison", "after_update", "final"
    :param debug_outputs: List to append debug output strings to
    :param kwargs: Stage-specific keyword arguments
    """
    if stage == "round_start":
        _debug_run_forward_round_start(
            kwargs["ht"],
            kwargs["round_num"],
            kwargs["region_expr"],
            kwargs["model_comparison_method"],
            debug_outputs,
            kwargs.get("oe_upper_method", "gamma"),
            kwargs.get("min_exp_mis", 1),
        )
    elif stage == "best_candidate":
        _debug_run_forward_best_candidate(
            kwargs["ht"],
            kwargs["ht2"],
            kwargs["round_num"],
            debug_outputs,
            kwargs.get("oe_upper_method", "gamma"),
        )
    elif stage == "model_comparison":
        _debug_run_forward_model_comparison(
            kwargs["ht"],
            kwargs["ht2"],
            kwargs["round_num"],
            kwargs["model_comparison_method"],
            kwargs["lrt_alpha"],
            kwargs["lrt_df_added"],
            kwargs["bonferroni_per_round"],
            kwargs["aic_weight_thresh"],
            debug_outputs,
        )
    elif stage == "after_update":
        _debug_run_forward_after_update(
            kwargs["ht"],
            kwargs["round_num"],
            kwargs["candidates_before_filter"],
            debug_outputs,
            kwargs.get("oe_upper_method", "gamma"),
        )
    elif stage == "final":
        _debug_run_forward_final(kwargs["ht"], debug_outputs)
    else:
        raise ValueError(f"Unknown debug stage: {stage}")


def debug_get_3d_residue(
    stage: str,
    dist_mat_expr: hl.expr.ArrayExpression,
    oe_expr: hl.expr.ArrayExpression,
    max_pae: Optional[float],
    min_plddt: Optional[float],
    pae_cutoff_method: str,
    plddt_cutoff_method: Optional[str],
    center_residue_index_expr: hl.expr.Int32Expression,
    debug_min_max: Optional[dict],
    debug_outputs_by_residue: dict,
    oe_expr_after_calc: Optional[hl.expr.ArrayExpression] = None,
):
    """
    Unified debug function for get_3d_residue.

    :param stage: Stage name - one of: "initial", "before_sorting", "after_sorting",
                  "pae_matrix_before_filter", "after_calculate_oe_upper", "final"
    :param dist_mat_expr: Distance matrix expression at current step
    :param oe_expr: Original OE expression
    :param max_pae: Maximum PAE cutoff
    :param min_plddt: Minimum pLDDT cutoff
    :param pae_cutoff_method: PAE cutoff method
    :param plddt_cutoff_method: pLDDT cutoff method
    :param center_residue_index_expr: Center residue index expression
    :param debug_min_max: Dict with min/max values for consistent coloring
    :param debug_outputs_by_residue: Dict to store debug outputs
    :param oe_expr_after_calc: OE expression after calculate_oe_upper (optional)
    """
    _debug_get_3d_residue(
        dist_mat_expr,
        oe_expr,
        max_pae,
        min_plddt,
        pae_cutoff_method,
        plddt_cutoff_method,
        center_residue_index_expr,
        debug_min_max,
        debug_outputs_by_residue,
        stage,
        oe_expr_after_calc,
    )
