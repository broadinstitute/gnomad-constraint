"""Hand-written layer of the constraint pipeline flowchart.

Everything here is prose the code cannot supply: what each stage is for, what a
step does to the data and what it costs, and the pitfalls worth knowing. The
structure (steps, flags, resources, call trees) is read from the code by
``generate_flowchart.py``, which also checks that every key below still names a
step, function, flag, or resource in the code.

Keys:

    - ``stages``: Bands of the map, in order. Every step is in exactly one.
    - ``tags``: Cost tags used by step ``ops``, with a label and one-line blurb.
    - ``key_lanes``: The key-shape ribbon. Each lane alternates ``key`` entries
      (a Table key and the Tables that carry it) with ``via`` entries (the
      operation that changes the key and its cost tag).
    - ``key_lanes_note``: Paragraph shown above the ribbon.
    - ``steps``: Keyed by step name (the key in ``get_constraint_resources``).
      ``summary`` is one line for the map; ``ops`` lists what the step does in
      order, each tagged with its cost; ``notes`` are pitfalls.
    - ``functions``: Keyed by function name, or dotted path when the name is
      not unique. ``notes`` are shown under the docstring.
    - ``resources``: Keyed by resource name. ``summary`` describes the Table and
      ``readme`` is its heading anchor in ``resources/README.md``.
    - ``flags``: Keyed by argparse destination. ``notes`` add to the help text.
    - ``downstream``: Consumers outside this repository.

Text may use ``code`` in double backticks.
"""

NOTES = {
    "stages": [
        {
            "id": "context",
            "name": "Context preparation",
            "note": (
                "Every possible SNV in the genome, annotated once with everything"
                " downstream needs: VEP, methylation, GERP, coverage, AN, gnomAD"
                " frequencies, and the regional ``adj_r`` corrections."
            ),
            "steps": [
                "prepare_context",
                "calculate_gerp_cutoffs",
                "preprocess_data",
            ],
        },
        {
            "id": "models",
            "name": "Mutation rate and models",
            "note": (
                "Fit the per-context mutation rate on rare intronic and intergenic genome"
                " variants, then calibrate it against synonymous exome variants with the"
                " plateau and coverage models."
            ),
            "steps": [
                "calculate_mutation_rate",
                "create_training_set",
                "build_models",
            ],
        },
        {
            "id": "apply",
            "name": "Expected counts",
            "note": (
                "Apply the models to every possible coding SNV and sum observed and"
                " expected counts per transcript and constraint group. Two routes reach"
                " the same constraint group Table: per variant (keeps the per-SNV"
                " expected counts) or aggregated (sums first, never writes per-variant"
                " rows)."
            ),
            "steps": [
                "apply_models_per_variant",
                "aggregate_per_variant_expected",
                "apply_models_aggregated",
                "aggregate_by_constraint_groups",
            ],
        },
        {
            "id": "metrics",
            "name": "Constraint metrics",
            "note": (
                "Per-transcript observed:expected ratios with confidence intervals,"
                " z-scores, pLI, constraint flags, gene quality metrics, and the LOEUF"
                " ranks and bins."
            ),
            "steps": [
                "compute_gene_quality_metrics",
                "compute_constraint_metrics",
            ],
        },
        {
            "id": "release",
            "name": "Release",
            "note": (
                "Reshape the metrics into the public release Table and TSVs, and release"
                " the mutation rate Table."
            ),
            "steps": [
                "prepare_release",
                "prepare_release_mutation_rate",
                "export_release_tsv",
            ],
        },
    ],
    "tags": {
        "ordered": {
            "label": "ordered join",
            "blurb": "Prefix key. Partition-local, no data movement.",
        },
        "shuffle": {
            "label": "shuffle",
            "blurb": "Re-key or non-prefix join. The whole left side moves.",
        },
        "agg": {
            "label": "aggregate",
            "blurb": "group_by. A shuffle, but onto far fewer rows.",
        },
        "explode": {
            "label": "explode",
            "blurb": "One row becomes several, e.g. one per transcript.",
        },
        "collect": {
            "label": "collect",
            "blurb": "Result comes back to the driver. Bounded by driver memory.",
        },
        "broadcast": {
            "label": "broadcast",
            "blurb": "Small Table or value sent to every worker as a literal.",
        },
        "filter": {
            "label": "filter",
            "blurb": "Row filter. Free, but changes what everything downstream sees.",
        },
        "rowexpr": {
            "label": "row expr",
            "blurb": "Pure annotation. No join, no action.",
        },
        "ckpt": {
            "label": "checkpoint",
            "blurb": "Written down, so the lineage above it stops replaying.",
        },
        "write": {
            "label": "write",
            "blurb": "A pipeline output is written.",
        },
    },
    "key_lanes": [
        {
            "tag": "Sites",
            "steps": [
                {
                    "key": "(locus, alleles)",
                    "note": "VEP context Table: every SNV",
                },
                {
                    "via": "split if needed, then keyed lookups on locus or the key",
                    "cost": "ordered",
                },
                {
                    "key": "(locus, alleles)",
                    "note": "annotated_context.ht",
                },
                {
                    "via": "row annotations only",
                    "cost": "rowexpr",
                },
                {
                    "key": "(locus, alleles)",
                    "note": "context.preprocessed.ht",
                },
            ],
        },
        {
            "tag": "Mutation rate",
            "steps": [
                {
                    "key": "(locus, alleles)",
                    "note": "context.preprocessed.ht",
                },
                {
                    "via": "group_by on the compute_mu counts",
                    "cost": "agg",
                },
                {
                    "key": "(context, ref, alt, methylation_level)",
                    "note": "mutation_rate.ht, about 100 rows",
                },
            ],
        },
        {
            "tag": "Models",
            "steps": [
                {
                    "key": "(locus, alleles)",
                    "note": "context.preprocessed.ht",
                },
                {
                    "via": "filter to training sites, group_by, mu broadcast",
                    "cost": "agg",
                },
                {
                    "key": (
                        "(context, ref, alt, methylation_level, cpg, transition, mutation_type,"
                        " mutation_type_model, genomic_region, build_model, exomes_coverage)"
                    ),
                    "note": "constraint_training.ht",
                },
                {
                    "via": "group_by, then aggregated to the driver",
                    "cost": "collect",
                },
                {
                    "key": "none",
                    "note": "plateau.he and coverage.he (expressions)",
                },
            ],
        },
        {
            "tag": "Per variant",
            "steps": [
                {
                    "key": "(locus, alleles)",
                    "note": "context.preprocessed.ht, one row per possible SNV",
                },
                {
                    "via": "explode transcript consequences",
                    "cost": "explode",
                },
                {
                    "key": "(locus, alleles), repeated",
                    "note": "per_variant_expected.ht, one row per consequence",
                },
                {
                    "via": "group-by",
                    "cost": "agg",
                },
                {
                    "key": (
                        "(annotation, modifier, gene, gene_id, transcript, canonical,"
                        " mane_select)"
                    ),
                    "note": "per_variant_expected.aggregated.ht",
                },
                {
                    "via": "group-by into constraint groups",
                    "cost": "agg",
                },
                {
                    "key": "(gene, gene_id, transcript, canonical, mane_select)",
                    "note": (
                        "constraint_group.ht, metrics and the release; the key never changes"
                        " again"
                    ),
                },
                {
                    "via": "flatten",
                    "cost": "rowexpr",
                },
                {
                    "key": "unkeyed",
                    "note": "constraint_metrics.tsv.bgz",
                },
            ],
        },
        {
            "tag": "Aggregated",
            "steps": [
                {
                    "key": "(locus, alleles)",
                    "note": "context.preprocessed.ht",
                },
                {
                    "via": "explode, then group-by",
                    "cost": "agg",
                },
                {
                    "key": (
                        "(context, ref, alt, methylation_level, mutation type, exomes_coverage,"
                        " apply_model, 7 VEP fields)"
                    ),
                    "note": "temporary; the models are applied here",
                },
                {
                    "via": "group-by",
                    "cost": "agg",
                },
                {
                    "key": (
                        "(annotation, modifier, gene, gene_id, transcript, canonical,"
                        " mane_select)"
                    ),
                    "note": (
                        "aggregated_expected.ht, which joins the per-variant lane at"
                        " constraint_group.ht"
                    ),
                },
            ],
        },
        {
            "tag": "Gene quality",
            "steps": [
                {
                    "key": "(locus, alleles)",
                    "note": "context.preprocessed.ht",
                },
                {
                    "via": "key_by locus, distinct",
                    "cost": "ordered",
                },
                {
                    "key": "(locus)",
                    "note": "joined onto GENCODE CDS bases",
                },
                {
                    "via": "explode transcript_id, group-by",
                    "cost": "agg",
                },
                {
                    "key": "(transcript)",
                    "note": (
                        "gene_quality_metrics.ht, looked up by transcript in"
                        " --compute-constraint-metrics"
                    ),
                },
            ],
        },
    ],
    "steps": {
        "prepare_context": {
            "summary": (
                "Annotate every possible SNV with VEP, methylation, GERP, coverage, AN,"
                " and frequencies"
            ),
            "ops": [
                {
                    "tag": "filter",
                    "text": (
                        "Coalesce the VEP context Table to 5,000 partitions (its partition count"
                        " broke Hail 0.2.133), then keep only chr20/X/Y or the test genes with"
                        " ``--test`` or ``--test-gene-list``."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Recompute the LOFTEE ``END_TRUNC`` filter on every transcript"
                        " consequence and reset ``lof`` to HC when no filter remains, LC"
                        " otherwise."
                    ),
                },
                {
                    "tag": "filter",
                    "text": (
                        "Split multiallelic rows if ``was_split`` is absent, keep chr1-22/X/Y,"
                        " trim the heptamer context to a trimer, and drop trimers that are not"
                        " pure ACGT."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Collapse strands (reverse-complement ``ref``, ``alt`` and ``context``"
                        " when ref is G or T), add ``cpg``, ``transition`` and ``mutation_type``,"
                        " and slim VEP to ``most_severe_consequence`` plus selected transcript"
                        " fields."
                    ),
                },
                {
                    "tag": "ordered",
                    "text": (
                        "Look up methylation (CpG only, 0-2 scale), GERP (missing set to 0),"
                        " exome and genome coverage, and AN by ``locus``. Look up exome and"
                        " genome frequencies and filters by ``(locus, alleles)``. All are"
                        " prefix-key lookups."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Add ``genomic_region`` and ``sfs_bin`` (from the exome adj AF). Keep"
                        " only mean and median coverage and the adj AN."
                    ),
                },
                {
                    "tag": "ordered",
                    "text": (
                        "Look up ``adj_r`` and ``syn_adj_r`` by locus in their interval-keyed"
                        " Tables, then index each dict by ``context``."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Write ``annotated_context.ht``: one row per possible SNV.",
                },
            ],
            "notes": [
                (
                    "``check_resource_existence`` does not cover three inputs. GERP is loaded"
                    " with ``hl.experimental.load_dataset``, and the ``adj_r`` and"
                    " ``syn_adj_r`` Tables are read from paths hard-coded in"
                    " ``resource_utils.py``."
                ),
                (
                    "Genome frequencies come from the v3.1.2 genomes sites Table even for"
                    " v4.x. It keeps the downsamplings that v4 dropped, and the mutation rate"
                    " needs them."
                ),
                (
                    "``adj_r``, ``syn_adj_r`` and ``sfs_bin`` are carried through to the"
                    " per-variant output, but no pipeline step applies them, so expected"
                    " counts are not adj_r-corrected. Both adj_r Tables are autosome builds"
                    " (``*_autosome.agg.ht``), so expect missing values on chrX and chrY."
                ),
            ],
        },
        "calculate_gerp_cutoffs": {
            "summary": (
                "Compute 5th/95th percentile GERP cutoffs for the mutation rate sites;"
                " writes nothing"
            ),
            "ops": [
                {
                    "tag": "filter",
                    "text": (
                        "Keep ``autosome_or_par`` rows of the annotated context Table."
                    ),
                },
                {
                    "tag": "collect",
                    "text": (
                        "``approx_quantiles`` of ``gerp`` at 0.05 and 0.95, returned to the"
                        " driver. Every possible SNV is a row, so each base counts three times,"
                        " and bases with no GERP score count as 0."
                    ),
                },
                {
                    "tag": "broadcast",
                    "text": (
                        "The two floats replace the ``--calculate-mutation-rate-gerp-*-cutoff``"
                        " values for the rest of this command only."
                    ),
                },
            ],
            "notes": [
                (
                    "The cutoffs take effect only when ``--preprocess-data`` runs in the same"
                    " command. They are logged and recorded in ``calculate_mu_globals`` of"
                    " the preprocessed Table, and nowhere else."
                ),
                (
                    "The percentiles are taken over autosome and PAR sites, but the mutation"
                    " rate filter keeps autosomes only (``in_autosome``, which excludes PAR)."
                ),
            ],
        },
        "preprocess_data": {
            "summary": (
                "Add mutation rate, observed/possible, and model-group annotations to"
                " every site"
            ),
            "ops": [
                {
                    "tag": "filter",
                    "text": "With ``--test``, keep only chr20/X/Y or the test genes.",
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "``compute_mu``: trim the genome frequency array to the full adj entry"
                        " plus the downsampling-1000 entries, then set observed/possible only at"
                        " mutation rate sites. Those are autosomal, have genome mean coverage"
                        " 15-60, GERP strictly between the cutoffs, and an intron or intergenic"
                        " most severe consequence, and are either absent from the genomes or PASS"
                        " with AC <= 5 in the global downsampling of 1,000."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "``exomes_coverage`` from ``--exome-coverage-metric``. The default,"
                        " ``AN_percent``, is AN as a percentage of all possible alleles, with"
                        " sex-aware counts on chrX and chrY non-PAR."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "``calibrate_mu``: exome observed/possible per frequency entry. Filtered"
                        " variants count as possible but not observed; PASS variants above"
                        " ``--max-af`` count as neither."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "``build_model`` (synonymous in an Ensembl MANE Select transcript only)"
                        " and ``apply_model`` (every site): a high or low coverage bucket plus a"
                        " ``{cpg, genomic_region}`` model group."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Drop ``freq`` and replace all globals with ``calculate_mu_globals``,"
                        " ``build_models_globals``, ``apply_models_globals`` and the exome"
                        " frequency metadata. The VEP, AN and frequency globals are dropped."
                    ),
                },
                {
                    "tag": "write",
                    "text": (
                        "Write ``context.preprocessed.ht``. A row-wise pass over every possible"
                        " SNV with no joins or shuffles."
                    ),
                },
            ],
            "notes": [
                (
                    "Most model settings are fixed in this step, not the steps that use them,"
                    " and are stored in the Table: the coverage cutoffs, the low coverage"
                    " filter, ``--max-af``, ``--exome-coverage-metric``,"
                    " ``--genetic-ancestry-groups``, ``--include-downsamplings``, the"
                    " mutation rate site filters, and ``--skip-coverage-model``. Changing any"
                    " of them means rerunning from here."
                ),
                (
                    "Unobserved exome variants get ``observed_variants`` of 0, not missing,"
                    " so a group with nothing observed sums to 0. Sites with undefined exome"
                    " coverage get missing counts and drop out of every sum."
                ),
                (
                    "``build_model`` falls back to a low cutoff of 0 (strict ``>``), so"
                    " zero-coverage sites never train. ``apply_model`` has no lower bound"
                    " unless ``--pipeline-low-coverage-filter`` is set."
                ),
                (
                    "``build_model`` keeps ``cpg`` only in the high-coverage bucket, so the"
                    " low-coverage training bucket is split by region alone. ``apply_model``"
                    " keeps ``cpg`` in both buckets."
                ),
            ],
        },
        "calculate_mutation_rate": {
            "summary": (
                "Fit the per-context mutation rate on rare intronic and intergenic genome"
                " variants"
            ),
            "ops": [
                {
                    "tag": "agg",
                    "text": (
                        "group_by ``(context, ref, alt, methylation_level)`` over the whole"
                        " preprocessed Table, summing the ``compute_mu`` observed arrays and"
                        " possible counts. Rows outside the mutation rate sites have missing"
                        " counts and add nothing."
                    ),
                },
                {
                    "tag": "ckpt",
                    "text": "Checkpoint the small grouped Table to a temp path.",
                },
                {
                    "tag": "collect",
                    "text": (
                        "``total_bases = sum(possible) // 3`` is returned to the driver."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Per frequency entry i, ``mu[i] = 1.2e-8 * total_bases * (observed[i] /"
                        " possible) / sum(observed[i])``. ``mu_snp`` is the entry for the global"
                        " downsampling of 1,000. Mutation type annotations are added."
                    ),
                },
                {
                    "tag": "write",
                    "text": (
                        "Repartition to ``--mutation-rate-partitions`` (default 1) and write"
                        " ``mutation_rate.ht``."
                    ),
                },
            ],
            "notes": [
                (
                    "The rate is fit on autosomes only (PAR excluded) but applied everywhere,"
                    " including chrX and chrY."
                ),
                (
                    "Only ``mu_snp`` is used downstream. The ``mu`` array also holds the"
                    " full-genome adj entry and any per-ancestry downsampling-1000 entries."
                ),
                (
                    "``total_bases`` assumes all three alleles at a qualifying base pass the"
                    " filters. Alleles dropped one at a time (AC > 5, failing filters) make"
                    " it a slight undercount of bases."
                ),
            ],
        },
        "create_training_set": {
            "summary": (
                "Count observed and possible synonymous exome variants per context and"
                " coverage"
            ),
            "ops": [
                {
                    "tag": "filter",
                    "text": (
                        "Flatten ``calibrate_mu`` and keep rows with a ``build_model`` bucket and"
                        " ``possible_variants > 0``: synonymous in an Ensembl MANE Select"
                        " transcript, in either coverage bucket."
                    ),
                },
                {
                    "tag": "ckpt",
                    "text": (
                        "Select only the grouping and count fields, then checkpoint to a temp"
                        " path to avoid memory and shuffle problems."
                    ),
                },
                {
                    "tag": "agg",
                    "text": (
                        "group_by context, ref, alt, methylation level, mutation type,"
                        " ``genomic_region``, ``build_model`` and ``exomes_coverage``, summing"
                        " observed arrays and possible counts."
                    ),
                },
                {
                    "tag": "broadcast",
                    "text": (
                        "Attach ``mu_snp``. The mutation rate Table is collected into a dict and"
                        " sent out as a literal instead of being joined."
                    ),
                },
                {
                    "tag": "write",
                    "text": (
                        "Repartition to ``--training-set-partition-hint``, checkpoint to"
                        " ``constraint_training.ht``, and export the TSV."
                    ),
                },
            ],
            "notes": [
                (
                    "chrX and chrY non-PAR sites are included. They get their own plateau"
                    " models through ``genomic_region`` in the model group."
                ),
                (
                    "``annotate_with_mu`` raises 'Missing mu' if any training row's"
                    " ``(context, ref, alt, methylation_level)`` is not in the mutation rate"
                    " Table, for example one built from a different subset."
                ),
            ],
        },
        "build_models": {
            "summary": (
                "Fit plateau models per CpG and region, and the low-coverage correction"
                " model"
            ),
            "ops": [
                {
                    "tag": "rowexpr",
                    "text": (
                        "Print ``build_models_globals`` and raise if ``--skip-coverage-model``"
                        " differs from the setting used in preprocessing."
                    ),
                },
                {
                    "tag": "agg",
                    "text": (
                        "group_by context, ref, alt, methylation level, mutation type and"
                        " ``build_model``. Exome coverage is kept as a grouping only for"
                        " low-coverage rows, so high-coverage rows collapse across coverage."
                    ),
                },
                {
                    "tag": "collect",
                    "text": (
                        "One aggregation returns everything to the driver. It fits a plateau"
                        " regression of observed/possible on ``[1, mu_snp]`` per ``{cpg,"
                        " genomic_region}`` and per exome frequency entry (weighted by possible"
                        " with ``--use-weights``). It also computes the high-coverage scale"
                        " factor and the per-coverage sums for low-coverage autosome and PAR"
                        " rows."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Coverage model, computed from the collected values with no extra pass."
                        " ``oe = observed / (scale * sum(possible * mu_snp))`` per coverage"
                        " level, then an unweighted regression of ``oe`` on ``[1, coverage]``"
                        " (log10 with ``--use-logarithmic-coverage-model``)."
                    ),
                },
                {
                    "tag": "write",
                    "text": (
                        "Write ``plateau.he``, and ``coverage.he`` unless the coverage model is"
                        " skipped."
                    ),
                },
            ],
            "notes": [
                (
                    "Coverage buckets come from the ``build_model`` annotation made in"
                    " ``--preprocess-data``. ``build_models``' own ``high_cov_definition``"
                    " and ``upper_cov_cutoff`` arguments are ignored because the pipeline"
                    " passes ``model_group_expr``."
                ),
                (
                    "chrX and chrY non-PAR get their own plateau models. The coverage model"
                    " and its scale factor use autosome and PAR sites only, and only the"
                    " full-dataset observed counts."
                ),
                (
                    "Only ``--skip-coverage-model`` is checked between steps. ``coverage.he``"
                    " does not record the log10 choice, so pass the same"
                    " ``--use-logarithmic-coverage-model`` to the apply steps."
                ),
            ],
        },
        "apply_models_per_variant": {
            "summary": "Expected count for every possible SNV x transcript consequence",
            "ops": [
                {
                    "tag": "rowexpr",
                    "text": (
                        "Promote the ``calibrate_mu`` fields (observed and possible counts,"
                        " ``build_model``, ``apply_model``) to the top level."
                    ),
                },
                {
                    "tag": "filter",
                    "text": (
                        "Keep rows with a defined ``apply_model`` and ``possible_variants > 0``."
                        " This drops sites in neither coverage model and variants with exome AF"
                        " above ``--max-af``."
                    ),
                },
                {
                    "tag": "explode",
                    "text": (
                        "Explode ``vep.transcript_consequences``: one row per (variant,"
                        " transcript consequence), keeping the (locus, alleles) key. Adds"
                        " ``annotation``, ``modifier``, ``gene``, ``gene_id``, ``transcript``,"
                        " ``canonical``, ``mane_select`` and drops the ``vep`` struct."
                    ),
                },
                {
                    "tag": "broadcast",
                    "text": (
                        "Look up ``mu_snp`` in a dict built from the mutation rate Table, keyed"
                        " by (context, ref, alt, methylation_level), so there is no join. A"
                        " context missing from that Table fails with ``Missing mu``."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Apply the models: ``ppo = mu_snp * slope + intercept`` per freq entry,"
                        " from the plateau model of the row's ``apply_model.model_group``."
                        " ``coverage_correction`` is 0 at zero coverage, the coverage model for"
                        " low-coverage rows, and 1 otherwise. ``expected_variants = ppo *"
                        " possible * coverage_correction``."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Store ``plateau_models``, ``coverage_model``, ``log10_coverage`` and"
                        " ``groupings`` (a tuple) in ``apply_models_globals``. Drop the promoted"
                        " copies, so observed and possible counts are only in ``calibrate_mu``."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Write ``per_variant_expected.ht``.",
                },
            ],
            "notes": [
                (
                    "Several rows share each (locus, alleles) key, one per transcript"
                    " consequence, including non-coding ones. The coding filter is applied"
                    " only in ``--aggregate-per-variant-expected``."
                ),
                (
                    "Observed and possible counts are only inside ``calibrate_mu``. Flatten"
                    " with ``ht.annotate(**ht.calibrate_mu)`` before summing."
                ),
                (
                    "``adj_r`` and ``syn_adj_r`` are carried on each row, but no step applies"
                    " or sums them. Expected counts have no regional correction."
                ),
                (
                    "``use_mane_select=True`` is hard-coded, so ``--custom-vep-annotation"
                    " worst_csq_by_gene`` raises a ValueError here."
                ),
                (
                    "With ``--skip-coverage-model``, sites below"
                    " ``--apply-model-high-cov-definition`` (default 90) have no apply model."
                    " They are dropped from both observed and expected counts."
                ),
            ],
        },
        "aggregate_per_variant_expected": {
            "summary": "Sum per-variant counts per transcript, consequence and modifier",
            "ops": [
                {
                    "tag": "collect",
                    "text": (
                        "Evaluate ``apply_models_globals.groupings`` on the driver to get the"
                        " seven group-by fields."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Promote ``calibrate_mu`` to recover the per-variant"
                        " ``observed_variants`` and ``possible_variants``."
                    ),
                },
                {
                    "tag": "filter",
                    "text": "Keep coding consequences only (``CSQ_CODING``).",
                },
                {
                    "tag": "ckpt",
                    "text": (
                        "Unkey, keep the groupings and ``AGGREGATE_SUM_FIELDS``, and checkpoint"
                        " to a temp file."
                    ),
                },
                {
                    "tag": "agg",
                    "text": (
                        "Group by the groupings and sum ``mu_snp``, ``mu``,"
                        " ``observed_variants``, ``possible_variants``,"
                        " ``predicted_proportion_observed``, ``coverage_correction`` and"
                        " ``expected_variants``. Arrays are summed per freq entry."
                    ),
                },
                {
                    "tag": "ckpt",
                    "text": "Checkpoint again, then ``naive_coalesce(1000)``.",
                },
                {
                    "tag": "write",
                    "text": "Write ``per_variant_expected.aggregated.ht``.",
                },
            ],
            "notes": [
                (
                    "The step declares the mutation rate Table and the models as inputs but"
                    " never reads them. ``check_resource_existence`` still requires them to"
                    " exist."
                ),
                (
                    "Summed ``mu_snp``, ``predicted_proportion_observed`` and"
                    " ``coverage_correction`` are totals over variants. Divide by"
                    " ``possible_variants`` for a group mean."
                ),
            ],
        },
        "apply_models_aggregated": {
            "summary": "Sum counts first, then apply the models; no per-variant rows",
            "ops": [
                {
                    "tag": "filter",
                    "text": (
                        "Same preparation as the per-variant route: promote ``calibrate_mu``,"
                        " then keep rows with a defined ``apply_model`` and ``possible_variants >"
                        " 0``."
                    ),
                },
                {
                    "tag": "explode",
                    "text": (
                        "Explode ``vep.transcript_consequences`` into the seven VEP grouping"
                        " fields."
                    ),
                },
                {
                    "tag": "filter",
                    "text": "Keep coding consequences only (``CSQ_CODING``).",
                },
                {
                    "tag": "agg",
                    "text": (
                        "Group by (context, ref, alt, methylation_level, cpg, transition,"
                        " mutation_type, mutation_type_model, exomes_coverage, apply_model) plus"
                        " the VEP groupings. Sum observed and possible counts, with partition"
                        " hint 100."
                    ),
                },
                {
                    "tag": "broadcast",
                    "text": "Look up ``mu_snp`` per group from the mutation rate dict.",
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Apply the plateau and coverage models to the group counts, then"
                        " checkpoint to a temp file."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Multiply ``mu_snp``, ``predicted_proportion_observed`` and"
                        " ``coverage_correction`` by ``possible_variants``, so the sums match the"
                        " per-variant route."
                    ),
                },
                {
                    "tag": "agg",
                    "text": (
                        "Group by the seven VEP groupings and sum ``AGGREGATE_SUM_FIELDS``."
                        " ``groupings`` is stored as an array in ``apply_models_globals``."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Write ``aggregated_expected.ht``.",
                },
            ],
            "notes": [
                (
                    "Expected counts equal the per-variant route's. Every variant in a"
                    " first-pass group shares ``mu_snp``, coverage and model bucket, and the"
                    " models are linear in those."
                ),
                (
                    "Same schema as ``per_variant_expected.aggregated.ht``, except"
                    " ``apply_models_globals.groupings`` is an array, not a tuple."
                ),
                (
                    "Writes no per-variant rows. Anything that needs per-SNV expected counts"
                    " needs the per-variant route."
                ),
                (
                    "Only ``--aggregate-by-constraint-groups --use-aggregated-expected``"
                    " reads this output."
                ),
            ],
        },
        "aggregate_by_constraint_groups": {
            "summary": "Per-transcript observed and expected sums for syn, mis and LoF",
            "ops": [
                {
                    "tag": "rowexpr",
                    "text": (
                        "Build one boolean filter per constraint group from ``annotation`` and"
                        " ``modifier``: ``syn``, ``mis``, ``lof_classic`` (stop_gained, splice"
                        " donor/acceptor), ``lof_hc_lc`` (LOFTEE HC or LC) and ``lof_hc`` (LOFTEE"
                        " HC). Store the definitions in ``constraint_group_meta``."
                    ),
                },
                {
                    "tag": "ckpt",
                    "text": "Checkpoint the filters to a temp file.",
                },
                {
                    "tag": "agg",
                    "text": (
                        "Group by (gene, gene_id, transcript, canonical, mane_select) and sum the"
                        " fields separately for each group, using ``hl.agg.array_agg`` with"
                        " ``hl.agg.filter``."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "``no_variants``: the adj observed count, summed over all groups, is 0."
                    ),
                },
                {
                    "tag": "filter",
                    "text": (
                        "Drop transcripts with no observed variants and no expected variants in"
                        " any group."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Reshape each group's per-freq arrays into an ``oe_info`` array of"
                        " structs."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Write ``constraint_group.ht``.",
                },
            ],
            "notes": [
                (
                    "The groups overlap: ``lof_hc`` is a subset of ``lof_hc_lc``, and"
                    " ``lof_classic`` goes by consequence and ignores LOFTEE. A variant"
                    " counts in every group it matches."
                ),
                (
                    "The output key is the input key filtered to ``RELEASE_KEY_ORDER``, so"
                    " both routes give (gene, gene_id, transcript, canonical, mane_select)."
                ),
                (
                    "Both routes write the same ``constraint_group.ht`` path. Running one"
                    " route overwrites the other's result."
                ),
            ],
        },
        "compute_gene_quality_metrics": {
            "summary": "Per-transcript CDS coverage, mapping quality, segdup and LCR",
            "ops": [
                {
                    "tag": "agg",
                    "text": (
                        "Build or read the per-locus GENCODE CDS Table (``get_gencode_cds_ht``):"
                        " CDS features without chrY PAR copies, one row per base, grouped by"
                        " locus. It is checkpointed under ``gs://gnomad-tmp`` and reused if it"
                        " already exists."
                    ),
                },
                {
                    "tag": "ordered",
                    "text": (
                        "Coverage: reduce the preprocessed context to one row per locus, join"
                        " ``exomes_coverage`` onto CDS bases by locus, and keep bases with"
                        " defined coverage."
                    ),
                },
                {
                    "tag": "agg",
                    "text": (
                        "Explode ``transcript_id`` and group by transcript: ``prop_bp_AN90`` is"
                        " the fraction of CDS bases with ``exomes_coverage >= 90``."
                    ),
                },
                {
                    "tag": "agg",
                    "text": (
                        "Mapping quality: take ``info.AS_MQ`` from exome-site SNVs, look up the"
                        " overlapping CDS transcripts by locus, explode, and average per"
                        " transcript (``mean_AS_MQ``)."
                    ),
                },
                {
                    "tag": "agg",
                    "text": (
                        "Region overlap: look up each CDS base in the segdup and LCR interval"
                        " Tables, explode, and take per-transcript fractions (``prop_segdup``,"
                        " ``prop_LCR``)."
                    ),
                },
                {
                    "tag": "ckpt",
                    "text": (
                        "Each of the three Tables is ``.cache()``d (written to the Hail tmp dir)"
                        " and joined on transcript."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Pack the metrics into ``gene_quality_metrics`` (prefixed ``exome_``) and"
                        " set ``gene_flags``: ``low_exome_mapping_quality`` (mean AS_MQ < 50) and"
                        " ``low_exome_coverage`` (``prop_bp_AN90`` < 0.1)."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Write ``gene_quality_metrics.ht``, keyed by transcript.",
                },
            ],
            "notes": [
                (
                    "``prop_bp_AN90`` compares ``exomes_coverage`` with 90. That is an AN"
                    " percent only with the default ``--exome-coverage-metric AN_percent``."
                ),
                (
                    "``gene_flags`` is missing, not empty, when either input metric is"
                    " missing. For example, a transcript with no exome SNV sites has no"
                    " ``mean_AS_MQ``."
                ),
                (
                    "Transcripts with no CDS base of defined coverage are absent: the"
                    " coverage Table is the left side of the join."
                ),
                (
                    "The GENCODE CDS, segdup and LCR Tables are read directly, not declared,"
                    " so ``check_resource_existence`` does not check them."
                ),
            ],
        },
        "compute_constraint_metrics": {
            "summary": "oe, CIs, z-scores, flags, pLI, then LOEUF ranks and bins",
            "ops": [
                {
                    "tag": "rowexpr",
                    "text": (
                        "Re-read ``constraint_group.ht`` into 10,000 partitions. For every group"
                        " and freq entry, compute ``oe = obs / exp``, gamma and"
                        " discretized-Poisson 90% CIs, and ``z_raw``."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Set per-group ``flags`` from the adj entry: ``no_exp_*`` (expected"
                        " missing or 0) and ``outlier_*`` (raw z beyond the thresholds). Then a"
                        " temp checkpoint."
                    ),
                },
                {
                    "tag": "collect",
                    "text": (
                        "Compute ``sd_raw_z`` per group in one aggregation, over transcripts"
                        " without ``no_variants`` and without that group's flags. Then ``z_score"
                        " = z_raw / sd``, and ``constraint_flags`` is the union of the syn, mis"
                        " and lof_hc flags. Temp checkpoint."
                    ),
                },
                {
                    "tag": "collect",
                    "text": (
                        "pLI, pNull and pRec by EM for ``lof_hc`` and ``lof_hc_lc``: one"
                        " aggregation per iteration, until the LI weight changes by less than"
                        " ``--min-diff-convergence``. Temp checkpoint."
                    ),
                },
                {
                    "tag": "shuffle",
                    "text": (
                        "Look up gene quality metrics, then GENCODE transcript and CDS"
                        " annotations, by transcript. These are non-prefix lookups on a"
                        " per-transcript Table, so the shuffles are small."
                    ),
                },
                {
                    "tag": "write",
                    "text": (
                        "Write ``constraint_metrics_pre_rank.ht``. Skipped with"
                        " ``--skip-pre-rank-metrics``."
                    ),
                },
                {
                    "tag": "shuffle",
                    "text": (
                        "Rank the upper bound of both CIs for every group (10 metrics, one global"
                        " sort each), over ENST MANE Select per gene, or canonical where a gene"
                        " has none. Bins are ``int(rank * k / n)``, where n counts the ranked"
                        " transcripts with a defined value. Each step is ``.cache()``d to the tmp"
                        " dir."
                    ),
                },
                {
                    "tag": "collect",
                    "text": (
                        "Compute percentile, decile and sextile thresholds of the gamma upper CI"
                        " for syn, mis and lof_hc with ``approx_quantiles`` (k=1000), on ENST"
                        " MANE Select transcripts without constraint flags. Stored in the"
                        " ``percentile_thresholds`` global. Every transcript gets threshold-based"
                        " ``constraint_bins``."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Write ``constraint_metrics.ht``.",
                },
            ],
            "notes": [
                (
                    "Two binnings exist side by side. Rank-based ``oe_ci_*_rank.upper_bin_*``"
                    " covers the ranked set; threshold-based ``constraint_bins`` uses"
                    " thresholds from flag-free MANE Select transcripts. The release keeps"
                    " only the rank-based LOEUF bins plus the LoF thresholds."
                ),
                (
                    "Raw z is positive when there are fewer variants than expected. The lower"
                    " thresholds flag excess variation; the syn upper threshold flags strong"
                    " depletion."
                ),
                (
                    "The pLI EM fits its mixture weights over every transcript row with"
                    " expected > 0, not only MANE Select or canonical."
                ),
                (
                    "Ranks are 0-based and bins run from 0 to k-1."
                    " ``--use-mane-select-over-canonical`` changes only which transcripts are"
                    " ranked, not the thresholds."
                ),
                (
                    "``get_gencode_ht`` re-imports the GTF and checkpoints it to a fixed"
                    " ``gs://gnomad-tmp/gencode_v39_with_versions.ht``, which is read if"
                    " present. It is not a declared input."
                ),
            ],
        },
        "prepare_release": {
            "summary": "Reshape the metrics into the public release Table",
            "ops": [
                {
                    "tag": "rowexpr",
                    "text": (
                        "Rename the GENCODE fields: ``transcript_id_version`` becomes"
                        " ``transcript_version``, and ``level`` becomes ``transcript_level``."
                    ),
                },
                {
                    "tag": "collect",
                    "text": (
                        "Evaluate ``constraint_group_meta``, ``exomes_freq_meta`` and"
                        " ``sd_raw_z`` on the driver to drive the reshaping."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Make each group a top-level struct at the adj freq entry, with ``mu``"
                        " (summed ``mu_snp``), ``possible``, ``obs``, ``exp``, ``oe``, ``z_raw``,"
                        " ``z_score`` and ``oe_ci`` (gamma lower and upper). ``lof_hc`` becomes"
                        " ``lof`` and also gets ``upper_rank`` and ``upper_bin_*``. ``lof`` and"
                        " ``lof_hc_lc`` get pLI, pNull and pRec."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "If downsamplings were computed, add ``gen_anc_obs`` and ``gen_anc_exp``:"
                        " per genetic ancestry, arrays in downsampling order."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Select the release fields: GENCODE scalars, ``gene_quality_metrics``,"
                        " ``gene_flags``, ``constraint_flags``, ``syn``, ``mis``, ``lof_hc_lc``"
                        " and ``lof``."
                    ),
                },
                {
                    "tag": "filter",
                    "text": (
                        "Drop transcripts with ``possible == 0`` in all four groups."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Replace the globals with ``*_params`` (pipeline globals minus internal"
                        " fields), a named ``sd_raw_z`` struct, ``loeuf_percentile_thresholds``,"
                        " ``downsamplings``, ``max_af`` and ``version``."
                    ),
                },
                {
                    "tag": "write",
                    "text": (
                        "Apply ``naive_coalesce(1000)`` and write the release Table."
                    ),
                },
            ],
            "notes": [
                (
                    "Release ``mu`` is the summed ``mu_snp``, without the coverage"
                    " correction. The internal coverage-corrected ``mu`` is dropped."
                ),
                (
                    "Not released: ``lof_classic``, per-group ``flags``, the"
                    " discretized-Poisson CI, ``constraint_bins``, mis and syn ranks and"
                    " thresholds, and genetic ancestry entries without downsamplings."
                ),
                (
                    "No step sets a ``version`` global, so the release has one only if"
                    " ``--release-version`` is passed."
                ),
            ],
        },
        "prepare_release_mutation_rate": {
            "summary": "Public mutation rate Table and TSV",
            "ops": [
                {
                    "tag": "rowexpr",
                    "text": (
                        "Keep ``mu`` (``mu_snp``, the scalar rate at the canonical downsampling),"
                        " ``cpg``, ``transition`` and ``mutation_type``. Drop the counts, the"
                        " per-freq ``mu`` array and ``mutation_type_model``."
                    ),
                },
                {
                    "tag": "rowexpr",
                    "text": (
                        "Replace the globals with ``calculate_mu_params`` (minus ``freq_meta``,"
                        " ``genetic_ancestry_groups`` and ``downsampling_idx``) and ``version``."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Write ``release/mutation_rate.ht``.",
                },
                {
                    "tag": "write",
                    "text": "Export the same rows, key fields included, as a TSV.",
                },
            ],
            "notes": [
                ("Without ``--release-version`` the Table has no ``version`` global."),
            ],
        },
        "export_release_tsv": {
            "summary": "Flat release TSV, LOEUF thresholds TSV, downsampling TSV",
            "ops": [
                {
                    "tag": "rowexpr",
                    "text": (
                        "``--export-release-tsv``: unkey, put the key fields first, and"
                        " ``flatten`` nested structs into dotted columns (for example"
                        " ``lof.oe_ci.upper``). Downsampling arrays stay as array columns."
                    ),
                },
                {
                    "tag": "write",
                    "text": "Export ``constraint_metrics.tsv.bgz``.",
                },
                {
                    "tag": "collect",
                    "text": (
                        "Evaluate ``loeuf_percentile_thresholds`` and parallelize one row per"
                        " (granularity, bin, threshold). Export"
                        " ``loeuf_percentile_thresholds.tsv``."
                    ),
                },
                {
                    "tag": "explode",
                    "text": (
                        "``--export-release-downsampling-tsv``: one row per (transcript, genetic"
                        " ancestry, downsampling), with obs and exp columns for ``syn``, ``mis``,"
                        " ``lof_hc_lc`` and ``lof``. Export the ``.downsampling.tsv.bgz``."
                    ),
                },
            ],
            "notes": [
                (
                    "``--export-release-downsampling-tsv`` needs a release built from a run"
                    " with ``--include-downsamplings``. Otherwise the ``downsamplings``"
                    " global is absent and the step fails."
                ),
                (
                    "In the thresholds TSV, ``bin`` k is the k-th boundary, counting from 1."
                    " A transcript whose LOEUF is at or above it is in bin k or higher, with"
                    " bins counted from 0."
                ),
                (
                    "The thresholds come from flag-free MANE Select transcripts, while"
                    " ``lof.oe_ci.upper_bin_*`` is rank-based over the ranked set. The two"
                    " need not agree for a given transcript."
                ),
            ],
        },
    },
    "functions": {
        "prepare_context_ht": {
            "notes": [
                (
                    "``adj_r`` and ``syn_adj_r`` dicts are indexed with ``context``, which is"
                    " already the strand-collapsed trimer here, so those Tables must be keyed"
                    " by collapsed contexts."
                ),
                (
                    "``sfs_bin`` uses the full-exome adj AF. Bin 0 holds both missing AF and"
                    " AF 0. AF above the last cutoff (0.5) is missing."
                ),
                ("Only the first (adj, all samples) element of each AN array is kept."),
            ],
        },
        "assemble_constraint_context_ht": {
            "notes": [
                (
                    "``methylation_level`` is 0-2 and set only for CpG variants; non-CpG is"
                    " 0. GRCh38 scores map with cutoffs 0 / 1-5 / 6+ on autosomes and PAR,"
                    " and 0 / 1-3 / 4+ on chrX and chrY non-PAR."
                ),
                (
                    "Missing GERP is set to 0, so with the default cutoffs those sites pass"
                    " the mutation rate GERP filter."
                ),
                (
                    "Frequencies and filters are joined on ``(locus, alleles)``; methylation,"
                    " GERP, coverage and AN on ``locus``."
                ),
            ],
        },
        "update_loftee_end_trunc_filter": {
            "notes": [
                (
                    "Recomputes ``END_TRUNC`` with a GERP distance cutoff of 0, then sets"
                    " ``lof`` to HC when no LOFTEE filter remains and LC otherwise. It can"
                    " move transcripts between the HC and LC groups."
                ),
            ],
        },
        "calculate_gerp_cutoffs": {
            "notes": [
                (
                    "The pipeline passes autosome and PAR rows of the annotated context"
                    " Table. Every possible SNV is a row, so each base counts three times,"
                    " and missing GERP is already 0."
                ),
            ],
        },
        "get_annotations_for_computing_mu": {
            "notes": [
                (
                    "``possible_variants`` comes from the full-genome adj entry (absent, or"
                    " AC > 0). ``observed_variants`` is AC > 0 per entry with no AF cap. The"
                    " AC <= 5 filter is applied to the global downsampling of 1,000, not the"
                    " full dataset."
                ),
                (
                    "Coverage is genome mean coverage, inclusive at both ends. GERP bounds"
                    " are strict. ``in_autosome`` excludes PAR."
                ),
            ],
        },
        "get_exome_coverage_expr": {
            "notes": [
                (
                    "The pipeline default is ``AN_percent`` (the function default is"
                    " ``median``): ``int(AN / possible alleles * 100)``. Possible alleles are"
                    " 2 x samples on autosomes and PAR, 2 x XX + XY on chrX non-PAR, and XY"
                    " on chrY non-PAR. Every coverage cutoff downstream is in these units."
                ),
            ],
        },
        "get_exomes_observed_and_possible": {
            "notes": [
                (
                    "With ``max_af``, ``observed_variants`` checks each entry's own AF, but"
                    " ``possible_variants`` checks only the full adj entry. A PASS variant"
                    " above ``max_af`` overall drops out of both counts."
                ),
            ],
        },
        "get_build_calibration_model_annotation": {
            "notes": [
                (
                    "A training site's transcript must have ``synonymous_variant`` as its own"
                    " most severe consequence and be an ``ENST`` MANE Select transcript (the"
                    " default; the pipeline never overrides it)."
                ),
                (
                    "A missing ``low_cov_cutoff`` becomes 0 here, so zero-coverage sites"
                    " never train. ``cpg`` is kept in the model group only for high-coverage"
                    " sites."
                ),
            ],
        },
        "calibration_model_group_expr": {
            "notes": [
                (
                    "High: ``high_cov_cutoff <= coverage <= upper_cov_cutoff``. Low:"
                    " ``low_cov_cutoff < coverage < high_cov_cutoff``, never with"
                    " ``skip_coverage_model``. Sites above ``upper_cov_cutoff`` get no group."
                ),
            ],
        },
        "calculate_mu_by_downsampling": {
            "notes": [
                (
                    "``mu[i] = 1.2e-8 * total_bases * (observed[i] / possible) /"
                    " sum(observed[i])`` with ``total_bases = sum(possible) // 3``. This"
                    " scales the rates so the three substitutions at a base sum to 1.2e-8 on"
                    " average."
                ),
                (
                    "``mu_snp`` is ``mu[calculate_mu_globals.downsampling_idx]`` (global,"
                    " 1,000 genomes) and is the only rate used downstream."
                ),
            ],
        },
        "create_training_set": {
            "notes": [
                (
                    "The grouping includes ``exomes_coverage``, so each distinct coverage"
                    " value is its own row. ``build_models`` collapses high-coverage rows"
                    " again."
                ),
            ],
        },
        "annotate_with_mu": {
            "notes": [
                (
                    "Collects the mutation rate Table into a dict and broadcasts it instead"
                    " of joining. A row whose ``(context, ref, alt, methylation_level)`` is"
                    " not in it fails with 'Missing mu'."
                ),
            ],
        },
        "build_models": {
            "notes": [
                (
                    "Plateau: for each ``{cpg, genomic_region}`` group, regress observed /"
                    " possible on ``[1, mu_snp]`` over high-coverage rows, once per exome"
                    " frequency entry. Each value is an ``[intercept, slope]`` pair."
                ),
                (
                    "Coverage: ``scale = sum(observed[0]) / sum(possible * mu_snp)`` over"
                    " high-coverage autosome and PAR rows. Then ``oe = sum(observed[0]) /"
                    " (scale * sum(possible * mu_snp))`` per low-coverage level, and an"
                    " unweighted regression of ``oe`` on ``[1, coverage]``."
                ),
                (
                    "``log10_coverage`` defaults to True here, but the pipeline passes"
                    " ``--use-logarithmic-coverage-model``, which is off by default."
                ),
            ],
        },
        "check_skip_coverage_model": {
            "notes": [
                (
                    "Checks only ``skip_coverage_model``. Coverage cutoffs, the log10 setting"
                    " and ``--use-weights`` are not compared between steps."
                ),
            ],
        },
        "_prepare_ht_for_apply_models": {
            "notes": [
                (
                    "Drops rows before exploding: ``apply_model`` missing (sites in neither"
                    " coverage model) or ``possible_variants == 0`` (exome AF above"
                    " ``max_af``)."
                ),
            ],
        },
        "create_per_variant_expected_ht": {
            "notes": [
                (
                    "At the end it drops the top-level copies of the ``calibrate_mu`` fields."
                    " Per-variant observed and possible counts survive only inside"
                    " ``calibrate_mu``."
                ),
                (
                    "Rows are not filtered to coding consequences, and several rows share"
                    " each (locus, alleles) key."
                ),
            ],
        },
        "apply_models": {
            "notes": [
                (
                    "``ppo = mu * slope + intercept`` (one per plateau model, so one per freq"
                    " entry). ``expected = ppo * possible * cov_corr``. ``mu = mu * possible"
                    " * cov_corr``."
                ),
                (
                    "``cov_corr`` is 0 at coverage 0, ``slope * cov + intercept`` (log10 cov"
                    " if requested) for low-coverage rows, and 1 otherwise. The pipeline"
                    " passes ``model_group_expr=apply_model``, so ``high_cov_definition`` is"
                    " unused."
                ),
            ],
        },
        "aggregate_per_variant_expected_ht": {
            "notes": [
                (
                    "Re-promotes ``calibrate_mu`` to get observed and possible counts,"
                    " filters to ``CSQ_CODING``, and checkpoints to temp files before and"
                    " after the ``group_by``."
                ),
                (
                    "``include_mu_annotations_in_grouping`` is never set by the pipeline."
                ),
            ],
        },
        "create_aggregated_expected_ht": {
            "notes": [
                (
                    "The first ``group_by`` uses ``partition_hint`` 100, which the pipeline"
                    " does not expose as a flag."
                ),
                (
                    "Multiplies ``mu_snp``, ``predicted_proportion_observed`` and"
                    " ``coverage_correction`` by ``possible_variants`` before the second sum,"
                    " so its output matches the per-variant route."
                ),
            ],
        },
        "aggregate_by_constraint_groups": {
            "notes": [
                (
                    "With the default groupings, the internal group names are ``syn``,"
                    " ``mis``, ``lof_classic``, ``lof_hc_lc`` and ``lof_hc``, in that order."
                    " The release drops ``lof_classic`` and renames ``lof_hc`` to ``lof``."
                ),
                (
                    "``no_variants`` looks only at the first (adj) observed count of each"
                    " group."
                ),
            ],
        },
        "_annotate_oe_ci_z": {
            "notes": [
                (
                    "The flag names are ``no_exp_{syn,mis,lof}`` and"
                    " ``outlier_{syn,mis,lof}``; the inline comment's ``z_raw_*`` is wrong."
                    " All three LoF groups share the ``lof`` suffix and thresholds."
                ),
            ],
        },
        "calculate_raw_z_score_sd": {
            "notes": [
                (
                    "With ``mirror_neg_raw_z``, the SD uses only negative raw z (more"
                    " variation than expected), mirrored. Constrained transcripts have"
                    ' positive raw z, so they never enter the SD. The "left tail" comment in'
                    " ``_compute_z_scores`` has the sign backwards."
                ),
            ],
        },
        "compute_pli": {
            "notes": [
                (
                    "The default expected values are the v2 ones (0.463, 0.089); the pipeline"
                    " passes 0.706 and 0.207 from ``--expectation-*``."
                ),
                (
                    "The mixture weights are fit over every row of the Table, with one"
                    " aggregation per EM iteration."
                ),
            ],
        },
        "oe_confidence_interval": {
            "notes": [
                (
                    "Gamma: shape ``obs + 1``, scale ``1 / exp``, 5% in each tail (a 90%"
                    " interval). Needs Hail >= 0.2.137."
                ),
                (
                    "Poisson searches oe in [0, 2) in steps of 0.001, so its upper bound is"
                    " capped at 1.999. Ranks, thresholds and the release use gamma."
                ),
            ],
        },
        "add_oe_upper_rank_and_bins": {
            "notes": [
                (
                    "Ranks 10 values (5 groups x 2 CIs), each with its own global sort and"
                    " ``.cache()`` writes to the tmp dir. Only the gamma rank of ``lof_hc``"
                    " is released."
                ),
            ],
        },
        "compute_constraint_percentile_bins": {
            "notes": [
                (
                    "The thresholds always use ENST MANE Select transcripts without"
                    " constraint flags (``mane_select_only=True``), whatever"
                    " ``use_mane_select_over_canonical`` is set to. The resulting"
                    " ``constraint_bins`` go to every transcript."
                ),
            ],
        },
        "compute_gene_quality_metrics": {
            "notes": [
                (
                    "``gene_flags`` is built with ``add_filters_expr``, so a missing metric"
                    " makes the whole set missing."
                ),
            ],
        },
        "add_gencode_transcript_annotations": {
            "notes": [
                (
                    "Inner-joins transcript features with summed CDS features, so transcripts"
                    " with no CDS get no GENCODE annotations."
                ),
            ],
        },
        "get_gencode_ht": {
            "notes": [
                (
                    "Imports the resource's default GENCODE version but always checkpoints to"
                    " ``gencode_v39_with_versions.ht``, read if present. A version change"
                    " would silently reuse the old checkpoint."
                ),
            ],
        },
        "_restructure_release_rows": {
            "notes": [
                (
                    "``RELEASE_CG_RENAME`` maps ``mu_snp`` to ``mu``, which replaces the"
                    " internal coverage-corrected ``mu``. The release ``mu`` is the summed"
                    " ``mu_snp``."
                ),
            ],
        },
        "_restructure_release_globals": {
            "notes": [
                (
                    "Of ``percentile_thresholds``, only ``lof`` is kept (as"
                    " ``loeuf_percentile_thresholds``). ``version`` is set only from"
                    " ``release_version`` or an existing global, which the pipeline never"
                    " sets."
                ),
            ],
        },
    },
    "resources": {
        "annotated_context_ht": {
            "summary": (
                "Every possible SNV on chr1-22/X/Y, strand-collapsed, with slimmed VEP,"
                " methylation, GERP, coverage, AN, exome and genome frequencies and"
                " filters, genomic region, SFS bin, and adj_r."
            ),
            "readme": "1-annotated_contextht",
        },
        "temp_preprocess_data_ht": {
            "summary": (
                "The annotated context plus ``compute_mu`` (genome counts for the"
                " mutation rate), ``calibrate_mu`` (exome counts and model buckets) and"
                " the parameter globals. Despite the ``temp_`` name it is written to the"
                " permanent ``preprocessed_data/`` path."
            ),
            "readme": "2-contextpreprocessedht",
        },
        "mutation_ht": {
            "summary": (
                "One row per context, ref, alt and methylation level, with summed genome"
                " counts, ``mu`` per frequency entry, and the scalar ``mu_snp`` used"
                " downstream."
            ),
            "readme": "3-mutation_rateht",
        },
        "train_ht": {
            "summary": (
                "Synonymous exome counts grouped by context, mutation type, region, model"
                " bucket and exome coverage, with ``mu_snp`` attached."
            ),
            "readme": "4-constraint_traininght",
        },
        "train_tsv": {
            "summary": "Flat TSV export of ``constraint_training.ht``.",
            "readme": "4-constraint_traininght",
        },
        "model_plateau": {
            "summary": (
                "Dict keyed by ``{cpg, genomic_region}``. Each value holds one"
                " ``[intercept, slope]`` pair per exome frequency entry. README §5 has"
                " the pair order reversed."
            ),
            "readme": "5-plateauhe-and-coveragehe",
        },
        "model_coverage": {
            "summary": (
                "``[intercept, slope]`` of the low-coverage correction. Linear in"
                " coverage unless ``--use-logarithmic-coverage-model`` was set; not"
                " written with ``--skip-coverage-model``."
            ),
            "readme": "5-plateauhe-and-coveragehe",
        },
        "per_variant_apply_ht": {
            "summary": (
                "Per-SNV expected counts, one row per (variant, transcript consequence),"
                " keyed by (locus, alleles) with repeated keys. Observed and possible"
                " counts are in ``calibrate_mu``."
            ),
            "readme": "6-per_variant_expectedht",
        },
        "apply_ht": {
            "summary": (
                "Per-variant counts summed per (annotation, modifier, gene, gene_id,"
                " transcript, canonical, mane_select), coding consequences only."
            ),
            "readme": "7-per_variant_expectedaggregatedht",
        },
        "aggregated_expected_ht": {
            "summary": (
                "The same sums, built by aggregating before applying the models. No"
                " per-variant rows."
            ),
            "readme": "8-aggregated_expectedht",
        },
        "constraint_group_ht": {
            "summary": (
                "Per-transcript sums for each constraint group (syn, mis, lof_classic,"
                " lof_hc_lc, lof_hc), with per-freq ``oe_info``."
            ),
            "readme": "9-constraint_groupht",
        },
        "gene_quality_metrics_ht": {
            "summary": (
                "Per-transcript CDS coverage, mapping quality, segdup and LCR overlap,"
                " plus ``gene_flags``. Keyed by transcript."
            ),
            "readme": "11-gene_quality_metricsht",
        },
        "pre_rank_constraint_metrics_ht": {
            "summary": (
                "Metrics before ranks and bins: oe, both CIs, z-scores, flags, pLI, gene"
                " quality and GENCODE. ``--skip-pre-rank-metrics`` reranks it."
            ),
        },
        "constraint_metrics_ht": {
            "summary": (
                "Internal per-transcript metrics with ranks, rank bins, threshold bins"
                " and ``percentile_thresholds``."
            ),
            "readme": "10-constraint_metricsht-internal",
        },
        "release_ht": {
            "summary": (
                "Public constraint metrics: one struct per released group (syn, mis,"
                " lof_hc_lc, lof) plus GENCODE, quality fields and flags."
            ),
            "readme": "12-releaseconstraint_metricsht",
        },
        "release_mutation_ht": {
            "summary": "Public per-context mutation rate (``mu``) with mutation type.",
            "readme": "13-releasemutation_rateht",
        },
        "release_mutation_tsv": {
            "summary": "TSV copy of the release mutation rate Table.",
            "readme": "13-releasemutation_rateht",
        },
        "release_tsv": {
            "summary": "Flattened release Table, bgzipped, with dotted column names.",
            "readme": "12-releaseconstraint_metricsht",
        },
        "release_downsampling_tsv": {
            "summary": (
                "Observed and expected counts per (transcript, genetic ancestry,"
                " downsampling) for syn, mis, lof_hc_lc and lof."
            ),
            "readme": "12-releaseconstraint_metricsht",
        },
        "lof_threshold_tsv": {
            "summary": "LOEUF bin boundaries, one row per (granularity, bin, threshold).",
            "readme": "12-releaseconstraint_metricsht",
        },
    },
    "flags": {
        "prepare_context_ht": {
            "notes": [
                (
                    "Also adds AN, exome and genome frequencies and filters,"
                    " ``genomic_region``, ``sfs_bin``, ``adj_r`` and ``syn_adj_r``, and"
                    " re-derives the LOFTEE ``END_TRUNC`` filter."
                ),
            ],
        },
        "calculate_gerp_cutoffs": {
            "notes": [
                (
                    "Bases with no GERP score count as 0 in the percentiles. The result is"
                    " recorded only in ``calculate_mu_globals`` of the preprocessed Table."
                ),
            ],
        },
        "calculate_mutation_rate_gerp_lower_cutoff": {
            "notes": [
                (
                    "This is the 5th percentile bound, not the 95th as the help says. Sites"
                    " must have GERP strictly above it."
                ),
            ],
        },
        "calculate_mutation_rate_gerp_upper_cutoff": {
            "notes": [
                (
                    "This is the 95th percentile bound, not the 5th as the help says. Sites"
                    " must have GERP strictly below it."
                ),
            ],
        },
        "calculate_mutation_rate_min_cov": {
            "notes": [
                (
                    "Genome mean coverage, inclusive, not the exome coverage metric. Applied"
                    " in ``--preprocess-data``."
                ),
            ],
        },
        "calculate_mutation_rate_max_cov": {
            "notes": [
                (
                    "Genome mean coverage, inclusive, not the exome coverage metric. Applied"
                    " in ``--preprocess-data``."
                ),
            ],
        },
        "exome_coverage_metric": {
            "notes": [
                (
                    "Sets the units of every coverage cutoff: high-coverage definitions,"
                    " upper cutoff, and low-coverage filter. With ``AN_percent`` they are"
                    " percent of possible alleles called."
                ),
            ],
        },
        "build_model_high_cov_definition": {
            "notes": [
                (
                    "Applied in ``--preprocess-data`` (the ``build_model`` bucket), not in"
                    " ``--build-models``. Inclusive (``>=``)."
                ),
            ],
        },
        "build_model_upper_cov_cutoff": {
            "notes": [
                (
                    "Sites above it are left out of the low-coverage model too, not only the"
                    " high-coverage Table. Applied in ``--preprocess-data``."
                ),
            ],
        },
        "apply_model_high_cov_definition": {
            "notes": [
                (
                    "Applied in ``--preprocess-data`` (the ``apply_model`` bucket). Rerun"
                    " from there before ``--apply-models-*`` to change it."
                ),
            ],
        },
        "pipeline_low_coverage_filter": {
            "notes": [
                (
                    "Works only through the ``build_model`` and ``apply_model`` buckets set"
                    " in ``--preprocess-data``. The bound is strict (``>``), so sites at the"
                    " cutoff are excluded too. When unset, training still excludes zero"
                    " coverage."
                ),
            ],
        },
        "max_af": {
            "notes": [
                (
                    "Applied in ``--preprocess-data`` to the exome counts only. Observed"
                    " checks each frequency entry's AF; possible checks the full adj AF. The"
                    " mutation rate uses AC <= 5 in the downsampling of 1,000 instead."
                ),
            ],
        },
        "genetic_ancestry_groups": {
            "notes": [
                (
                    "Takes effect in ``--preprocess-data``. It adds one exome frequency entry"
                    " per group (plus per-group downsamplings with"
                    " ``--include-downsamplings``) and one genome downsampling-1000 entry per"
                    " group for the mutation rate."
                ),
            ],
        },
        "include_downsamplings": {
            "notes": [
                (
                    "Exome counts only. Uses the reduced ``DOWNSAMPLINGS['v4']`` list unless"
                    " ``--genetic-ancestry-groups`` is set, in which case it uses every"
                    " downsampling in the exome frequency metadata."
                ),
            ],
        },
        "skip_coverage_model": {
            "notes": [
                (
                    "Pass the same value to ``--preprocess-data``, ``--build-models`` and the"
                    " apply steps. The later steps compare it with the stored globals and"
                    " raise on a mismatch."
                ),
                (
                    "With it, sites below ``--apply-model-high-cov-definition`` belong to no"
                    " model and drop out of both observed and expected counts."
                ),
            ],
        },
        "use_weights": {
            "notes": [
                (
                    "Weights only the plateau regressions. The coverage regression is always"
                    " unweighted."
                ),
            ],
        },
        "use_logarithmic_coverage_model": {
            "notes": [
                (
                    "Off by default, so the coverage model is linear in coverage, even though"
                    " ``build_models`` itself defaults to log10. Not checked between steps:"
                    " pass the same value to ``--build-models`` and ``--apply-models-*``."
                ),
            ],
        },
        "custom_vep_annotation": {
            "notes": [
                (
                    "Both apply steps hard-code ``use_mane_select=True``, so"
                    " ``worst_csq_by_gene`` raises a ValueError there."
                ),
            ],
        },
        "use_aggregated_expected": {
            "notes": [
                (
                    "Changes only the input. Both routes write the same"
                    " ``constraint_group.ht`` path."
                ),
            ],
        },
        "use_mane_select_over_canonical": {
            "notes": [
                (
                    "Affects only which transcripts are ranked. The percentile thresholds"
                    " always use ENST MANE Select transcripts without constraint flags."
                ),
            ],
        },
        "raw_z_outlier_threshold_lower_lof": {
            "notes": [
                (
                    "Applies to all three LoF groups. Raw z is positive for depletion, so a"
                    " lower threshold flags excess variation."
                ),
            ],
        },
        "raw_z_outlier_threshold_upper_syn": {
            "notes": [
                (
                    "Raw z is positive for depletion, so this flags strongly depleted"
                    " synonymous counts."
                ),
            ],
        },
        "release_version": {
            "notes": [
                (
                    "No pipeline step sets a ``version`` global, so without this flag the"
                    " release Tables have none."
                ),
            ],
        },
        "export_release_downsampling_tsv": {
            "notes": [
                (
                    "Needs a release from a run with ``--include-downsamplings``; otherwise"
                    " the ``downsamplings`` global is missing and the step fails."
                ),
            ],
        },
    },
    "downstream": [],
    "key_lanes_note": (
        "Every Table up to the per-variant expected Table is keyed by ``(locus,"
        " alleles)``. After that each group-by cuts the row count by orders of"
        " magnitude (6.3 billion per-variant rows, 3.3 million per consequence"
        " and transcript, 310,685 transcripts), and from ``constraint_group.ht``"
        " on the key is the transcript and never changes again."
    ),
}
