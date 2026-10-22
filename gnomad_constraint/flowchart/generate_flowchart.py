r"""Generate the data behind the interactive constraint pipeline flowchart.

The page, ``constraint_pipeline.html`` in this directory, is one self-contained
file. Three blocks inside it are generated, each between a pair of markers:

    - CODE: the pipeline steps, CLI flags, input and output resources, and the
      function call tree of every step, read from the code by ``--extract``.
      Needs neither a cluster nor GCS.
    - NOTES: the hand-written layer in ``notes.py``, validated against CODE by
      ``--extract`` and ``--check``.
    - TABLES: the schema, row count, globals and a short preview of every
      pipeline input and output, read on Dataproc by ``--collect`` and spliced
      in by ``--inject-tables``.

Usage::

    # Locally, from the repository root, with gnomad_methods and gnomad_qc
    # checkouts on PYTHONPATH.
    python -m gnomad_constraint.flowchart.generate_flowchart --extract

    # On Dataproc (absolute script path, single-zip pyfiles).
    hailctl dataproc submit <cluster> \
      <repo>/gnomad_constraint/flowchart/generate_flowchart.py \
      --pyfiles <zip> -- --collect --out gs://gnomad-tmp-4day/<path>/tables.json

    # Locally again.
    python -m gnomad_constraint.flowchart.generate_flowchart \
      --inject-tables tables.json
    python -m gnomad_constraint.flowchart.generate_flowchart --check

``--extract`` records the commit of both repositories so every source link on
the page is a permalink. It refuses to run when any file it reads has
uncommitted changes, unless ``--allow-dirty`` is passed, in which case the page
shows a warning banner.

Function call trees follow functions in gnomad_constraint all the way down and
stop at the first gnomad_methods function, which is shown as a leaf.
"""

import argparse
import ast
import copy
import datetime
import gzip
import importlib.util
import json
import logging
import pathlib
import re
import subprocess
import types
from typing import Any, Dict, List, Optional, Set, Tuple
from unittest import mock

import hail as hl
from gnomad.resources.resource_utils import BaseResource

import gnomad_constraint.pipeline.constraint_pipeline as constraint_pipeline
import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.flowchart.notes import NOTES
from gnomad_constraint.resources.constants import CURRENT_VERSION

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("generate_flowchart")
logger.setLevel(logging.INFO)

HTML_PATH = str(pathlib.Path(__file__).with_name("constraint_pipeline.html"))
"""The page the generated blocks are spliced into."""

PIPELINE_MODULE = "gnomad_constraint.pipeline.constraint_pipeline"
"""Module whose ``main()`` defines the pipeline steps."""

CUSTOM_VEP_ANNOTATION = "transcript_consequences"
"""VEP annotation the documented resource paths are built for (the CLI default)."""

BUILD = "GRCh38"
"""Reference build of every version the pipeline still supports (v4.x)."""

REPOS = {
    "gnomad_constraint": "gnomad_constraint.",
    "gnomad_methods": "gnomad.",
}
"""Repositories whose functions are documented, keyed by the module prefix."""

REPO_ANCHORS = {
    "gnomad_constraint": PIPELINE_MODULE,
    "gnomad_methods": "gnomad.utils.constraint",
}
"""A module in each repository, used to find its checkout (``gnomad`` itself is a
namespace package with no file of its own)."""

REPO_SLUGS = {
    "gnomad_constraint": "broadinstitute/gnomad-constraint",
    "gnomad_methods": "broadinstitute/gnomad_methods",
}
"""GitHub repository of each, used when a package is installed rather than checked
out (as in CI), so there is no git remote to read."""

RESOURCE_READ_ATTRS = {"ht", "he", "mt"}
"""Methods that read a resource object, e.g. ``seg_dup_intervals.ht()``."""

HAIL_OPS = {
    "write": "write",
    "checkpoint": "checkpoint",
    "export": "export",
    "repartition": "repartition",
    "naive_coalesce": "naive_coalesce",
    "write_expression": "write",
    "_set_flags": "shuffle_flag",
}
"""Hail method calls recorded on each step, mapped to the tag the page shows."""

GLOBALS_METHOD_PREFIXES = ("annotate_", "select_", "drop_", "transmute_")
"""Attribute prefixes that mark a method call rather than a read of a global."""

PREVIEW_ROWS = 5
"""Rows in each table preview."""

PREVIEW_FIELD_CAP = 40
"""Row fields carried in a preview; the schema block lists all of them."""

ARRAY_CAP = 12
"""Longer arrays of scalars are summarised rather than dumped."""

GLOBAL_STRUCT_ARRAY_CAP = 128
"""Arrays of structs in globals are inventories, so they get a far higher cap."""

PREVIEW_STRUCT_ARRAY_CAP = 4
"""Arrays of structs in preview rows (e.g. VEP transcript consequences)."""

STRING_CAP = 240
"""Strings longer than this are cut. VEP config globals run to kilobytes."""

TSV_PREVIEW_LINES = 6
"""Lines of each exported TSV shown on the page (header plus five rows)."""

MARKERS = {
    "CODE": ("/* CODE-DATA-START */", "/* CODE-DATA-END */"),
    "NOTES": ("/* NOTES-DATA-START */", "/* NOTES-DATA-END */"),
    "TABLES": ("/* TABLES-DATA-START */", "/* TABLES-DATA-END */"),
}
"""Start and end markers of each generated block in the page."""

VOLATILE_KEYS = {"line", "lines", "sha", "dirty", "generated", "branch"}
"""Keys ``--check`` ignores, so unrelated edits that shift lines do not fail it."""

UPSTREAM_DOC_KEYS = {"summary", "param_docs", "return_doc"}
"""Docstring keys ``--check`` ignores on gnomad_methods functions, so an upstream
documentation edit does not fail this repository's CI. Signatures still count."""


# ---------------------------------------------------------------------------
# Repositories
# ---------------------------------------------------------------------------


def _git(cwd: str, *args: str) -> str:
    """
    Run a git command and return its stripped stdout.

    :param cwd: Directory to run the command in.
    :param args: Arguments to git.
    :return: Stripped stdout.
    """
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


def get_repo_info(repo: str) -> Dict[str, Any]:
    """
    Describe where a documented repository's code was imported from.

    :param repo: Repository key from ``REPOS``.
    :return: Dictionary with the import ``root`` (the directory holding the
        top-level package), GitHub ``slug``, and the ``branch`` and ``sha`` of the
        checkout, which are None when the package is installed rather than
        checked out.
    """
    module = REPO_ANCHORS[repo]
    origin = pathlib.Path(importlib.util.find_spec(module).origin).resolve()
    root = str(origin.parents[len(module.split(".")) - 1])
    info = {"root": root, "slug": REPO_SLUGS[repo], "branch": None, "sha": None}
    try:
        toplevel = _git(root, "rev-parse", "--show-toplevel")
    except subprocess.CalledProcessError:
        return info
    if pathlib.Path(toplevel).resolve() != pathlib.Path(root):
        return info
    remote = _git(root, "remote", "get-url", "origin")
    info["slug"] = re.sub(r"(\.git)?$", "", re.split(r"github\.com[:/]", remote)[-1])
    info["branch"] = _git(root, "rev-parse", "--abbrev-ref", "HEAD")
    info["sha"] = _git(root, "rev-parse", "HEAD")
    return info


def get_dirty_files(root: str, files: List[str]) -> List[str]:
    """
    Return which of ``files`` differ from the commit checked out at ``root``.

    :param root: Repository root.
    :param files: Repository-relative paths.
    :return: The subset of ``files`` with uncommitted changes.
    """
    if not files:
        return []
    # Not _git: stripping the output would cut the status column of the first line.
    changed = subprocess.run(
        ["git", "status", "--porcelain", "--", *files],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return sorted(line[3:] for line in changed.splitlines())


# ---------------------------------------------------------------------------
# Source parsing
# ---------------------------------------------------------------------------


class ModuleSource:
    """
    Parsed source of one Python module: its functions and how it names imports.

    :param module: Dotted module name.
    """

    def __init__(self, module: str) -> None:
        """Parse ``module`` from the file it is imported from."""
        self.module = module
        self.file = importlib.util.find_spec(module).origin
        with open(self.file) as f:
            self.tree = ast.parse(f.read())

        # Names bound to modules (``import x.y as z``) and to symbols
        # (``from x import y``), mapped to the dotted path they stand for.
        self.module_aliases: Dict[str, str] = {}
        self.symbols: Dict[str, str] = {}
        self.functions: Dict[str, ast.FunctionDef] = {}
        for node in self.tree.body:
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.asname:
                        self.module_aliases[a.asname] = a.name
                    else:
                        self.module_aliases[a.name.split(".")[0]] = a.name.split(".")[0]
            elif isinstance(node, ast.ImportFrom) and node.module:
                for a in node.names:
                    self.symbols[a.asname or a.name] = f"{node.module}.{a.name}"
            elif isinstance(node, ast.FunctionDef):
                self.functions[node.name] = node

    def resolve(self, func: ast.expr) -> Tuple[Optional[str], Optional[str]]:
        """
        Resolve the callee of a call to a dotted path.

        :param func: The ``func`` of an ``ast.Call``.
        :return: Tuple of (kind, dotted path). Kind is "call" for a function call,
            "read" for a read of a resource object (``x.ht()``), or None when the
            callee is not a module-level name.
        """
        if isinstance(func, ast.Name):
            if func.id in self.functions:
                return "call", f"{self.module}.{func.id}"
            if func.id in self.symbols:
                return "call", self.symbols[func.id]
        elif isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
            name = func.value.id
            if name in self.module_aliases:
                return "call", f"{self.module_aliases[name]}.{func.attr}"
            if name in self.symbols and func.attr in RESOURCE_READ_ATTRS:
                return "read", self.symbols[name]
        return None, None


_SOURCES: Dict[str, ModuleSource] = {}


def get_source(module: str) -> Optional[ModuleSource]:
    """
    Return the parsed source of ``module``, caching it.

    :param module: Dotted module name.
    :return: Parsed module, or None if it cannot be found.
    """
    if module not in _SOURCES:
        try:
            _SOURCES[module] = ModuleSource(module)
        except (ImportError, AttributeError, TypeError, ValueError):
            return None
    return _SOURCES[module]


def get_repo(path: str) -> Optional[str]:
    """
    Return which documented repository a dotted path belongs to.

    :param path: Dotted path of a module or function.
    :return: Repository key from ``REPOS``, or None.
    """
    if path.startswith("gnomad_constraint.flowchart"):
        return None
    for repo, prefix in REPOS.items():
        if path.startswith(prefix):
            return repo
    return None


def find_function(path: str) -> Tuple[Optional[ModuleSource], Optional[ast.AST]]:
    """
    Find the definition of a function from its dotted path.

    :param path: Dotted path, e.g. ``gnomad.utils.constraint.apply_models``.
    :return: Tuple of the parsed module and the ``FunctionDef``, or (None, None)
        when the path is not a function defined at module level.
    """
    module, _, name = path.rpartition(".")
    src = get_source(module)
    if src is None or name not in src.functions:
        return None, None
    return src, src.functions[name]


def parse_docstring(doc: Optional[str]) -> Dict[str, Any]:
    """
    Split a Sphinx-style docstring into its summary, parameters, and return.

    :param doc: Cleaned docstring, as returned by ``ast.get_docstring``.
    :return: Dictionary with ``summary``, ``params`` (name to description), and
        ``returns``.
    """
    summary: List[str] = []
    params: Dict[str, str] = {}
    returns = ""
    # Which field a continuation line belongs to: "summary", "return", "skip",
    # or the name of the parameter being described.
    section = "summary"
    for line in (doc or "").splitlines():
        stripped = line.strip()
        param = re.match(r":param (\w+):\s*(.*)", stripped)
        if param:
            section = param.group(1)
            params[section] = param.group(2)
        elif re.match(r":returns?:", stripped):
            section = "return"
            returns = stripped.split(":", 2)[2].strip()
        elif re.match(r":(raises|rtype|type)\b", stripped):
            section = "skip"
        elif section == "summary":
            summary.append(line)
        elif stripped and section == "return":
            returns += " " + stripped
        elif stripped and section != "skip":
            params[section] += " " + stripped
    return {"summary": "\n".join(summary).strip(), "params": params, "returns": returns}


def describe_signature(node: ast.FunctionDef) -> List[Dict[str, Optional[str]]]:
    """
    List a function's parameters with their annotations and defaults.

    :param node: Function definition.
    :return: One dictionary per parameter with ``name``, ``annotation``, and
        ``default`` (source text, or None).
    """
    a = node.args
    positional = a.posonlyargs + a.args
    defaults: List[Optional[ast.expr]] = [None] * (
        len(positional) - len(a.defaults)
    ) + list(a.defaults)
    params = []
    for arg, default in list(zip(positional, defaults)) + list(
        zip(a.kwonlyargs, a.kw_defaults)
    ):
        params.append(
            {
                "name": arg.arg,
                "annotation": ast.unparse(arg.annotation) if arg.annotation else None,
                "default": ast.unparse(default) if default is not None else None,
            }
        )
    for prefix, arg in (("*", a.vararg), ("**", a.kwarg)):
        if arg:
            params.append(
                {"name": prefix + arg.arg, "annotation": None, "default": None}
            )
    return params


def get_globals_written(node: ast.AST) -> List[str]:
    """
    List the globals a function sets with ``annotate_globals``/``select_globals``.

    :param node: Function definition.
    :return: Sorted global names.
    """
    names: Set[str] = set()
    for n in ast.walk(node):
        if (
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute)
            and n.func.attr in {"annotate_globals", "select_globals"}
        ):
            names.update(k.arg for k in n.keywords if k.arg)
    return sorted(names)


def get_globals_read(node: ast.AST) -> List[str]:
    """
    List the ``*_globals`` fields a block of code reads.

    :param node: Any AST node.
    :return: Sorted global names.
    """
    return sorted(
        {
            n.attr
            for n in ast.walk(node)
            if isinstance(n, ast.Attribute)
            and n.attr.endswith("_globals")
            and not n.attr.startswith(GLOBALS_METHOD_PREFIXES)
        }
    )


def get_calls(
    node: ast.AST, src: ModuleSource
) -> Tuple[List[Tuple[str, ast.Call]], List[Tuple[str, int]]]:
    """
    List the documented functions and resource objects a block of code uses.

    :param node: Any AST node.
    :param src: Parsed module the node belongs to.
    :return: Tuple of (calls, reads). Calls are (dotted path, call node) pairs in
        source order, one per call site, limited to functions defined in a
        documented repository. Reads are (dotted path, line) pairs of resource
        objects read with ``.ht()``/``.he()``/``.mt()``.
    """
    calls, reads = [], []
    for n in sorted(
        (n for n in ast.walk(node) if isinstance(n, ast.Call)),
        key=lambda n: (n.lineno, n.col_offset),
    ):
        kind, path = src.resolve(n.func)
        if kind is None or get_repo(path) is None:
            continue
        if kind == "read":
            reads.append((path, n.lineno))
        elif find_function(path)[1] is not None:
            calls.append((path, n))
    return calls, reads


def describe_function(path: str, expand: bool) -> Dict[str, Any]:
    """
    Describe one function for the page.

    :param path: Dotted path of the function.
    :param expand: Whether to record the functions it calls. False for
        gnomad_methods functions, which are leaves of the call tree.
    :return: The function's record.
    """
    src, node = find_function(path)
    repo = get_repo(path)
    doc = parse_docstring(ast.get_docstring(node))
    record = {
        "name": node.name,
        "module": src.module,
        "repo": repo,
        "file": src.file,
        "lines": [node.lineno, node.end_lineno],
        "params": describe_signature(node),
        "return_annotation": ast.unparse(node.returns) if node.returns else None,
        "summary": doc["summary"],
        "param_docs": doc["params"],
        "return_doc": doc["returns"],
        "globals_written": get_globals_written(node),
        "globals_read": get_globals_read(node),
        "calls": None,
        "reads": [],
    }
    if expand:
        calls, reads = get_calls(node, src)
        record["calls"] = list(dict.fromkeys(p for p, _ in calls if p != path))
        record["reads"] = list(dict.fromkeys(p for p, _ in reads))
    return record


def collect_functions(roots: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Describe every function reachable from ``roots``.

    gnomad_constraint functions are followed through all their calls;
    gnomad_methods functions are recorded but not followed.

    :param roots: Dotted paths of the functions the pipeline steps call.
    :return: Function records keyed by dotted path.
    """
    functions: Dict[str, Dict[str, Any]] = {}
    queue = list(roots)
    while queue:
        path = queue.pop(0)
        if path in functions:
            continue
        expand = get_repo(path) == "gnomad_constraint"
        functions[path] = describe_function(path, expand)
        queue.extend(functions[path]["calls"] or [])
    return functions


# ---------------------------------------------------------------------------
# Pipeline structure
# ---------------------------------------------------------------------------


class _NameEvents(ast.NodeVisitor):
    """Record loads and stores of plain names in evaluation order."""

    def __init__(self) -> None:
        """Start with no events."""
        self.events: List[Tuple[str, str]] = []

    def visit_Assign(self, node: ast.Assign) -> None:
        """Visit the value before the targets it is stored in."""
        self.visit(node.value)
        for target in node.targets:
            self.visit(target)

    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        """Visit the value before the target it is stored in."""
        if node.value is not None:
            self.visit(node.value)
        self.visit(node.target)

    def visit_Name(self, node: ast.Name) -> None:
        """Record the name and whether it is loaded or stored."""
        kind = "store" if isinstance(node.ctx, ast.Store) else "load"
        self.events.append((node.id, kind))


def name_events(nodes: List[ast.AST]) -> Tuple[Set[str], Set[str]]:
    """
    Split the names used in a block into those it reads first and those it sets.

    :param nodes: Statements of the block.
    :return: Tuple of (names whose first use is a load, names stored anywhere).
    """
    visitor = _NameEvents()
    for n in nodes:
        visitor.visit(n)
    first: Dict[str, str] = {}
    stored = set()
    for name, kind in visitor.events:
        first.setdefault(name, kind)
        if kind == "store":
            stored.add(name)
    return {n for n, k in first.items() if k == "load"}, stored


def get_flag_refs(node: ast.AST, aliases: Dict[str, Set[str]]) -> Set[str]:
    """
    Return the argparse destinations an expression depends on.

    :param node: Any AST node.
    :param aliases: Local names in ``main()`` bound to argparse destinations.
    :return: Destinations read directly (``args.x``) or through an alias.
    """
    refs = set()
    for n in ast.walk(node):
        if (
            isinstance(n, ast.Attribute)
            and isinstance(n.value, ast.Name)
            and n.value.id == "args"
        ):
            refs.add(n.attr)
        elif isinstance(n, ast.Name) and n.id in aliases:
            refs.update(aliases[n.id])
    return refs


def get_main_structure(src: ModuleSource) -> Dict[str, Any]:
    """
    Read the step blocks out of the pipeline's ``main()``.

    Every ``if args.<flag>:`` statement inside the ``try`` block of ``main()`` is
    one pipeline step.

    :param src: Parsed pipeline module.
    :return: Dictionary with ``blocks`` (one per step, in source order),
        ``aliases`` (local names bound to argparse destinations), and
        ``resource_kwargs`` (keyword arguments of the ``get_constraint_resources``
        call that come from flags).
    """
    main = src.functions["main"]
    aliases: Dict[str, Set[str]] = {}
    resource_kwargs: Dict[str, Set[str]] = {}
    try_node = None
    for stmt in main.body:
        if isinstance(stmt, ast.Assign):
            _, path = (
                src.resolve(stmt.value.func)
                if isinstance(stmt.value, ast.Call)
                else (None, None)
            )
            if path and path.endswith(".get_constraint_resources"):
                # Not an alias: every step reads ``resources``, and the flags it
                # is built from are reported once, as ``resource_flags``.
                for i, arg in enumerate(stmt.value.args):
                    resource_kwargs[f"#{i}"] = get_flag_refs(arg, aliases)
                for kw in stmt.value.keywords:
                    resource_kwargs[kw.arg] = get_flag_refs(kw.value, aliases)
                continue
            refs = get_flag_refs(stmt.value, aliases)
            for target in stmt.targets:
                if isinstance(target, ast.Name) and refs:
                    aliases[target.id] = refs
        elif isinstance(stmt, ast.Try):
            try_node = stmt

    blocks = []
    for stmt in try_node.body:
        if not isinstance(stmt, ast.If):
            continue
        triggers = sorted(get_flag_refs(stmt.test, {}))
        step_ids = [
            n.attr
            for n in ast.walk(stmt)
            if isinstance(n, ast.Attribute)
            and isinstance(n.value, ast.Name)
            and n.value.id == "resources"
        ]
        loads, stores = name_events(stmt.body)
        blocks.append(
            {
                "node": stmt,
                "step_id": step_ids[0],
                "triggers": triggers,
                "loads": loads,
                "stores": stores,
            }
        )
    return {"blocks": blocks, "aliases": aliases, "resource_kwargs": resource_kwargs}


def get_hail_ops(node: ast.AST) -> List[Dict[str, Any]]:
    """
    List the Hail writes, exports, repartitions, and shuffle flags in a block.

    :param node: Any AST node.
    :return: One dictionary per call with ``op`` and ``line``.
    """
    ops = []
    for n in ast.walk(node):
        if (
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute)
            and n.func.attr in HAIL_OPS
        ):
            op = HAIL_OPS[n.func.attr]
            if op == "shuffle_flag":
                values = [k.value for k in n.keywords]
                if not values or not (
                    isinstance(values[0], ast.Constant) and values[0].value
                ):
                    continue
            ops.append({"op": op, "line": n.lineno})
    return sorted(ops, key=lambda o: o["line"])


def get_call_arguments(
    call: ast.Call, path: str, aliases: Dict[str, Set[str]]
) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """
    Map the arguments of a call to the flags they come from, or to literals.

    :param call: The call node.
    :param path: Dotted path of the called function.
    :param aliases: Local names in ``main()`` bound to argparse destinations.
    :return: Tuple of (flags, fixed). Flags maps each parameter whose argument
        depends on a flag to the sorted argparse destinations. Fixed maps each
        parameter given a literal to its source text.
    """
    _, node = find_function(path)
    names = [p["name"] for p in describe_signature(node)]
    args = [(names[i] if i < len(names) else None, a) for i, a in enumerate(call.args)]
    args += [(kw.arg, kw.value) for kw in call.keywords]
    flags, fixed = {}, {}
    for name, value in args:
        if name is None:
            continue
        refs = get_flag_refs(value, aliases)
        if refs:
            flags[name] = sorted(refs)
        elif not any(isinstance(n, ast.Name) for n in ast.walk(value)) and isinstance(
            value, (ast.Constant, ast.List, ast.Tuple, ast.Dict)
        ):
            fixed[name] = ast.unparse(value)
    return flags, fixed


def describe_resource(name: str, resource: Any) -> Dict[str, str]:
    """
    Describe one pipeline resource.

    :param name: Attribute name the pipeline uses for it.
    :param resource: TableResource, ExpressionResource, other resource, or a path.
    :return: Dictionary with ``name``, ``path``, and ``kind`` (``ht``, ``he``, or
        ``file``).
    """
    path = resource if isinstance(resource, str) else resource.path
    if path.endswith(".ht"):
        kind = "ht"
    elif path.endswith(".he"):
        kind = "he"
    else:
        kind = "file"
    return {"name": name, "path": path, "kind": kind}


def load_pipeline_steps(**kwargs: Any) -> Dict[str, Any]:
    """
    Build the pipeline's resource collection without reading any table.

    ``get_constraint_resources`` reads the VEP context Table only to get its
    reference build. Every supported version is GRCh38, so that read is replaced
    with ``BUILD``, which lets this run without GCS access.

    :param kwargs: Keyword arguments to ``get_constraint_resources`` that differ
        from the defaults.
    :return: The collection's steps keyed by step name.
    """
    real_get_vep_context_ht = constraint_res.get_vep_context_ht

    def _get_vep_context_ht(version: str) -> BaseResource:
        res = copy.copy(real_get_vep_context_ht(version))
        res.ht = mock.MagicMock()
        return res

    params = {
        "version": CURRENT_VERSION,
        "custom_vep_annotation": CUSTOM_VEP_ANNOTATION,
        "overwrite": False,
        "test": False,
        **kwargs,
    }
    with mock.patch.object(
        constraint_res, "get_vep_context_ht", _get_vep_context_ht
    ), mock.patch.object(
        constraint_res,
        "get_reference_genome",
        return_value=types.SimpleNamespace(name=BUILD),
    ):
        return constraint_res.get_constraint_resources(**params).pipeline_steps


def describe_step_resources(steps: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Describe the declared inputs, outputs, and upstream steps of every step.

    :param steps: Steps of a resource collection, keyed by step name.
    :return: Step name to a dictionary with ``name`` (the collection's label),
        ``deps``, ``inputs``, and ``outputs``.
    """
    by_label = {s.pipeline_step: k for k, s in steps.items()}
    out = {}
    for key, step in steps.items():
        # Inputs from previous steps are named in that step's outputs; inputs
        # added directly are named only by the attribute set on this step.
        names = {
            id(r): n
            for prev in step.previous_steps or []
            for n, r in prev._output_resource_dict.items()
        }
        names.update(
            {
                id(r): n
                for n, r in vars(step).items()
                if (hasattr(r, "path") or isinstance(r, str))
                and n not in {"pipeline_step", "pipeline_name"}
            }
        )
        inputs = []
        for source, resources in step.input_resources.items():
            for r in resources:
                res = describe_resource(names.get(id(r), "?"), r)
                res["from"] = by_label.get(source, source)
                inputs.append(res)
        out[key] = {
            "name": step.pipeline_step,
            "deps": sorted(
                {by_label[p.pipeline_step] for p in step.previous_steps or []}
            ),
            "inputs": inputs,
            "outputs": [
                describe_resource(n, r) for n, r in step._output_resource_dict.items()
            ],
        }
    return out


def get_flags() -> List[Dict[str, Any]]:
    """
    Describe every argument of the pipeline's parser.

    :return: One dictionary per argument, in parser order.
    """
    parser = constraint_pipeline.get_script_argument_parser()
    groups: Dict[int, List[str]] = {}
    for group in parser._action_groups:
        for action in group._group_actions:
            groups.setdefault(id(action), []).append(group.title)
    flags = []
    for action in parser._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        switch = action.nargs == 0 or isinstance(action, argparse.BooleanOptionalAction)
        flags.append(
            {
                "flag": action.option_strings[0],
                "aliases": action.option_strings[1:],
                "dest": action.dest,
                "groups": groups.get(id(action), []),
                "help": action.help,
                "default": action.default,
                "type": None if action.type is None else action.type.__name__,
                "choices": list(action.choices) if action.choices else None,
                "nargs": action.nargs,
                "switch": switch,
            }
        )
    return flags


def get_resource_modes(
    resource_kwargs: Dict[str, Set[str]], flags: Dict[str, Dict[str, Any]]
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    List the alternative resource layouts the pipeline can be run with.

    Each boolean keyword of ``get_constraint_resources`` that ``main()`` sets from
    a flag gives one alternative. ``--skip-coverage-model`` reaches the call
    through ``models`` instead, so it is added explicitly.

    :param resource_kwargs: Keyword arguments of the ``get_constraint_resources``
        call in ``main()`` mapped to the flags they come from.
    :param flags: Flag records keyed by argparse destination.
    :return: (flag, keyword arguments) pairs.
    """
    modes = [("--skip-coverage-model", {"models": ["plateau"]})]
    for kwarg, dests in resource_kwargs.items():
        if kwarg.startswith("#") or len(dests) != 1:
            continue
        dest = next(iter(dests))
        if flags.get(dest, {}).get("switch") and flags[dest]["default"] is False:
            modes.append((flags[dest]["flag"], {kwarg: True}))
    return modes


def extract(allow_dirty: bool = False) -> Dict[str, Any]:
    """
    Read the pipeline's structure out of the code.

    :param allow_dirty: Whether to continue when a file that is read has
        uncommitted changes.
    :return: The CODE block of the page.
    """
    src = get_source(PIPELINE_MODULE)
    structure = get_main_structure(src)
    aliases = structure["aliases"]
    flags = get_flags()
    flags_by_dest = {f["dest"]: f for f in flags}

    default_steps = describe_step_resources(load_pipeline_steps())
    modes = {}
    for flag, kwargs in get_resource_modes(structure["resource_kwargs"], flags_by_dest):
        for key, step in describe_step_resources(load_pipeline_steps(**kwargs)).items():
            if step != default_steps[key]:
                modes.setdefault(key, []).append({"flag": flag, **step})

    steps = []
    roots = []
    for i, block in enumerate(structure["blocks"]):
        node = block["node"]
        calls, reads = get_calls(node, src)
        call_records = []
        for path, call in calls:
            flags_in, fixed = get_call_arguments(call, path, aliases)
            call_records.append(
                {"fn": path, "line": call.lineno, "flags": flags_in, "fixed": fixed}
            )
            roots.append(path)

        # Values a block reads that an earlier block set, e.g. the GERP cutoffs.
        values_in = {}
        for earlier in structure["blocks"][:i]:
            shared = sorted(block["loads"] & earlier["stores"])
            if shared:
                values_in[earlier["step_id"]] = shared

        refs = get_flag_refs(node, aliases)
        steps.append(
            {
                "id": block["step_id"],
                "triggers": block["triggers"],
                "flags": sorted(refs),
                "file": src.file,
                "lines": [node.lineno, node.end_lineno],
                "calls": call_records,
                "reads": [{"path": p, "line": n} for p, n in reads],
                "ops": get_hail_ops(node),
                "globals_read": get_globals_read(node),
                "values_in": [{"from": k, "names": v} for k, v in values_in.items()],
                **default_steps[block["step_id"]],
                "modes": modes.get(block["step_id"], []),
            }
        )

    missing = set(default_steps) - {s["id"] for s in steps}
    if missing:
        raise ValueError(f"Resource steps with no block in main(): {sorted(missing)}")

    functions = collect_functions(list(dict.fromkeys(roots)))

    repos = {}
    for repo in REPOS:
        info = get_repo_info(repo)
        files = {
            str(pathlib.Path(s.file).resolve().relative_to(info["root"]))
            for s in _SOURCES.values()
            if get_repo(f"{s.module}.") == repo
        }
        info["dirty"] = (
            get_dirty_files(info["root"], sorted(files)) if info["sha"] else []
        )
        if info["dirty"] and not allow_dirty:
            raise SystemExit(
                f"Uncommitted changes in {info['root']}: {info['dirty']}. Commit them"
                " or pass --allow-dirty."
            )
        repos[repo] = info

    # Store repository-relative paths so links can be built from slug and sha.
    roots_by_repo = {r: i.pop("root") for r, i in repos.items()}
    for record in list(functions.values()) + steps:
        repo = record.get("repo", "gnomad_constraint")
        record["file"] = str(
            pathlib.Path(record["file"]).resolve().relative_to(roots_by_repo[repo])
        )

    globals_set = {
        g: sorted(p for p, f in functions.items() if g in f["globals_written"])
        for g in sorted({g for f in functions.values() for g in f["globals_written"]})
    }

    return {
        "meta": {
            "generated": datetime.date.today().isoformat(),
            "version": CURRENT_VERSION,
            "custom_vep_annotation": CUSTOM_VEP_ANNOTATION,
            "build": BUILD,
            "repos": repos,
            "resource_flags": sorted(
                {d for dests in structure["resource_kwargs"].values() for d in dests}
            ),
        },
        "flags": flags,
        "steps": steps,
        "functions": functions,
        "globals": globals_set,
    }


# ---------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------


def resolve_function_key(key: str, functions: Dict[str, Any]) -> Optional[str]:
    """
    Resolve a notes key to a function's dotted path.

    :param key: Dotted path, or a bare name that is unique among the functions.
    :param functions: Function records keyed by dotted path.
    :return: Dotted path, or None if the key matches no function or several.
    """
    if key in functions:
        return key
    matches = [p for p, f in functions.items() if f["name"] == key]
    return matches[0] if len(matches) == 1 else None


def validate_notes(code: Dict[str, Any]) -> List[str]:
    """
    Check that every key in the notes still names something in the code.

    :param code: The CODE block.
    :return: Problems found; empty when the notes are consistent.
    """
    problems = []
    step_ids = {s["id"] for s in code["steps"]}
    staged = [s for stage in NOTES["stages"] for s in stage.get("steps", [])]
    for step in sorted(step_ids - set(staged)):
        problems.append(f"step {step} is in no stage")
    for step in sorted(set(staged) - step_ids):
        problems.append(f"stage lists unknown step {step}")
    for step in sorted(set(NOTES["steps"]) - step_ids):
        problems.append(f"notes for unknown step {step}")
    for key in NOTES["functions"]:
        if resolve_function_key(key, code["functions"]) is None:
            problems.append(f"notes for unknown or ambiguous function {key}")
    dests = {f["dest"] for f in code["flags"]}
    for key in NOTES.get("flags", {}):
        if key not in dests:
            problems.append(f"notes for unknown flag {key}")
    outputs = {o["name"] for s in code["steps"] for o in s["outputs"] + s["inputs"]} | {
        o["name"] for s in code["steps"] for m in s["modes"] for o in m["outputs"]
    }
    for key in NOTES.get("resources", {}):
        if key not in outputs:
            problems.append(f"notes for unknown resource {key}")
    grouped = [m for g in NOTES["source_groups"] for m in g["members"]]
    sources = get_source_keys(code)
    for key in sorted(sources - set(grouped)):
        problems.append(f"source {key} is in no source group")
    for key in sorted(set(grouped) - sources):
        problems.append(f"source group lists unknown source {key}")
    return problems


def get_source_keys(code: Dict[str, Any]) -> Set[str]:
    """
    Name every input the pipeline reads but does not produce.

    Matches how the page builds its Sources band: declared inputs from outside the
    pipeline, resource getters called in a step block, and resource objects read
    with ``.ht()``.

    :param code: The CODE block.
    :return: Resource names, getter function names, and the last part of each
        resource object's dotted path.
    """
    step_ids = {s["id"] for s in code["steps"]}
    produced = {
        o["path"]
        for s in code["steps"]
        for o in s["outputs"] + [o for m in s["modes"] for o in m["outputs"]]
    }
    getter_prefix = f"{constraint_res.__name__}."
    keys = set()
    for step in code["steps"]:
        keys.update(
            i["name"]
            for i in step["inputs"]
            if i["from"] not in step_ids and i["path"] not in produced
        )
        for call in step["calls"]:
            name = call["fn"][len(getter_prefix) :]
            if call["fn"].startswith(getter_prefix) and re.fullmatch(
                r"get_.*_ht", name
            ):
                keys.add(name)
        keys.update(r["path"].split(".")[-1] for r in step["reads"])
    return keys


def resolved_notes(code: Dict[str, Any]) -> Dict[str, Any]:
    """
    Return the notes with function keys replaced by dotted paths.

    :param code: The CODE block.
    :return: Notes ready to embed.
    """
    notes = copy.deepcopy(NOTES)
    notes["functions"] = {
        resolve_function_key(k, code["functions"]): v
        for k, v in NOTES["functions"].items()
    }
    return notes


# ---------------------------------------------------------------------------
# Table inspection (Dataproc)
# ---------------------------------------------------------------------------


def jsonable(value: Any, struct_cap: int = GLOBAL_STRUCT_ARRAY_CAP) -> Any:
    """
    Convert a collected Hail value to something ``json.dumps`` accepts.

    Long arrays, dictionaries, and strings are summarised rather than dropped.

    :param value: A collected value.
    :param struct_cap: Longest array of structs kept whole.
    :return: JSON-safe structure.
    """
    if value is None or isinstance(value, (bool, int)):
        return value
    if isinstance(value, float):
        return value if value == value and abs(value) != float("inf") else str(value)
    if isinstance(value, str):
        return value if len(value) <= STRING_CAP else value[:STRING_CAP] + "…"
    if isinstance(value, hl.Struct):
        return {k: jsonable(v, struct_cap) for k, v in value.items()}
    if isinstance(value, hl.genetics.Locus):
        return f"{value.contig}:{value.position}"
    if isinstance(value, hl.Interval):
        return str(value)
    if isinstance(value, (dict, hl.utils.frozendict)):
        # Struct keys (e.g. the plateau model groups) become compact JSON text.
        items = [
            (k if isinstance(k, str) else json.dumps(jsonable(k)), v)
            for k, v in value.items()
        ]
        if len(items) > max(ARRAY_CAP, struct_cap):
            return {
                "__truncated__": len(items),
                "head": [f"{k}: {jsonable(v, struct_cap)}" for k, v in items[:6]],
                "tail": [f"{k}: {jsonable(v, struct_cap)}" for k, v in items[-3:]],
            }
        return {"__dict__": [[k, jsonable(v, struct_cap)] for k, v in items]}
    if isinstance(value, (set, frozenset)):
        value = sorted(value, key=str)
    if isinstance(value, (list, tuple)):
        structured = bool(value) and isinstance(
            value[0], (hl.Struct, dict, hl.utils.frozendict)
        )
        cap = struct_cap if structured else ARRAY_CAP
        if len(value) > cap:
            return {
                "__truncated__": len(value),
                "head": [jsonable(v, struct_cap) for v in value[:cap]],
                "tail": [jsonable(v, struct_cap) for v in value[-2:]],
            }
        return [jsonable(v, struct_cap) for v in value]
    return str(value)


def read_metadata(path: str) -> Dict[str, Any]:
    """
    Read a Hail Table's ``metadata.json.gz``.

    :param path: Table path.
    :return: Parsed metadata.
    """
    with hl.hadoop_open(f"{path}/metadata.json.gz", "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw.decode("utf-8"))


def modification_time(path: str) -> Optional[str]:
    """
    Return when a file was last written.

    :param path: File path.
    :return: UTC modification time in ISO format, or None.
    """
    try:
        seconds = float(hl.hadoop_stat(path)["modification_time"])
    except Exception:  # pylint: disable=broad-except
        return None
    return datetime.datetime.fromtimestamp(seconds, tz=datetime.timezone.utc).isoformat(
        timespec="minutes"
    )


def inspect_table(path: str) -> Dict[str, Any]:
    """
    Read one Table's schema, row count, globals, and first rows.

    Row counts come from the per-partition counts Hail stores at write time, so
    nothing is scanned.

    :param path: Table path.
    :return: The page's record for the Table.
    """
    meta = read_metadata(path)
    counts = (meta.get("components") or {}).get("partition_counts") or {}
    counts = counts.get("counts")
    ht = hl.read_table(path)
    row_fields = list(ht.row)
    key_fields = list(ht.key)
    preview_fields = row_fields
    if len(row_fields) > PREVIEW_FIELD_CAP:
        rest = [f for f in row_fields if f not in key_fields]
        preview_fields = key_fields + rest[: PREVIEW_FIELD_CAP - len(key_fields)]
    head = (
        ht.head(PREVIEW_ROWS)
        .select(*[f for f in preview_fields if f not in key_fields])
        .collect()
    )
    lines: List[str] = []
    ht.describe(handler=lines.append)
    glob = ht.globals.collect()[0]
    return {
        "kind": "ht",
        "written": modification_time(f"{path}/metadata.json.gz"),
        "hail_version": meta.get("hail_version"),
        "rows": sum(counts) if counts else None,
        "partitions": len(counts) if counts else None,
        "key": key_fields,
        "n_fields": len(row_fields),
        "describe": "\n".join(lines),
        "preview_fields": preview_fields,
        "preview": [
            {f: jsonable(r[f], PREVIEW_STRUCT_ARRAY_CAP) for f in preview_fields}
            for r in head
        ],
        "globals": [
            {"name": k, "type": str(ht[k].dtype), "value": jsonable(v)}
            for k, v in glob.items()
        ],
    }


def inspect_expression(path: str) -> Dict[str, Any]:
    """
    Read a written expression (the models) in full.

    :param path: Expression path.
    :return: The page's record for the expression.
    """
    expr = hl.experimental.read_expression(path)
    return {
        "kind": "he",
        "written": modification_time(f"{path}/metadata.json.gz"),
        "type": str(expr.dtype),
        "value": jsonable(hl.eval(expr)),
    }


def inspect_file(path: str) -> Dict[str, Any]:
    """
    Read the first lines of an exported TSV.

    :param path: File path.
    :return: The page's record for the file.
    """
    head = []
    with hl.hadoop_open(path, "r") as f:
        for line in f:
            head.append(line.rstrip("\n")[: STRING_CAP * 8])
            if len(head) >= TSV_PREVIEW_LINES:
                break
    return {"kind": "file", "written": modification_time(path), "head": head}


def collect(out_path: str) -> None:
    """
    Inspect every pipeline input and output and write the JSON.

    :param out_path: Where to write, GCS or local.
    :return: None.
    """
    hl.init(log="/generate_flowchart.log", tmp_dir="gs://gnomad-tmp-4day")
    hl.default_reference(BUILD)
    resources = {}
    for kwargs in [{}, {"skip_pre_rank_metrics": True}]:
        for key, step in describe_step_resources(load_pipeline_steps(**kwargs)).items():
            for r in step["inputs"] + step["outputs"]:
                resources[r["path"]] = r["kind"]

    inspectors = {"ht": inspect_table, "he": inspect_expression, "file": inspect_file}
    tables = {}
    for path, kind in resources.items():
        logger.info("Inspecting %s...", path)
        try:
            tables[path] = inspectors[kind](path)
        except Exception as err:  # pylint: disable=broad-except
            logger.warning("Could not inspect %s: %s", path, err)
            tables[path] = {"kind": kind, "error": str(err).splitlines()[0][:500]}

    payload = {
        "collected": datetime.date.today().isoformat(),
        "hail_version": hl.version(),
        "tables": tables,
    }
    with hl.hadoop_open(out_path, "w") as f:
        f.write(json.dumps(payload))
    logger.info("Wrote %d records to %s.", len(tables), out_path)


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------


def _block_pattern(name: str) -> re.Pattern:
    """
    Return the pattern matching one generated block, markers included.

    :param name: Block name, a key of ``MARKERS``.
    :return: Compiled pattern.
    """
    start, end = MARKERS[name]
    return re.compile(re.escape(start) + ".*?" + re.escape(end), re.DOTALL)


def write_block(html: str, name: str, payload: Any) -> str:
    """
    Replace one generated block of the page.

    :param html: Page source.
    :param name: Block name, a key of ``MARKERS``.
    :param payload: Value to embed as ``var <name> = ...;``.
    :return: Updated page source.
    """
    start, end = MARKERS[name]
    pattern = _block_pattern(name)
    if not pattern.search(html):
        raise SystemExit(f"{HTML_PATH} has no {start} ... {end} block.")
    # Escape "</" so no string in the data can close the script element.
    body = json.dumps(payload, separators=(",", ":"), default=str).replace("</", "<\\/")
    block = f"{start}\nvar {name} = {body};\n{end}"
    return pattern.sub(lambda _: block, html, count=1)


def read_block(html: str, name: str) -> Any:
    """
    Read one generated block back out of the page.

    :param html: Page source.
    :param name: Block name, a key of ``MARKERS``.
    :return: The embedded value.
    """
    match = _block_pattern(name).search(html)
    body = match.group(0).split("=", 1)[1].rsplit(";", 1)[0]
    return json.loads(body.replace("<\\/", "</"))


def strip_volatile(value: Any) -> Any:
    """
    Drop the keys that change without the pipeline changing.

    :param value: CODE block, or any part of it.
    :return: Copy without ``VOLATILE_KEYS``, and without ``UPSTREAM_DOC_KEYS`` on
        gnomad_methods functions.
    """
    if isinstance(value, dict):
        skip = VOLATILE_KEYS
        if value.get("repo") == "gnomad_methods":
            skip = VOLATILE_KEYS | UPSTREAM_DOC_KEYS
        return {k: strip_volatile(v) for k, v in value.items() if k not in skip}
    if isinstance(value, list):
        return [strip_volatile(v) for v in value]
    return value


def describe_drift(old: Dict[str, Any], new: Dict[str, Any]) -> List[str]:
    """
    Summarise how two CODE blocks differ.

    :param old: CODE block in the page.
    :param new: CODE block extracted now.
    :return: One line per differing step, function, or flag.
    """
    lines = []
    for section, key in (("steps", "id"), ("flags", "dest")):
        a = {r[key]: r for r in old[section]}
        b = {r[key]: r for r in new[section]}
        for k in sorted(set(a) | set(b)):
            if a.get(k) != b.get(k):
                lines.append(f"{section[:-1]} {k} changed")
    for k in sorted(set(old["functions"]) | set(new["functions"])):
        if old["functions"].get(k) != new["functions"].get(k):
            lines.append(f"function {k} changed")
    if old.get("globals") != new.get("globals"):
        lines.append("globals changed")
    return lines


def main(args: argparse.Namespace) -> None:
    """Dispatch to the requested mode."""
    if args.collect:
        collect(args.out)
        return

    with open(args.html) as f:
        html = f.read()

    if args.inject_tables:
        with open(args.inject_tables) as f:
            html = write_block(html, "TABLES", json.load(f))
        with open(args.html, "w") as f:
            f.write(html)
        logger.info("Injected table records into %s.", args.html)
        return

    code = extract(allow_dirty=args.allow_dirty or args.check)
    problems = validate_notes(code)

    if args.check:
        drift = describe_drift(
            strip_volatile(read_block(html, "CODE")), strip_volatile(code)
        )
        embedded_notes = read_block(html, "NOTES")
        if embedded_notes != json.loads(json.dumps(resolved_notes(code))):
            drift.append("notes.py differs from the embedded notes")
        for line in problems + drift:
            logger.error(line)
        if problems or drift:
            raise SystemExit("The flowchart is stale. Re-run with --extract.")
        logger.info("The flowchart matches the code.")
        return

    if problems:
        raise SystemExit("Notes do not match the code:\n  " + "\n  ".join(problems))
    html = write_block(html, "CODE", code)
    html = write_block(html, "NOTES", resolved_notes(code))
    with open(args.html, "w") as f:
        f.write(html)
    logger.info(
        "Wrote %d steps, %d flags, and %d functions to %s.",
        len(code["steps"]),
        len(code["flags"]),
        len(code["functions"]),
        args.html,
    )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--extract",
        help="Read the pipeline structure from the code and write it to the page.",
        action="store_true",
    )
    mode.add_argument(
        "--check",
        help="Fail if the page no longer matches the code or the notes.",
        action="store_true",
    )
    mode.add_argument(
        "--collect",
        help="Inspect every pipeline input and output Table. Run on Dataproc.",
        action="store_true",
    )
    mode.add_argument(
        "--inject-tables",
        help="Local JSON written by --collect to splice into the page.",
        type=str,
    )
    parser.add_argument(
        "--out",
        help="Where --collect writes its JSON.",
        type=str,
    )
    parser.add_argument(
        "--allow-dirty",
        help="Let --extract read files with uncommitted changes.",
        action="store_true",
    )
    parser.add_argument(
        "--html",
        help="Page to read and rewrite.",
        type=str,
        default=HTML_PATH,
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    if args.collect and not args.out:
        parser.error("--collect needs --out.")
    main(args)
