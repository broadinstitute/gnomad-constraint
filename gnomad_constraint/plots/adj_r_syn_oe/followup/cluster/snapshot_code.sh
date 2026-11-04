# Sourced by the launchers after they set $ZIP. Builds the --pyfiles zip from the code
# as it is now, and points $P (pipeline scripts) and $S (cluster scripts) at a copy
# taken at the same moment. Jobs are submitted from those copies, so code edited while
# a chain runs doesn't reach the jobs it hasn't submitted yet: a script from disk and
# utilities from an older zip don't match.
REPO=/Users/jgoodric/PycharmProjects/gnomad-constraint
SNAP=$(mktemp -d /tmp/rr_code.XXXXXX)
cp -R "$REPO/gnomad_constraint/pipeline" "$SNAP/pipeline"
cp -R "$REPO/gnomad_constraint/plots/adj_r_syn_oe/followup/cluster" "$SNAP/cluster"
P=$SNAP/pipeline
S=$SNAP/cluster
rm -f "$ZIP"
(cd "$REPO" && zip -rq "$ZIP" gnomad_constraint/ -x 'gnomad_constraint/experimental/*' 'gnomad_constraint/plots/*' 'gnomad_constraint/flowchart/*' 'gnomad_constraint/renv/*' '*.pyc' '*__pycache__*' '*.RData' '*.DS_Store' '*.png' '*.pdf' '*.ipynb' '*.Rproj.user*')
(cd /Users/jgoodric/PycharmProjects/gnomad_qc && zip -rq "$ZIP" gnomad_qc/ -x '*.pyc' '*__pycache__*' '*.DS_Store' '*.ipynb')
(cd /Users/jgoodric/PycharmProjects/gnomad_methods && zip -rq "$ZIP" gnomad/ -x '*.pyc' '*__pycache__*' '*.DS_Store')
echo "=== code snapshot in $SNAP, zip $ZIP"
