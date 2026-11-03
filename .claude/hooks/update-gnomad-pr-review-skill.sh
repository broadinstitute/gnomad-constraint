#!/bin/sh
# SessionStart hook: install Mike Wilson's gnomad-pr-review skill from gnomad_qc
# as a personal skill (~/.claude/skills/gnomad-pr-review) and keep it current.
# The skill reviews PRs in gnomad_qc, gnomad-constraint, and gnomad_methods, so
# one copy serves all three and this repo ships none.
#
# Only the skill's directory is fetched, into a shallow sparse clone kept under
# ~/.claude, and the skill is a copy of it marked as installed by this hook. A
# gnomad-pr-review skill without that marker is the user's own and is left
# alone; edits to a marked copy are replaced on the next update. The hook
# checks upstream at most once a day, is silent when offline or already
# current, and prints only when there is something to tell the user. Set
# GNOMAD_SKIP_REVIEW_HOOK=1 to turn it off.

[ -n "$GNOMAD_SKIP_REVIEW_HOOK" ] && exit 0
[ -n "$HOME" ] || exit 0

repo="https://github.com/broadinstitute/gnomad_qc.git"
# The skill is in review in gnomad_qc#798. Until it merges, install it from
# that branch; once the branch is gone, the hook falls back to main.
ref="mw/qc_review_skill"
path=".claude/skills/gnomad-pr-review"
dir="$HOME/.claude/skills/gnomad-pr-review"
marker=".installed-by-gnomad-hook"
src="$HOME/.claude/.gnomad-pr-review-skill-src"
stamp="$HOME/.claude/.gnomad-pr-review-skill-last-check"
lock="$HOME/.claude/.gnomad-pr-review-skill.lock"
reload_skills='"hookSpecificOutput": {"hookEventName": "SessionStart", "reloadSkills": true}'

# Let one hook act at a time; see update-review-skill.sh. A lock older than
# five minutes was left by a hook that was killed, so it is cleared.
mkdir -p "$HOME/.claude" 2>/dev/null
now=$(date +%s)
if ! mkdir "$lock" 2>/dev/null; then
    since=$(cat "$lock/since" 2>/dev/null)
    if [ -z "$since" ] || [ $((now - since)) -le 300 ]; then
        exit 0
    fi
    rm -rf "$lock"
    mkdir "$lock" 2>/dev/null || exit 0
fi
echo "$now" > "$lock/since"
trap 'rm -rf "$lock"' EXIT
trap 'exit 1' HUP INT TERM

# Plain stdout from a SessionStart hook only reaches Claude, so messages for
# the user go in the JSON systemMessage field.
json_escape() {
    printf '%s' "$*" | sed 's/\\/\\\\/g; s/"/\\"/g'
}
notify() {
    printf '{"systemMessage": "%s"}\n' "$(json_escape "$@")"
}

# Check upstream at most once a day.
today=$(date +%Y%m%d)
[ "$(cat "$stamp" 2>/dev/null)" = "$today" ] && exit 0

if [ -e "$dir" ] && [ ! -f "$dir/$marker" ]; then
    echo "$today" > "$stamp"
    notify "$dir was not installed by this repo's hook, so the hook does not" \
        "update the gnomad-pr-review skill there. To get automatic updates," \
        "move that directory aside; the next session installs a fresh copy."
    exit 0
fi

# Fetch the skill's directory at the newest commit of the branch, or of main
# once the branch is gone.
if [ ! -d "$src/.git" ]; then
    rm -rf "$src"
    if ! { git init -q "$src" \
        && git -C "$src" remote add origin "$repo" \
        && git -C "$src" sparse-checkout set "$path"; } >/dev/null 2>&1; then
        rm -rf "$src"
        exit 0
    fi
fi
fetched=""
for r in "$ref" main; do
    if git -C "$src" fetch -q --depth 1 --filter=blob:none origin "$r" \
        >/dev/null 2>&1; then
        fetched="$r"
        break
    fi
done
# Offline: say nothing, and try again next session.
[ -n "$fetched" ] || exit 0
git -C "$src" checkout -q --detach FETCH_HEAD >/dev/null 2>&1 || exit 0
echo "$today" > "$stamp"
[ -f "$src/$path/SKILL.md" ] || exit 0

sha=$(git -C "$src" rev-parse HEAD)
[ "$(cut -d' ' -f3 "$dir/$marker" 2>/dev/null)" = "$sha" ] && exit 0

# Copy into place through a temporary directory, so an interrupted copy can't
# leave a partial skill behind.
first_install=""
[ -e "$dir" ] || first_install=1
rm -rf "$dir.tmp"
if mkdir -p "$(dirname "$dir")" \
    && cp -R "$src/$path" "$dir.tmp" \
    && echo "$repo $fetched $sha" > "$dir.tmp/$marker" \
    && rm -rf "$dir" \
    && mv "$dir.tmp" "$dir"; then
    if [ -n "$first_install" ]; then
        printf '{"systemMessage": "%s", %s}\n' "$(json_escape \
            "Installed the gnomad-pr-review skill from gnomad_qc ($fetched) into" \
            "$dir as a personal skill, so it is available in every project.")" \
            "$reload_skills"
    else
        printf '{%s}\n' "$reload_skills"
    fi
else
    rm -rf "$dir.tmp"
fi

exit 0
