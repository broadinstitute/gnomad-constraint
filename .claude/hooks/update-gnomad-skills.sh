#!/bin/sh
# SessionStart hook: make sure the shared gnomAD skills
# (https://github.com/broadinstitute/gnomad-skills) are installed as personal
# skills and keep them current, so /gnomad-pr-review, gnomad-commit,
# gnomad-pr-description and gnomad-clean-git-history are available here.
#
# If the skills are already linked into ~/.claude/skills from a checkout of
# gnomad-skills, that checkout is updated. Otherwise the hook clones one into
# $GNOMAD_SKILLS_DIR (default ~/.claude/gnomad-skills). Either way it runs the
# checkout's scripts/link-skills.sh, which links each skill and leaves alone any
# skill directory it did not make. The hook checks upstream at most once a day,
# is silent when offline or already current, and prints only when there is
# something to tell the user. Set GNOMAD_SKIP_SKILLS_HOOK=1 to turn it off.

[ -n "$GNOMAD_SKIP_SKILLS_HOOK" ] && exit 0
[ -n "$HOME" ] || exit 0

repo="https://github.com/broadinstitute/gnomad-skills.git"
link="$HOME/.claude/skills/gnomad-pr-review"
stamp="$HOME/.claude/.gnomad-skills-last-check"
reload_skills='"hookSpecificOutput": {"hookEventName": "SessionStart", "reloadSkills": true}'

# Plain stdout from a SessionStart hook only reaches Claude, so messages for
# the user go in the JSON systemMessage field.
json_escape() {
    printf '%s' "$*" | tr '\n' ' ' | sed 's/\\/\\\\/g; s/"/\\"/g'
}
notify() {
    printf '{"systemMessage": "%s"%s}\n' "$(json_escape "$1")" "${2:+, $2}"
}

# link-skills.sh links ~/.claude/skills/<name> to <checkout>/skills/<name>, so
# an existing link points back to the checkout.
if [ -L "$link" ]; then
    src=$(dirname "$(dirname "$(readlink "$link")")")
else
    src="${GNOMAD_SKILLS_DIR:-$HOME/.claude/gnomad-skills}"
fi
today=$(date +%Y%m%d)

if [ ! -d "$src/.git" ]; then
    # Clone into a temporary directory and move it into place, so a clone cut
    # off by the hook's timeout can't leave a partial checkout behind. Offline:
    # say nothing, and try again next session.
    rm -rf "$src.partial"
    mkdir -p "$(dirname "$src")" 2>/dev/null
    git clone -q "$repo" "$src.partial" >/dev/null 2>&1 \
        && mv "$src.partial" "$src" || { rm -rf "$src.partial"; exit 0; }
    echo "$today" > "$stamp"
    skipped=$(sh "$src/scripts/link-skills.sh" 2>&1 >/dev/null)
    notify "Installed the shared gnomAD skills (/gnomad-pr-review, gnomad-commit, gnomad-pr-description, gnomad-clean-git-history) from $repo into $src and linked them into ~/.claude/skills. ${skipped}" "$reload_skills"
    exit 0
fi

# Check upstream at most once a day.
[ "$(cat "$stamp" 2>/dev/null)" = "$today" ] && exit 0
before=$(git -C "$src" rev-parse HEAD 2>/dev/null)
git -C "$src" pull -q --ff-only >/dev/null 2>&1 || exit 0
echo "$today" > "$stamp"

# Relink every day: it picks up added or renamed skills, and it reports a skill
# directory that blocks a link, such as an old hand-installed copy.
skipped=$(sh "$src/scripts/link-skills.sh" 2>&1 >/dev/null)
if [ -n "$skipped" ]; then
    notify "gnomad-skills: $skipped" "$reload_skills"
elif [ "$(git -C "$src" rev-parse HEAD 2>/dev/null)" != "$before" ]; then
    printf '{%s}\n' "$reload_skills"
fi

exit 0
