#!/bin/sh
# SessionStart hook: make sure the shared gnomAD skills
# (https://github.com/broadinstitute/gnomad-skills) are installed as personal
# skills and keep them current, so /gnomad-pr-review, gnomad-commit,
# gnomad-pr-description and gnomad-clean-git-history are available here.
#
# The skills live in a checkout of gnomad-skills: the one ~/.claude/skills is
# already linked to, if any, or else $GNOMAD_SKILLS_DIR (default
# ~/.claude/gnomad-skills), cloned on first use. Every session runs the
# checkout's scripts/link-skills.sh, which links each skill into ~/.claude/skills
# and leaves alone any skill directory it did not make. The checkout is pulled
# at most once a day. The hook is silent when offline or already current, and
# prints only when there is something to tell the user. Set
# GNOMAD_SKIP_SKILLS_HOOK=1 to turn it off.

[ -n "$GNOMAD_SKIP_SKILLS_HOOK" ] && exit 0
[ -n "$HOME" ] || exit 0
mkdir -p "$HOME/.claude" 2>/dev/null

repo="https://github.com/broadinstitute/gnomad-skills.git"
skills="$HOME/.claude/skills"
link="$skills/gnomad-pr-review"
stamp="$HOME/.claude/.gnomad-skills-last-check"
reload_skills='"hookSpecificOutput": {"hookEventName": "SessionStart", "reloadSkills": true}'
export GIT_TERMINAL_PROMPT=0

# Plain stdout from a SessionStart hook only reaches Claude, so messages for
# the user go in the JSON systemMessage field.
json_escape() {
    printf '%s' "$*" | tr '\n' ' ' | sed 's/\\/\\\\/g; s/"/\\"/g'
}
notify() {
    printf '{"systemMessage": "%s"%s}\n' "$(json_escape "$1")" "${2:+, $2}"
}
is_checkout() {
    [ -d "$1/.git" ] && [ -f "$1/scripts/link-skills.sh" ]
}
# The skill links and where they point (link-skills.sh rewrites every link, so
# their timestamps change on each run).
links() {
    for f in "$skills"/*; do
        [ -L "$f" ] && printf '%s %s\n' "$f" "$(readlink "$f")"
    done
}

# link-skills.sh links ~/.claude/skills/<name> to <checkout>/skills/<name>, so
# a link into a checkout points back to it. Any other link is ignored.
src="${GNOMAD_SKILLS_DIR:-$HOME/.claude/gnomad-skills}"
if [ -L "$link" ]; then
    linked=$(dirname "$(dirname "$(readlink "$link")")")
    is_checkout "$linked" && src="$linked"
fi
today=$(date +%Y%m%d)
checked=""
[ "$(cat "$stamp" 2>/dev/null)" = "$today" ] || checked=1

if ! is_checkout "$src"; then
    if [ -e "$src" ]; then
        # Never clone onto something that is already there; say so once a day.
        [ -n "$checked" ] || exit 0
        echo "$today" > "$stamp"
        notify "The shared gnomAD skills are not installed: $src exists but is not a gnomad-skills checkout. Move it aside, or set GNOMAD_SKILLS_DIR to another location."
        exit 0
    fi
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

# Pull at most once a day, and count the day as checked even if the pull fails
# (offline, or a checkout that can't fast-forward), so it isn't retried every
# session. A pull that brings in commits can change a skill's text without
# changing any link, so it asks for a reload too.
changed=""
if [ -n "$checked" ]; then
    echo "$today" > "$stamp"
    head=$(git -C "$src" rev-parse HEAD 2>/dev/null)
    git -C "$src" pull -q --ff-only >/dev/null 2>&1
    [ "$(git -C "$src" rev-parse HEAD 2>/dev/null)" = "$head" ] || changed=1
fi

# Link on every session start: it is local and fast, picks up added or renamed
# skills, and links a skill as soon as whatever blocked it is moved aside.
before=$(links)
skipped=$(sh "$src/scripts/link-skills.sh" 2>&1 >/dev/null)
[ "$(links)" = "$before" ] || changed=1

if [ -n "$skipped" ] && [ -n "$checked" ]; then
    notify "gnomad-skills: $skipped" "${changed:+$reload_skills}"
elif [ -n "$changed" ]; then
    printf '{%s}\n' "$reload_skills"
fi

exit 0
