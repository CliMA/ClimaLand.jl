#!/bin/bash
# PreToolUse guard for the unattended weekend loop.
# Backstop for EXCLUDED commands (git/gh/qsub run outside the OS sandbox) and a
# coarse net for destructive rm reaching outside the repo/scratch. The sandbox
# is the real filesystem boundary for every non-excluded command; this hook adds
# workflow/policy guarantees the sandbox does not cover.
# Contract: exit 0 always. Print the deny JSON to stdout to block; print nothing
# to let the command proceed. Runs even under --dangerously-skip-permissions.

input=$(cat)

deny() {
  printf '{"hookSpecificOutput":{"hookEventName":"PreToolUse","permissionDecision":"deny","permissionDecisionReason":"%s"}}\n' "$1"
  exit 0
}

# --- git / gh policy (these run unsandboxed via excludedCommands) ---
grep -Eq 'push[[:space:]]+([^|;&]*[[:space:]])?(--force|-f|--force-with-lease)' <<<"$input" \
  && deny "Force-push is blocked by weekend-loop policy."
grep -Eq 'push[[:space:]][^|;&]*\b(main|master)\b' <<<"$input" \
  && deny "Pushing to main/master is blocked by weekend-loop policy."
grep -Eq 'reset[[:space:]]+--hard' <<<"$input" \
  && deny "git reset --hard is blocked by weekend-loop policy."
grep -Eq 'clean[[:space:]]+-[a-zA-Z]*f' <<<"$input" \
  && deny "git clean -f is blocked by weekend-loop policy."
grep -Eq '\brebase\b' <<<"$input" \
  && deny "git rebase is blocked by weekend-loop policy."
grep -Eq 'gh[[:space:]]+pr[[:space:]]+(merge|close)' <<<"$input" \
  && deny "gh pr merge/close is blocked by weekend-loop policy."

# --- coarse rm safety net for the excluded-command compound hole ---
# Block rm/rmdir/unlink whose target is root, the home dir, or an absolute path
# that is NOT under the repo or the designated scratch tree.
grep -Eq '\b(rm|rmdir|unlink)[[:space:]]+(-[a-zA-Z]*[[:space:]]+)*(/|~|\$HOME)([[:space:]]|/|"|$)' <<<"$input" \
  && ! grep -Eq '(/glade/u/home/arenchon/ClimaLand\.jl|/glade/derecho/scratch/arenchon)' <<<"$input" \
  && deny "rm targeting root/home (outside repo & scratch) is blocked."

grep -Eq '\b(rm|rmdir|unlink)[[:space:]]+(-[a-zA-Z]*[[:space:]]+)*/glade/(campaign|work)' <<<"$input" \
  && deny "rm targeting /glade/campaign or /glade/work is blocked."

exit 0
