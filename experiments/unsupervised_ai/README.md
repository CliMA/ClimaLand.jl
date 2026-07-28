# Unsupervised AI loop — a reproducible workflow

This document explains the **unattended Claude Code loop** we use to make steady,
committed progress on a research/engineering task on an HPC system (NCAR Derecho),
with the AI running autonomously (no per-action approvals) but inside a
**defense-in-depth safety cage**. It is written so another team can copy the
pattern, point it at their own task, and tune the safeguards and cadence.

The concrete example throughout is the ClimaLand PR #1815 loop (climate-responsive
P-model inputs). Anything that is **project-specific is marked `CHANGE THIS`** —
paths, the compute environment, the PBS allocation, and the task prompt.

---

## 1. What it is, in one picture

```
  tmux (survives disconnect)
   └── claude --dangerously-skip-permissions   # no approval prompts...
        │   settings.local.json enables:
        │     • OS sandbox (bubblewrap): filesystem + network allowlist
        │     • permissions.deny: hard bans (force-push, reset --hard, pr merge…)
        │     • PreToolUse hook (guard.sh): policy backstop for git/gh/rm
        │
        └── /loop 30m follow <PROMPT_FILE>       # a recurring cron re-fires the prompt
             every 30 min → one iteration:
               STEP 0 orient (read log, qstat, git)
               STEP 1 one small change  ──edit──►  STEP 2 test on PBS (never on login node)
               STEP 3 commit+push (WIP, [skip ci]) → report on the PR → prune → stop
             state persists in:
               • UNSUPERVISED_LOG.md  (source of truth, in-repo)
               • the GitHub PR        (human-facing summary + rolling comments)
```

The AI is **not** trusted because it is unrestricted; it is unrestricted *because*
it is boxed. The box is four independent layers (Section 4). Heavy compute never
runs where the AI runs — it is submitted to the batch scheduler.

---

## 2. Prerequisites — what to install

| Tool | Why | Notes |
|---|---|---|
| **Claude Code CLI** | the agent | `claude --version` (we ran 2.1.x). Install per Anthropic docs. |
| **tmux** | keep the session alive after you disconnect | any recent version |
| **git** | version control | remote may be **SSH** (`git@github.com:…`) — see gotchas |
| **gh** (GitHub CLI) | PR comments / body updates | file-auth (`~/.config/gh/hosts.yml`), **not** a token in env |
| **bubblewrap** (`bwrap`) | the OS filesystem/network sandbox | verified 0.11.0 on Derecho |
| **socat** | bwrap's network proxy (only if you allow network in-sandbox) | optional; see Section 4.1 |
| compute env | to actually run the science | `CHANGE THIS`: on Derecho it's `module load climacommon` + a pinned Julia; batch = **PBS** (`qsub`) |

Auth, done once, out of band (never by the loop):
- **GitHub**: `gh auth login` to file (scope `repo`), so `gh` and `git` work without
  a token in the environment. The loop **never prints or commits a token**.
- **git push**: if the remote is SSH, your key/agent must be available. (We keep
  `git` *outside* the sandbox so SSH push works — see Section 4.)

---

## 3. The files that make up the workflow

All live in this repo so they are versioned and reproducible:

| File | Role |
|---|---|
| `experiments/integrated/generic_site/unsupervised_loop_prompt.txt` | **The loop task prompt** — the full instructions each iteration follows. `CHANGE THIS` per task. |
| `experiments/integrated/generic_site/UNSUPERVISED_LOG.md` | **Source of truth** — job table, findings, decisions, plan. The AI trusts this over memory. |
| `.claude/unsupervised-loop.settings.json` | Staged **safety config** template (sandbox + deny + hook). Copied to `settings.local.json` to activate. |
| `.claude/hooks/guard.sh` | **PreToolUse backstop** — blocks force-push, reset --hard, rebase, `gh pr merge/close`, dangerous `rm`. Runs even under `--dangerously-skip-permissions`. |
| `.claude/settings.local.json` | The **active** config (git-ignored / local only — do **not** commit to a shared repo). |

> Keep `settings.local.json` local. If committed to a shared repo,
> `sandbox.enabled + failIfUnavailable` would break contributors who lack `bwrap`.

---

## 4. The safety architecture (defense in depth)

Four independent layers. Any one failing does not open the cage.

### 4.1 OS sandbox (bubblewrap) — the real boundary
From `.claude/settings.local.json`:

```jsonc
"sandbox": {
  "enabled": true,
  "failIfUnavailable": true,      // hard-fail if bwrap missing — NO silent unsandboxed fallback
  "allowUnsandboxedCommands": false,
  "filesystem": {
    "allowWrite": [ ".", "~/.julia", "~/.cache",
                    "/glade/derecho/scratch/arenchon", "/tmp" ],   // CHANGE THIS
    "denyWrite": [ ".claude" ]     // the AI cannot edit its own safety config
  },
  "network": {
    "allowedDomains": [ "github.com", "*.github.com", "*.githubusercontent.com",
                        "*.julialang.org", "julialang-s3.julialang.org" ]  // CHANGE THIS
  },
  "credentials": { "files": [ {"path":"~/.ssh","mode":"deny"},
                              {"path":"~/.aws","mode":"deny"} ] },
  "excludedCommands": [ "git","git *","gh","gh *",
                        "qsub","qsub *","qstat","qstat *","qdel","qdel *" ]
}
```

Key decisions:
- **Writable allowlist must include the language depot.** A repo-only sandbox
  *breaks Julia* — it needs `~/.julia` (compile cache + artifacts). Add scratch
  and `/tmp`. `CHANGE THIS` for your language/project. Verify field names against
  your installed Claude Code version.
- **`excludedCommands`** run *outside* the sandbox. We exclude `git`/`gh` (SSH push
  + gh file-auth don't work sandboxed) and the PBS commands (`qsub`/`qstat`/`qdel`).
  Because these bypass the sandbox, they are re-covered by layers 4.2–4.3.
- **`socat`**: bwrap's network proxy needs it *only if you allow in-sandbox
  network*. If you can't install it, either (a) build it into scratch, or (b) drop
  the `network` block and let package downloads happen inside PBS jobs (compute
  nodes have internet) instead of on the login node.
- `denyWrite: [".claude"]` — the AI cannot rewrite its own guardrails.

### 4.2 `permissions.deny` — declarative hard bans
```jsonc
"permissions": { "deny": [
  "Bash(git push --force*)", "Bash(git push -f*)", "Bash(git push --force-with-lease*)",
  "Bash(git reset --hard*)", "Bash(git clean -f*)", "Bash(git rebase*)",
  "Bash(gh pr merge*)", "Bash(gh pr close*)", "Bash(sudo *)"
]}
```

### 4.3 PreToolUse hook (`guard.sh`) — policy backstop
A tiny script that runs before every Bash tool call and **runs even under
`--dangerously-skip-permissions`**. It re-covers the excluded (unsandboxed)
commands: blocks force-push, push to `main`/`master`, `reset --hard`, `clean -f`,
`rebase`, `gh pr merge/close`, and `rm`/`unlink` targeting root/home or paths
outside the repo & scratch. Contract: exit 0 always; print a deny-JSON to block,
print nothing to allow. Wire it up:
```jsonc
"hooks": { "PreToolUse": [ { "matcher": "Bash",
  "hooks": [ { "type": "command", "command": "<abs-path>/.claude/hooks/guard.sh" } ] } ] }
```

### 4.4 Prompt-level HARD RULES + PBS-only compute
The task prompt (Section 5) ends with **HARD RULES** the model is told never to
violate (no force-push, no push to main, never run heavy compute on the login
node, ≤2 batch jobs in flight, charge the right allocation, etc.). The most
important operational rule: **all heavy compute goes through the batch scheduler
(PBS), never in the AI's own shell.** The AI edits + submits jobs; the jobs run
elsewhere. This caps the blast radius of anything the AI does directly.

---

## 5. The loop prompt (the task definition)

`unsupervised_loop_prompt.txt` is the whole task, written for an agent that wakes
with no memory of the last run. Its skeleton (copy it for a new task, `CHANGE THIS`
throughout):

- **Framing**: "You run UNATTENDED. Each wake-up make ONE small, safe, committed
  increment, then stop. Correctness > speed."
- **PRIMARY GOAL** + how success is measured (for us: reproduce MODIS LAI → drive a
  RMSE metric down).
- **STEP 0 — ORIENT** (every iteration): read the log (source of truth), `qstat`
  your jobs and reconcile them, `git status`/`log`.
- **STEP 1 — ONE SMALL STEP**: work the current backlog item only; one logical
  change; keep edits modular and in-style.
- **STEP 2 — TEST VIA THE SCHEDULER**: submit at most **one** new batch job; a full
  implement→test→verify cycle spans 2–3 iterations because the scheduler is async.
  Include the exact job-script template + how to set up the compute env.
- **STEP 3 — COMMIT, PUSH, RECORD, REPORT**: format; **one commit per iteration**
  with `[skip ci]`; push to the working branch (never main); update the PR (body =
  living summary, prune thread comments to a rolling window); update the log.
- **BACKLOG** (ordered), **HARD RULES**, **WHEN-STUCK** (3-strike rule: document a
  clean blocker and move on rather than thrash).

Tuning the task = editing this one file. The loop picks up changes on the next
wake (the prompt is re-read from disk each iteration).

---

## 6. Launching the loop

```bash
# 1. Get on the machine and into a persistent session (tmux survives disconnect).
ssh derecho7                      # CHANGE THIS — tmux is per-login-node; note which one
tmux new -s climaland             # or: tmux attach -t climaland

# 2. Activate the safety config (local only — never commit to a shared repo).
cd /glade/u/home/arenchon/GitHub/ClimaLand.jl          # CHANGE THIS
cp .claude/unsupervised-loop.settings.json .claude/settings.local.json

# 3. Start Claude with no approval prompts (safe because of the 4 layers above).
claude --dangerously-skip-permissions

# 4. Inside Claude, start the recurring loop:
/loop 30m follow experiments/integrated/generic_site/unsupervised_loop_prompt.txt
```

`/loop <interval> <prompt>` schedules a **recurring cron** that re-fires the prompt
at the interval (session-only: it lives as long as the Claude session / tmux does).
`follow <file>` just means "read that file and do what it says." The first
iteration runs immediately; subsequent ones fire on the cron.

To **reattach** after disconnecting: `ssh <same-login-node>` → `tmux attach -t climaland`.

---

## 7. How one iteration works (and where state lives)

Because each wake is memoryless, all cross-iteration state is on disk:

- **`UNSUPERVISED_LOG.md`** (in-repo, committed) — the **source of truth**: a job
  table (`job_id | item | purpose | status | result`), findings, decisions, and
  the plan for the current item. The AI reads this first every iteration and
  trusts it over its own recollection.
- **The GitHub PR** — the human-facing surface:
  - **Opening comment (PR body) = living summary**: backlog checklist, current
    metric table, key findings. Regenerated each iteration.
  - **Thread comments = rolling window** (we keep ~10; the loop deletes its own
    oldest comment after posting a new one). The full history stays in the log, so
    pruning loses nothing.
  - Each comment also reports **iteration wall-time** and the **actual CPU-hours**
    of any job that finished (`qstat -xf <id>` walltime × ncpus).

Async by design: submitting a job in one iteration and reading its result in a
later one is normal. `[skip ci]` keeps WIP pushes from triggering CI; CI is run
deliberately only when the PR is ready (for us, also gated by a `launch buildkite`
label).

---

## 8. Tuning knobs (what to change, and where)

| Want to change… | Do this |
|---|---|
| **The task** | Edit `unsupervised_loop_prompt.txt` (goal, backlog, test recipe, hard rules). Re-read automatically next wake. |
| **Loop frequency** | Re-run `/loop <interval> …` (e.g. `30m`, `1h`), or delete/recreate the cron. Match the interval to your job turnaround. |
| **Work-longer vs wait** | The AI is gated by async job turnaround, not its own runtime; a shorter interval mostly means "check sooner", not "more work per wake". |
| **Writable/network allowlist** | `settings.local.json` → `sandbox.filesystem.allowWrite` / `network.allowedDomains`. Keep it minimal. |
| **Hard bans** | `permissions.deny` + `.claude/hooks/guard.sh`. Add patterns; the hook is a regex backstop. |
| **Commit/comment policy** | One-commit-per-iteration and rolling-comment count are set in the prompt's STEP 3. |
| **Compute target / allocation** | STEP 2 job-script template in the prompt (queue, `ncpus`, account, module env). `CHANGE THIS`. |
| **CI gating** | `[skip ci]` on WIP; drop it (and add your CI trigger) only when ready. |

---

## 9. Adapting to your own project — checklist

1. Write your own `*_loop_prompt.txt`: goal + measurable success + ordered backlog
   + a concrete **STEP 2 test recipe for your compute env** + hard rules.
2. Create a `*_LOG.md` source-of-truth stub (job table + sections) and commit it.
3. Copy `.claude/unsupervised-loop.settings.json` → adjust `allowWrite` (your
   depot/cache), `allowedDomains`, `excludedCommands`, and the hook path.
4. Review `.claude/hooks/guard.sh` — it is regex/best-effort; add your own bans.
5. Confirm out-of-band auth (git push + gh) works **without** secrets in the env.
6. Open the PR once, by hand; the loop only ever **updates** it (never opens a new one).
7. Start in tmux with `--dangerously-skip-permissions`, `/loop`, and watch the
   first few iterations before leaving it unattended.

---

## 10. Operating notes & gotchas (learned the hard way)

- **tmux is per-login-node.** Reattach from the *same* node you started on.
- **`gh` may not be on PATH** in a bare shell — we prefix calls with the module
  load: `bash -lc 'module load gh; gh pr comment …'`. `CHANGE THIS` for your host.
- **Never run heavy compute in the AI's shell.** Login-node-safe = formatting,
  `git`, `qsub`, quick syntax/`using` load checks. Everything that integrates the
  model → a batch job.
- **Don't edit the driver/src while a batch job is running** — the job loads the
  *live* file per site, so a mid-run edit corrupts or breaks later sites. Wait for
  the job to finish (it's an explicit rule in the prompt).
- **History is disposable.** Many `[skip ci]` WIP commits are fine — GitHub
  **"Squash and merge"** collapses them into one clean commit at the end. Do **not**
  rewrite/`rebase`/force-push the shared branch (it's a hard rule and needs a
  force-push).
- **Secrets**: file-based auth only; never echo/commit a token. `~/.ssh`/`~/.aws`
  are denied in the sandbox.
- **Verify the settings schema** against your installed Claude Code version — field
  names (e.g. sandbox filesystem keys) can differ between releases.

---

## 11. Stopping / pausing

- **Pause the loop**: cancel the cron (in Claude, the loop can `CronDelete` it; or
  just stop the Claude session).
- **Stop entirely**: exit Claude / kill the tmux session. In-repo state (the log,
  the branch, the PR) persists; you can resume later by relaunching and `/loop`-ing
  the same prompt.
- **Revert the safety config**: `rm .claude/settings.local.json` (returns to your
  normal interactive Claude settings).

---

*Maintained alongside the ClimaLand unsupervised loop (PR #1815). The live example
of everything above is in `experiments/integrated/generic_site/`.*
