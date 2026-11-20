# irt-cap-stratified-safety-lead-in

A Claude Code skill for stress-testing a small, stratified safety lead-in and finding the **maximum sample size to allow in IRT** so every arm reaches its minimum enrollment.

**Example question:** A safety lead-in needs 10 subjects per arm at 1:1:1 with block size 3. With stratification factors, enrollment splits into small strata with partial blocks, so a perfect 10/10/10 is not guaranteed. What total N should IRT allow so that every arm reaches at least 10?

## How it works

1. You provide arms, stratification factors, block size, expected stratum enrollment mix, and the per-arm minimum.
2. The skill simulates stratified block randomization at increasing total N (1,000 runs per candidate by default).
3. It reports the smallest N where every run has at least the minimum per arm. That N is the recommended IRT enrollment cap. At that N some arms will exceed the minimum.

The result depends on the assumed stratum enrollment mix (`stratum_weights`). Unequal or uncertain weights raise the required N.

## Files

- `SKILL.md`: the skill definition (parameters, steps, Objective 1).
- `simulate_trial.R`: base-R reference implementation. Edit the parameters at the top and run `Rscript simulate_trial.R`. Outputs (CSVs, coverage plot, log) are written to `Outputs/`.

## Install

Copy `SKILL.md` into a folder named `irt-cap-stratified-safety-lead-in` under your Claude Code skills directory (for example `~/.claude/skills/`), keeping `simulate_trial.R` alongside it. Start a new session; the skill triggers on phrases such as "IRT cap", "safety lead-in", "minimum per arm", and "stratified randomization simulation".

## Disclaimer

Illustrative methodology. Validate against your IRT vendor's actual randomization behavior before setting an enrollment cap for a live study.
