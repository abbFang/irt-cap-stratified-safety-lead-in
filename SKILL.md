---
name: irt-cap-stratified-safety-lead-in
description: Use when a study team needs the maximum sample size to allow in IRT so every arm reaches a minimum enrollment (e.g., 10 per arm at 1:1:1) in a small stratified safety lead-in, run-in, or dose-escalation cohort. Simulates stratified block randomization with partial blocks and unequal stratum enrollment, and finds the smallest N where every simulated run meets the per-arm minimum. Triggers on "IRT cap", "enrollment cap", "safety lead-in", "minimum per arm", "block size", "stratified randomization simulation", "partial blocks", "randomization list".
---

# Skill: Right-Sizing the IRT Cap for Safety Lead-ins

## Purpose

Stress-test a small, stratified safety lead-in (e.g., 10 subjects per arm, 1:1:1) and find the maximum sample size to allow in the IRT system so every arm is guaranteed the minimum enrollment. Because stratification factors create small strata and partial blocks, a perfect split (e.g., 10/10/10) is not guaranteed even with `block_size = n_arms`.

Under the hood, the skill simulates participant randomization using stratified block randomization. It produces a randomization list that is balanced within each stratum and across treatment arms, mimicking real-world clinical trial practice (e.g., as used in IRT/IVRS systems), then searches for the smallest total N that meets the per-arm minimum in every simulation run (Objective 1). That N is the recommended IRT enrollment cap.

## Parameters

| Parameter | Description | Example |
|---|---|---|
| `sample_size` | Total number of participants to randomize | `120` |
| `treatment_arms` | Number and names of treatment arms (including control) | `c("Placebo", "Low Dose", "High Dose")` |
| `stratification_factors` | Named factors and their levels | `list(sex = c("M", "F"), age_group = c("<65", ">=65"))` |
| `block_size` | Block size within each stratum; must be a positive multiple of number of arms | `6` |
| `stratum_weights` | Relative enrollment proportions per stratum, in the same order as `expand.grid()` produces strata. Need not sum to 1 — they are normalized internally. Use `NULL` for equal allocation. | `c(1, 9, 9, 20)` |
| `min_per_arm` | *(Objective 1)* Minimum required subjects per treatment arm | `20` |
| `random_alloc` | If `TRUE`, draw stratum sizes stochastically from a multinomial distribution each simulation run (mimics real enrollment variability). If `FALSE` (default), use deterministic largest-remainder proportional allocation. | `FALSE` |
| `seed` | Random seed for reproducibility | `42` |

## Real-World Design Notes

- Stratified block randomization is the standard in Phase II/III trials. Each stratum gets its own independent randomization sequence.
- Participants are assigned sequentially within their stratum as they enroll — they are not pre-allocated across strata.
- **Stratum enrollment is inherently unequal.** In real trials, strata fill at different rates determined by patient population prevalence, not by design. A rare stratum (e.g., young females with a specific biomarker) may accumulate only 1–2 subjects while a common stratum reaches 20+. Use `stratum_weights` to reflect this; equal allocation (`NULL`) is a simplifying assumption suitable only for preliminary exploration.
- **Block size selection depends on blinding.** In blinded (double-blind) trials, block size is kept confidential and should be at least `2 × n_arms` to prevent deduction of final-position assignments. In **open-label trials**, allocation concealment is not required; `block_size = n_arms` (one copy of each arm per block) is acceptable and preferred when sample sizes are small, as it guarantees the tightest within-stratum balance after each complete block.
- Partial blocks at the end of a stratum are expected and allowed — enrollment rarely fills every stratum perfectly.
- Block sizes are typically kept confidential from site staff in blinded studies; only the sponsor/IRT system holds the full master list.

---

## Instructions

When the user invokes this skill, follow these steps:

### 1. Confirm Parameters
Ask the user to confirm or provide all parameters above. Validate:
- `block_size` is a positive multiple of the number of treatment arms.
- Stratification factor levels are mutually exclusive and exhaustive.
- If `stratum_weights` is provided, its length equals the number of strata (`prod(lengths(stratification_factors))`).
- If `min_per_arm` is provided, note that Objective 1 will be run.

### 2. Enumerate Strata
Generate all combinations of stratification factor levels using `expand.grid()`. Each unique combination is one stratum. Label strata sequentially (`S1`, `S2`, …).

### 3. Compute Stratum Allocation
For a given total `N`, distribute subjects across strata using `stratum_weights`. Two modes are supported:

**Deterministic mode (`random_alloc = FALSE`, default):**
- Normalize weights to proportions: `props <- stratum_weights / sum(stratum_weights)`
- Initial allocation: `alloc <- floor(props * N)`
- Distribute the remainder (1 subject at a time) to the strata with the largest fractional remainders until `sum(alloc) == N`
- Produces the same stratum sizes every simulation run given the same `N`
- When `stratum_weights = NULL`, treat as equal weights

**Random mode (`random_alloc = TRUE`):**
- Draw stratum sizes from a multinomial distribution: `alloc <- as.vector(rmultinom(1, N, prob = weights / sum(weights)))`
- Stratum sizes vary stochastically across simulation runs, mimicking real-world enrollment where patients arrive in unpredictable proportions
- The expected allocation matches `stratum_weights`, but any given run may deviate — small strata may receive zero or many subjects by chance
- When `stratum_weights = NULL`, treat as equal weights (uniform multinomial)

### 4. Generate Randomization Sequences Per Stratum
For each stratum independently:
- Determine the number of complete blocks needed: `ceiling(alloc[s] / block_size)`
- Generate each block by randomly permuting the treatment arm labels so each arm appears `block_size / n_arms` times — use `sample()` within each block
- Concatenate all blocks; trim to `alloc[s]` subjects (partial final block retained as-is)
- Tag each row with `Stratum_ID`, `Block_ID`, and within-stratum position

### 5. Assemble the Master Randomization List
Combine all strata into a single data frame, assign global `Subject_ID` sequentially, and include:

| Column | Description |
|---|---|
| `Subject_ID` | Global sequential ID (001, 002, …) |
| *(one col per factor)* | Stratum-defining factor values |
| `Stratum_ID` | Stratum label (S1, S2, …) |
| `Block_ID` | Block number within stratum |
| `Position_in_Block` | Position within block (1 to block_size) |
| `Treatment_Assignment` | Assigned treatment arm label |
| `Is_Partial_Block` | `TRUE` if this subject falls in a partial final block |

### 6. Summarize Balance
Print two summaries:
1. **Overall balance** — count and percentage per treatment arm across all subjects.
2. **Within-stratum balance** — count per treatment arm within each stratum, including stratum size and highlighting any arms with counts below the overall average (due to partial blocks or small strata).

### 7. Write the R Script

Write a self-contained R script (base R only, no external packages unless the user requests `tidyverse`) to `<project_root>/run_trial_sim.R`. The script must:
- Accept all parameters as variables at the top (easy to modify).
- Produce the master randomization list as a `data.frame`.
- Print the randomization list and balance summaries using `table()` and `prop.table()`.
- Write all outputs (CSV files, PNG coverage curve, log) to `<project_root>/Outputs/`.
- Be runnable end-to-end via `Rscript`.

### 8. Locate and Execute R

After writing the script, find `Rscript.exe` and run it. Search in this order:

1. **Registry** (most reliable on Windows):
   ```
   reg query "HKLM\SOFTWARE\R-core\R" /v InstallPath
   reg query "HKLM\SOFTWARE\WOW6432Node\R-core\R" /v InstallPath
   reg query "HKCU\SOFTWARE\R-core\R" /v InstallPath
   ```
   Append `\bin\Rscript.exe` to the returned `InstallPath`.

2. **Common install locations** (fallback):
   - `C:\Program Files\R\` — list subdirectories, pick the highest version, append `\bin\Rscript.exe`
   - `C:\Users\<username>\AppData\Local\R\`

3. **PATH**:
   ```
   where Rscript
   ```

Once found, execute:
```
"<path_to_Rscript.exe>" "<path_to_run_trial_sim.R>"
```

If R cannot be found, tell the user and ask them to run `! "<path_to_Rscript.exe>" "<script_path>"` in the prompt, or open the script in RStudio and press **Ctrl+Shift+S** to source it.

### 9. Report Results

After execution, read `Outputs/simulation_log.txt` and report back:
- The minimum total N found
- The coverage table (N vs. coverage %)
- Paths to the saved output files

---

## Objective 1: Minimum Sample Size for Guaranteed Per-Arm Coverage

**Goal:** Find the smallest total `sample_size` such that, across all realizations of the randomization, every treatment arm is guaranteed to receive at least `min_per_arm` subjects.

**Approach (simulation-based):**

1. For a candidate `N`, run the stratified block randomization `n_sim` times (default: 1000) with different seeds, using the specified `stratum_weights` for allocation.
2. In each simulation, record the minimum count across all treatment arms (`min_arm_count`).
3. Compute the proportion of simulations where `min_arm_count >= min_per_arm`. This is the empirical coverage probability.
4. Increment `N` by 1 and repeat until coverage reaches the user-specified threshold (default: 100%).
5. If the upper bound `N_hi` is reached before coverage = 100%, raise an error — do not silently return `NA`.
6. Report the minimum `N` found and a coverage curve (N vs. proportion of simulations meeting the criterion).

**Note:** With unequal `stratum_weights`, very small strata introduce more partial-block imbalance. The minimum N will generally be higher than under equal allocation, and the coverage curve will rise more gradually.

**Note on `random_alloc = TRUE`:** Stochastic stratum allocation adds a second source of variance — beyond block-level randomization, each simulation also draws different stratum sizes from the multinomial. This typically increases the minimum N required to reach 100% coverage and causes the coverage curve to rise more slowly. The tail risk is that a rare stratum draws zero subjects in some runs, making per-arm guarantees impossible to satisfy without a sufficiently large N to keep such draws improbable.

**R Output:**
- A coverage curve plot (`plot()` or `ggplot2` if available).
- A table of `N` vs. coverage probability.
- The recommended minimum `N` printed clearly.

---

## Adding Future Objectives

Additional objectives will be appended as new `## Objective N` sections below. Each objective should specify its goal, analytical approach, and expected R output.
