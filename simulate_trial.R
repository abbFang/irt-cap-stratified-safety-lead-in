# ============================================================
# Stratified Block Randomization — Clinical Trial Simulation
# ============================================================

# ---- Parameters (edit here) --------------------------------
treatment_arms  <- c("Arm_A", "Arm_B", "Arm_C")
strat_factors   <- list(
  Factor1 = c("F1_L1", "F1_L2"),
  Factor2 = c("F2_L1", "F2_L2")
)
block_size      <- 3     # = n_arms: appropriate for open-label; use n_arms*2+ for blinded
stratum_weights <- c(1, 9, 9, 20)  # relative enrollment per stratum (NULL = equal)
min_per_arm     <- 10    # Objective 1: minimum subjects per arm
n_sim           <- 1000  # simulations used for Objective 1
seed_example    <- 42    # seed for the illustrative randomization list
random_alloc    <- FALSE # TRUE = multinomial stratum draw each sim; FALSE = deterministic largest-remainder
# ------------------------------------------------------------

# ---- Output directory --------------------------------------
out_dir <- "Outputs"   # relative to the working directory; change if needed
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(out_dir, "simulation_log.txt")
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({ sink(type = "output"); sink(type = "message"); close(log_con) }, add = TRUE)
# ------------------------------------------------------------

n_arms         <- length(treatment_arms)
reps_per_block <- block_size / n_arms

stopifnot("block_size must be a multiple of number of arms" =
            block_size %% n_arms == 0)

# ---- Enumerate Strata ----
strata            <- expand.grid(strat_factors, stringsAsFactors = FALSE)
strata$Stratum_ID <- paste0("S", seq_len(nrow(strata)))
n_strata          <- nrow(strata)

if (!is.null(stratum_weights)) {
  stopifnot("stratum_weights length must equal number of strata" =
              length(stratum_weights) == n_strata)
  stopifnot("stratum_weights must all be positive" =
              all(stratum_weights > 0))
}

cat("=== Strata ===\n")
print(strata, row.names = FALSE)

# ---- Stratum Allocation Function ---------------------------
# Distributes N subjects across strata proportional to weights.
# Uses largest-remainder method to ensure sum(alloc) == N exactly.
compute_alloc <- function(N, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, n_strata)
  props   <- weights / sum(weights)
  alloc   <- floor(props * N)
  rem     <- N - sum(alloc)
  if (rem > 0) {
    frac_parts        <- props * N - alloc
    top_idx           <- order(frac_parts, decreasing = TRUE)[seq_len(rem)]
    alloc[top_idx]    <- alloc[top_idx] + 1
  }
  alloc
}

# ---- Core Simulation Function ------------------------------
# Generates a master randomization list for N total subjects.
# Partial final blocks within a stratum are retained as-is.
# When rand_alloc = TRUE, stratum sizes are drawn from a multinomial
# distribution so that enrollment proportions vary across simulation runs.
simulate_rand <- function(N, weights = stratum_weights, seed = NULL,
                          rand_alloc = random_alloc) {
  if (!is.null(seed)) set.seed(seed)

  if (rand_alloc) {
    w <- if (is.null(weights)) rep(1, n_strata) else weights
    alloc <- as.vector(rmultinom(1, N, prob = w / sum(w)))
  } else {
    alloc <- compute_alloc(N, weights)
  }

  rows <- vector("list", n_strata)

  for (s in seq_len(n_strata)) {
    n_s      <- alloc[s]
    if (n_s == 0L) { rows[[s]] <- NULL; next }  # stratum drew no subjects (possible under random_alloc)
    n_blocks <- ceiling(n_s / block_size)

    # Each block is a random permutation of treatment arms
    seq_full <- unlist(lapply(seq_len(n_blocks), function(b) {
      sample(rep(treatment_arms, reps_per_block))
    }))

    # Trim to stratum allocation (allows partial final block)
    seq_s        <- seq_full[seq_len(n_s)]
    block_id     <- ceiling(seq_len(n_s) / block_size)
    pos_in_block <- ((seq_len(n_s) - 1L) %% block_size) + 1L
    is_partial   <- block_id == n_blocks & (n_s %% block_size != 0L)

    rows[[s]] <- data.frame(
      Factor1              = strata$Factor1[s],
      Factor2              = strata$Factor2[s],
      Stratum_ID           = strata$Stratum_ID[s],
      Stratum_N            = n_s,
      Block_ID             = block_id,
      Position_in_Block    = pos_in_block,
      Treatment_Assignment = seq_s,
      Is_Partial_Block     = is_partial,
      stringsAsFactors     = FALSE
    )
  }

  out            <- do.call(rbind, rows)
  out$Subject_ID <- sprintf("%03d", seq_len(nrow(out)))
  rownames(out)  <- NULL
  out[, c("Subject_ID", "Factor1", "Factor2", "Stratum_ID", "Stratum_N",
          "Block_ID", "Position_in_Block", "Treatment_Assignment",
          "Is_Partial_Block")]
}

# ============================================================
# Objective 1: Minimum N for All Arms >= min_per_arm Subjects
# ============================================================

cat("\n\n=== Objective 1: Minimum N for Each Arm >=", min_per_arm, "Subjects ===\n")

N_lo <- n_arms * min_per_arm               # analytic lower bound
N_hi <- N_lo + n_strata * block_size * 10  # upper bound — error if reached

coverage_results <- data.frame(N = integer(0), coverage = numeric(0))

for (N_cand in N_lo:N_hi) {
  n_pass <- sum(vapply(seq_len(n_sim), function(sim) {
    rl <- simulate_rand(N_cand, seed = sim, rand_alloc = random_alloc)
    min(table(rl$Treatment_Assignment)) >= min_per_arm
  }, logical(1)))

  cov <- n_pass / n_sim
  coverage_results <- rbind(coverage_results,
                            data.frame(N = N_cand, coverage = cov))
  if (cov >= 1.0) break
}

# Guard: error if upper bound reached without achieving full coverage
if (!any(coverage_results$coverage >= 1.0)) {
  stop(sprintf(
    paste("Upper bound N_hi = %d was reached without achieving 100%% coverage",
          "(best coverage = %.1f%%).",
          "Increase N_hi in the script or reduce min_per_arm."),
    N_hi, max(coverage_results$coverage) * 100
  ))
}

cat("\nCoverage by N:\n")
print(coverage_results, row.names = FALSE)

min_N <- coverage_results$N[which(coverage_results$coverage >= 1.0)[1]]
cat(sprintf(
  "\n>>> Minimum total N to guarantee all %d arms have >= %d subjects: N = %d <<<\n",
  n_arms, min_per_arm, min_N
))

# Save coverage table
write.csv(coverage_results,
          file.path(out_dir, "obj1_coverage_by_N.csv"),
          row.names = FALSE)
cat("Saved: obj1_coverage_by_N.csv\n")

# ---- Coverage Curve Plot ----
png(file.path(out_dir, "obj1_coverage_curve.png"),
    width = 800, height = 500, res = 110)
plot(
  coverage_results$N, coverage_results$coverage * 100,
  type = "b", pch = 19, col = "steelblue", lwd = 2,
  xlab = "Total Sample Size (N)",
  ylab = "Coverage (% of simulations where all arms pass)",
  main = sprintf(
    "Min N for All Arms >= %d Subjects\n%d strata (unequal enrollment) | %d arms | block size = %d",
    min_per_arm, n_strata, n_arms, block_size
  ),
  ylim = c(0, 105),
  xaxt = "n"
)
axis(1, at = coverage_results$N)
abline(h = 100, lty = 2, col = "firebrick", lwd = 1.5)
abline(v = min_N, lty = 2, col = "darkgreen", lwd = 1.5)
text(min_N + 0.15, 55,
     paste0("Min N = ", min_N),
     col = "darkgreen", pos = 4, cex = 0.9)
grid(nx = NA, ny = NULL, lty = "dotted", col = "gray85")
dev.off()
cat("Saved: obj1_coverage_curve.png\n")

# ============================================================
# Illustrative Randomization List at min_N
# ============================================================

cat(sprintf("\n\n=== Master Randomization List (N = %d, seed = %d) ===\n",
            min_N, seed_example))
rand_list <- simulate_rand(min_N, seed = seed_example, rand_alloc = FALSE)  # deterministic for reproducible illustration
print(rand_list, row.names = FALSE)

write.csv(rand_list,
          file.path(out_dir, "randomization_list.csv"),
          row.names = FALSE)
cat("Saved: randomization_list.csv\n")

# ---- Stratum sizes at min_N ----
cat("\n=== Stratum Allocation at N =", min_N, "===\n")
alloc_summary <- data.frame(
  Stratum_ID = strata$Stratum_ID,
  Factor1    = strata$Factor1,
  Factor2    = strata$Factor2,
  Weight     = if (is.null(stratum_weights)) rep(1, n_strata) else stratum_weights,
  N_allocated = compute_alloc(min_N, stratum_weights)
)
alloc_summary$Pct <- round(alloc_summary$N_allocated / min_N * 100, 1)
print(alloc_summary, row.names = FALSE)

# ---- Overall Balance ----
cat("\n=== Overall Treatment Balance ===\n")
overall_tab <- table(rand_list$Treatment_Assignment)
print(overall_tab)
cat("Proportions:\n")
print(round(prop.table(overall_tab), 3))

# ---- Within-Stratum Balance ----
cat("\n=== Within-Stratum Balance ===\n")
stratum_balance <- as.data.frame.matrix(
  table(rand_list$Stratum_ID, rand_list$Treatment_Assignment)
)
stratum_balance$Stratum_ID  <- rownames(stratum_balance)
stratum_balance$Stratum_N   <- alloc_summary$N_allocated
rownames(stratum_balance)   <- NULL
stratum_balance <- stratum_balance[, c("Stratum_ID", "Stratum_N", treatment_arms)]
print(stratum_balance, row.names = FALSE)

write.csv(stratum_balance,
          file.path(out_dir, "balance_by_stratum.csv"),
          row.names = FALSE)
cat("Saved: balance_by_stratum.csv\n")

# ---- Partial Block Report ----
if (any(rand_list$Is_Partial_Block)) {
  cat("\n=== Partial Block Subjects ===\n")
  partial <- rand_list[rand_list$Is_Partial_Block, ]
  print(partial, row.names = FALSE)
  write.csv(partial,
            file.path(out_dir, "partial_block_subjects.csv"),
            row.names = FALSE)
  cat("Saved: partial_block_subjects.csv\n")
} else {
  cat("\nNo partial blocks at N =", min_N, "\n")
}

cat("\n=== All outputs saved to:", out_dir, "===\n")
