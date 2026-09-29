# ============================================================
# BOIN-Hill simulation configuration
#
# Users can modify the simulation/design settings in this file.
#
# Input files:
#   ../data/시나리오_UBOIN.csv   : true toxicity/efficacy scenarios
#   ../data/파라미터_UBOIN.csv   : prior specifications for BOIN-Hill
#
# Main analysis:
#   UBLRM_Paper.R
# ============================================================


# ============================================================
# 1. Reproducibility
# ============================================================

# Random seed
seed <- 42
set.seed(seed)


# ============================================================
# 2. Number of simulation replicates
# ============================================================

# For a quick reproducibility check in Code Ocean:
sim <- 10

# For reproducing the manuscript simulation:
# sim <- 2000


# ============================================================
# 3. Trial design parameters
# ============================================================

para <- list(
  seed = seed,

  # Planned number of patients in Stage 1
  S1 = 12,

  # Planned number of patients in Stage 2
  S2 = 24,

  # Maximum total sample size
  N = 54,

  # Cohort size
  cohort = 3,


  # ----------------------------------------------------------
  # Follow-up settings
  # ----------------------------------------------------------

  # DLT assessment window
  window = 1,


  # ----------------------------------------------------------
  # Stage 1 / Backfill stopping settings
  # ----------------------------------------------------------

  # Maximum number of patients allowed at a dose
  n_cap = 12,

  # Stop escalation when the current dose reaches this number
  # under the stay decision
  n_stop = 9,

  # Maximum total number of patients used during escalation
  N_esc = 30,


  # ----------------------------------------------------------
  # Toxicity / efficacy thresholds
  # ----------------------------------------------------------

  # Target DLT probability
  DLT = 0.25,

  # Maximum acceptable toxicity probability
  T_admissible = 0.30,

  # Minimum acceptable efficacy probability
  E_admissible = 0.20,


  # ----------------------------------------------------------
  # Posterior decision thresholds
  # ----------------------------------------------------------

  # Toxicity cutoff
  CT = 0.95,

  # Efficacy cutoff
  CE = 0.90
)


# ============================================================
# 4. Utility score
# ============================================================

# Joint toxicity/efficacy utilities
#
# Order:
#   (T=0,E=0),
#   (T=0,E=1),
#   (T=1,E=0),
#   (T=1,E=1)
#
score <- c(
  30,
  100,
  0,
  30
)


# ============================================================
# 5. BOIN-Hill prior settings
# ============================================================

# BOIN-Hill is evaluated under both prior specifications.
#
# prior 1 = strong informative prior
# prior 2 = weak informative prior
#
boin_hill_prior_ids <- c(1, 2)

# Other methods do not use the BOIN-Hill prior.
# We use prior 1 only as a placeholder simulation object so that
# those methods are evaluated only once.
default_prior_id <- 1


# ============================================================
# 6. Input files
# ============================================================

# Prior specifications for BOIN-Hill
prior_file <- file.path(
  "..",
  "data",
  "파라미터_UBOIN.csv"
)

# True toxicity and efficacy scenarios
scenario_file <- file.path(
  "..",
  "data",
  "시나리오_UBOIN.csv"
)


# ============================================================
# 7. Read prior specifications
# ============================================================

prior_table <- readr::read_csv(
  prior_file,
  show_col_types = FALSE
)


# ============================================================
# 8. Read simulation scenarios
# ============================================================

scenario_table <- readr::read_csv(
  scenario_file,
  show_col_types = FALSE
)


# ============================================================
# 9. Construct scenarios
# ============================================================

scenarios <- list()

for (i in seq(2, nrow(scenario_table), by = 2)) {

  scenario_name <- as.character(
    scenario_table[[1]][i]
  )

  # True toxicity probability at each dose
  pi_T <- as.numeric(
    unlist(
      scenario_table[
        i,
        3:ncol(scenario_table)
      ],
      use.names = FALSE
    )
  )

  # True efficacy probability at each dose
  pi_E <- as.numeric(
    unlist(
      scenario_table[
        i + 1,
        3:ncol(scenario_table)
      ],
      use.names = FALSE
    )
  )

  scenarios[[scenario_name]] <- list(

    # True toxicity probabilities
    pi_T = pi_T,

    # True efficacy probabilities
    pi_E = pi_E,

    # Dependence parameter used for joint T/E generation
    coff = 1,

    # Utility score
    score = score
  )
}


# ============================================================
# 10. Dose information
# ============================================================

# Number of dose levels
para$J <- ncol(scenario_table) - 2

# Actual dose values
para$doses <- as.numeric(
  unlist(
    scenario_table[
      1,
      3:ncol(scenario_table)
    ],
    use.names = FALSE
  )
)

# Reference dose
para$ref_dose <- as.numeric(
  scenario_table[[ncol(scenario_table)]][1]
)