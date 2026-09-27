# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "MASS", "shiny", "rhandsontable", "LaplacesDemon",
  "ggplot2", "readxl", "rjags", "future",
  "doParallel", "foreach", "parallel", "BayesLogit",
  "VGAM", "BOIN", "UBCRM", "reshape2",
  "dplyr", "readr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages(
  invisible(lapply(required_packages, library, character.only = TRUE))
)

set.seed(42)


# ============================================================
# 1. Source functions
# ============================================================

source_files <- c(
  "Performance.R",
  "Scenario.R",
  "Make_data.R",
  "load_dependencies.R",
  "Plot.R",
  "UBLRM_Prior.R",
  "UBLRM_Model.R",
  "Backfill_Stage1_Run.R",
  "Stage1_Run.R",
  "Stage2_Run.R"
)

for (file in source_files) {
  source(file.path("Function_ver2", file))
}


# ============================================================
# 2. Trial parameters
# ============================================================

para <- list(
  seed = 123,

  J = 5,
  S1 = 12,
  S2 = 24,
  N = 54,
  cohort = 3,

  window = 1,

  n_cap = 12,
  n_stop = 9,
  N_esc = 30,

  DLT = 0.25,
  T_admissible = 0.30,
  E_admissible = 0.20,

  CT = 0.95,
  CE = 0.90
)


# ============================================================
# 3. Prior
# ============================================================

prior_table <- readr::read_csv(
  "파라미터_UBOIN.csv",
  show_col_types = FALSE
)


get_priors <- function(
  prior_table,
  scenario_name,
  candidate_num = 1
) {

  mean_col <- paste0(
    "Cand",
    candidate_num,
    "_mean(min,a)"
  )

  sd_col <- paste0(
    "Cand",
    candidate_num,
    "_sd(max,b)"
  )

  sub_df <- dplyr::filter(
    prior_table,
    Scenario == scenario_name
  )

  E_prior <- list()
  T_prior <- list()

  for (i in seq_len(nrow(sub_df))) {

    param <- sub_df$Param[i]
    dist <- sub_df$Dist[i]

    mean_val <- sub_df[[mean_col]][i]
    sd_val <- sub_df[[sd_col]][i]

    # Efficacy priors
    if (param %in% c(
      "delta0", "delta1",
      "b0", "b1"
    )) {

      if (dist == "Beta") {

        E_prior[[paste0(param, "_a")]] <- mean_val
        E_prior[[paste0(param, "_b")]] <- sd_val

      } else if (dist == "Normal") {

        E_prior[[paste0(param, "_mean")]] <- mean_val
        E_prior[[paste0(param, "_sd")]] <- sd_val
      }
    }

    # Toxicity priors
    if (param %in% c("a0", "a1")) {

      if (dist == "Uniform") {

        T_prior[[paste0(param, "_min")]] <- mean_val
        T_prior[[paste0(param, "_max")]] <- sd_val

      } else if (dist == "Normal") {

        T_prior[[paste0(param, "_mean")]] <- mean_val
        T_prior[[paste0(param, "_sd")]] <- sd_val
      }
    }
  }

  list(
    E_prior = E_prior,
    T_prior = T_prior
  )
}


# ============================================================
# 4. Scenario setup
# ============================================================

scenario_table <- readr::read_csv(
  "시나리오_UBOIN.csv",
  show_col_types = FALSE
)

score <- c(30, 100, 0, 30)

scenarios <- list()


for (i in seq(2, nrow(scenario_table), by = 2)) {

  scenario_name <- as.character(
    scenario_table[[1]][i]
  )

  pi_T <- as.numeric(
    unlist(
      scenario_table[
        i,
        3:ncol(scenario_table)
      ],
      use.names = FALSE
    )
  )

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
    pi_T = pi_T,
    pi_E = pi_E,
    coff = 1,
    score = score
  )
}


# Dose information
para$J <- ncol(scenario_table) - 2

para$doses <- as.numeric(
  unlist(
    scenario_table[
      1,
      3:ncol(scenario_table)
    ],
    use.names = FALSE
  )
)

para$ref_dose <- as.numeric(
  scenario_table[
    [ncol(scenario_table)]
  ][1]
)


# ============================================================
# 5. Plot scenarios
# ============================================================

par(mfrow = c(2, 4))

for (i in seq_along(scenarios)) {

  scenario_name <- names(scenarios)[i]
  scenario <- scenarios[[scenario_name]]

  prob <- calculate_joint_probabilities(
    scenario
  )

  score_prob <- t(
    apply(
      prob,
      1,
      function(row) row * scenario$score
    )
  )

  utility <- rowSums(score_prob)

  plot_scenario(
    para,
    scenario,
    name = scenario_name,
    utility = utility,
    show_legend = (i == 1)
  )
}


# ============================================================
# 6. Simulation settings
# ============================================================

sim <- 10
prior_ids <- c(1, 2)


# ============================================================
# 7. Generate simulation data
# ============================================================

generated_data <- setNames(
  vector(
    "list",
    length(scenarios)
  ),
  names(scenarios)
)

for (scenario_name in names(scenarios)) {

  generated_data[[scenario_name]] <-
    generate_toxicity_efficacy_time(
      scenarios[[scenario_name]],
      sim * 100
    )
}


# ============================================================
# 8. Build simulation objects
# ============================================================

sim_objects <- setNames(
  vector(
    "list",
    length(prior_ids)
  ),
  paste0("prior", prior_ids)
)


for (p in prior_ids) {

  prior_key <- paste0("prior", p)

  sim_objects[[prior_key]] <- setNames(
    vector(
      "list",
      length(scenarios)
    ),
    names(scenarios)
  )

  for (scenario_name in names(scenarios)) {

    scenario <- scenarios[[scenario_name]]

    priors <- get_priors(
      prior_table,
      scenario_name,
      candidate_num = p
    )

    sim_objects[
      [prior_key]
    ][
      [scenario_name]
    ] <- list(

      OBD = OBD_find(
        para,
        scenario
      ),

      sim_num = sim,

      T_prior_blrm = priors$T_prior,

      E_prior_blrm = priors$E_prior,

      scenario = scenario,

      prior_id = p,

      data = generated_data[
        [scenario_name]
      ]
    )
  }
}


# ============================================================
# 9. Model combinations
# ============================================================

combinations <- list(

  "3+3" = list(
    stage1_model = dose3P3,
    stage2_model = NULL
  ),

  "BF_BOIN" = list(
    stage1_model = Backfill.UBOIN.stage1,
    stage2_model = NULL
  ),

  "BF_BOIN+UBOIN" = list(
    stage1_model = Backfill.UBOIN.stage1,
    stage2_model = UBOIN.stage2
  ),

  "BF_BOIN+UBLRM" = list(
    stage1_model = Backfill.UBOIN.stage1,
    stage2_model = UBLRM.stage2
  )
)


# ============================================================
# 10. Prior label
# ============================================================

get_prior_label <- function(prior_id) {

  labels <- c(
    "1" = "strong",
    "2" = "weak",
    "3" = "non"
  )

  label <- labels[
    as.character(prior_id)
  ]

  if (is.na(label)) {
    label <- paste0(
      "prior",
      prior_id
    )
  }

  label
}


# ============================================================
# 11. Run one scenario
# ============================================================

process_scenario <- function(
  scenario_name,
  prior_id,
  sim_objects
) {

  prior_key <- paste0(
    "prior",
    prior_id
  )

  prior_label <- get_prior_label(
    prior_id
  )

  sim_obj <- sim_objects[
    [prior_key]
  ][
    [scenario_name]
  ]

  sim_data <- sim_obj$data

  scenario_results <- data.frame()


  for (
    combined_model_name
    in names(combinations)
  ) {

    stage1_model <-
      combinations[
        [combined_model_name]
      ]$stage1_model

    stage2_model <-
      combinations[
        [combined_model_name]
      ]$stage2_model


    message(
      "Scenario: ", scenario_name,
      " | Model: ", combined_model_name,
      " | Prior: ", prior_label
    )


    result <- run_simulation_stage(
      para,
      sim_data,
      sim_obj,
      stage1_model,
      stage2_model
    )


    output <- check_performance(
      para,
      sim_obj,
      result
    )


    result_row <- data.frame(

      Scenario = scenario_name,

      Model = combined_model_name,

      Prior = prior_label,

      PCS = output$PCS,

      TR_mean = output$TR_mean,

      TR_std = output$TR_std,

      early_stop =
        output$early_stop,

      select_OBD_stage1 =
        paste(
          output$stage1_counts,
          collapse = ", "
        ),

      select_OBD_stage2 =
        if (!is.null(
          output$stage2_counts
        )) {
          paste(
            output$stage2_counts,
            collapse = ", "
          )
        } else {
          "NULL"
        },

      BF_num =
        output$BF_num,

      n_stage1 =
        output$stage1_n,

      n_stage2 =
        if (!is.null(
          output$stage2_n
        )) {
          output$stage2_n
        } else {
          "NULL"
        },

      TR_per_dose =
        paste(
          output$TR_per_dose,
          collapse = ", "
        ),

      OBD =
        paste(
          sim_obj$OBD,
          collapse = ", "
        ),

      stringsAsFactors = FALSE
    )


    scenario_results <- rbind(
      scenario_results,
      result_row
    )
  }

  scenario_results
}


# ============================================================
# 12. Run all simulations
# ============================================================

start_time <- Sys.time()

scenario_set <- names(scenarios)

all_results <- data.frame()

result_filename <- paste0(
  sim,
  "_Test2.csv"
)


for (scenario_name in scenario_set) {

  for (prior_id in prior_ids) {

    scenario_results <-
      process_scenario(
        scenario_name,
        prior_id,
        sim_objects
      )

    all_results <- rbind(
      all_results,
      scenario_results
    )

    # Save progress after each run
    write.csv(
      all_results,
      result_filename,
      row.names = FALSE
    )
  }
}


end_time <- Sys.time()

cat(
  "Total execution time:",
  format(end_time - start_time),
  "\n"
)