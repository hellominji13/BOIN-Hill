# ============================================================
# UBLRM_Paper.R
#
# Main simulation script for comparison of:
#   1. 3+3
#   2. BF_BOIN
#   3. UBOIN
#   4. BOIN_Hill
#
# Trial settings, scenarios, and prior specifications
# are defined in config.R.
#
# Input files are read in config.R:
#   ../data/시나리오_UBOIN.csv
#   ../data/파라미터_UBOIN.csv
#
# Output files are written to:
#   ../results/
# ============================================================



# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "MASS",
  "shiny",
  "rhandsontable",
  "LaplacesDemon",
  "ggplot2",
  "readxl",
  "rjags",
  "future",
  "doParallel",
  "foreach",
  "parallel",
  "BayesLogit",
  "VGAM",
  "BOIN",
  "UBCRM",
  "reshape2",
  "dplyr",
  "readr",
  "copula"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0) {

  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", ")
  )
}


suppressPackageStartupMessages(
  invisible(
    lapply(
      required_packages,
      library,
      character.only = TRUE
    )
  )
)



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

  source(
    file.path(
      "Function_ver2",
      file
    )
  )
}



# ============================================================
# 2. Load simulation configuration
# ============================================================

# config.R contains:
#
#   sim
#     Number of simulation replicates
#
#   para
#     Trial-design parameters
#
#   prior_table
#     Prior specifications read from:
#       ../data/파라미터_UBOIN.csv
#
#   scenarios
#     True toxicity / efficacy scenarios read from:
#       ../data/시나리오_UBOIN.csv
#
#   boin_hill_prior_ids
#     Prior specifications evaluated for BOIN_Hill
#     e.g. c(1, 2)
#
#   default_prior_id
#     Placeholder prior used by methods that do not
#     use the BOIN_Hill prior
#
source("config.R")


# Defaults, in case these are not explicitly provided
if (!exists("boin_hill_prior_ids")) {
  boin_hill_prior_ids <- c(1, 2)
}


if (!exists("default_prior_id")) {
  default_prior_id <- boin_hill_prior_ids[1]
}


set.seed(para$seed)



# ============================================================
# 3. Extract BOIN-Hill priors
# ============================================================

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


  if (nrow(sub_df) == 0) {

    stop(
      "No prior specification found for scenario: ",
      scenario_name
    )
  }


  if (!(mean_col %in% names(sub_df))) {

    stop(
      "Prior column not found: ",
      mean_col
    )
  }


  if (!(sd_col %in% names(sub_df))) {

    stop(
      "Prior column not found: ",
      sd_col
    )
  }


  E_prior <- list()
  T_prior <- list()


  for (i in seq_len(nrow(sub_df))) {

    param <- sub_df$Param[i]
    dist  <- sub_df$Dist[i]

    mean_val <- sub_df[[mean_col]][i]
    sd_val   <- sub_df[[sd_col]][i]


    # --------------------------------------------------------
    # Efficacy-model priors
    # --------------------------------------------------------

    if (
      param %in%
        c(
          "delta0",
          "delta1",
          "b0",
          "b1"
        )
    ) {

      if (dist == "Beta") {

        E_prior[[paste0(param, "_a")]] <- mean_val

        E_prior[[paste0(param, "_b")]] <- sd_val

      } else if (dist == "Normal") {

        E_prior[[paste0(param, "_mean")]] <- mean_val

        E_prior[[paste0(param, "_sd")]] <- sd_val
      }
    }


    # --------------------------------------------------------
    # Toxicity-model priors
    # --------------------------------------------------------

    if (
      param %in%
        c(
          "a0",
          "a1"
        )
    ) {

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
# 4. Create results directory
# ============================================================

results_dir <- file.path(
  "..",
  "results"
)


if (!dir.exists(results_dir)) {

  dir.create(
    results_dir,
    recursive = TRUE
  )
}



# ============================================================
# 5. Plot true scenarios
# ============================================================

scenario_plot_file <- file.path(
  results_dir,
  "scenario_plot.pdf"
)


pdf(
  scenario_plot_file,
  width = 12,
  height = 6
)


par(
  mfrow = c(2, 4)
)


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
      function(row) {
        row * scenario$score
      }
    )
  )


  utility <- rowSums(
    score_prob
  )


  plot_scenario(
    para,
    scenario,
    name = scenario_name,
    utility = utility,
    show_legend = (i == 1)
  )
}


dev.off()



# ============================================================
# 6. Generate simulation data
# ============================================================
#
# A single simulated dataset is generated for each scenario
# and shared across all methods and prior settings.
#
# This ensures that methods are compared using the same
# simulated patient outcomes.
#
# `sim * 100` is retained from the original simulation code.
# ============================================================

generated_data <- setNames(
  vector(
    "list",
    length(scenarios)
  ),
  names(scenarios)
)


for (
  scenario_name
  in names(scenarios)
) {

  generated_data[[scenario_name]] <-

    generate_toxicity_efficacy_time(
      scenarios[[scenario_name]],
      sim * 100
    )
}



# ============================================================
# 7. Build simulation objects
# ============================================================
#
# BOIN_Hill uses two prior specifications:
#
#   prior 1 = strong
#   prior 2 = weak
#
# The same simulation objects can also be used by
# 3+3, BF_BOIN, and UBOIN. Those methods are run only once,
# using default_prior_id as a placeholder because they do not
# depend on the BOIN-Hill prior specification.
# ============================================================

prior_ids <- boin_hill_prior_ids


sim_objects <- setNames(
  vector(
    "list",
    length(prior_ids)
  ),
  paste0(
    "prior",
    prior_ids
  )
)


for (p in prior_ids) {

  prior_key <- paste0(
    "prior",
    p
  )


  sim_objects[[prior_key]] <- setNames(

    vector(
      "list",
      length(scenarios)
    ),

    names(scenarios)
  )


  for (
    scenario_name
    in names(scenarios)
  ) {

    scenario <- scenarios[[scenario_name]]


    priors <- get_priors(
      prior_table = prior_table,
      scenario_name = scenario_name,
      candidate_num = p
    )


    sim_objects[[prior_key]][[scenario_name]] <- list(

      OBD = OBD_find(
        para,
        scenario
      ),

      sim_num = sim,

      T_prior_blrm =
        priors$T_prior,

      E_prior_blrm =
        priors$E_prior,

      scenario =
        scenario,

      prior_id =
        p,

      data =
        generated_data[[scenario_name]]
    )
  }
}



# ============================================================
# 8. Methods
# ============================================================
#
# Display names used in the final output:
#
#   3+3
#   BF_BOIN
#   UBOIN
#   BOIN_Hill
#
# Internal function names are kept unchanged.
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


  "UBOIN" = list(
    stage1_model = Backfill.UBOIN.stage1,
    stage2_model = UBOIN.stage2
  ),


  "BOIN_Hill" = list(
    stage1_model = Backfill.UBOIN.stage1,
    stage2_model = UBLRM.stage2
  )
)



# ============================================================
# 9. Prior labels
# ============================================================

get_prior_label <- function(
  prior_id
) {

  labels <- c(
    "1" = "strong",
    "2" = "weak"
  )


  prior_id_chr <- as.character(
    prior_id
  )


  if (
    prior_id_chr
    %in%
    names(labels)
  ) {

    return(
      unname(
        labels[
          prior_id_chr
        ]
      )
    )
  }


  paste0(
    "prior",
    prior_id
  )
}



# ============================================================
# 10. Run one method for one scenario
# ============================================================

process_scenario <- function(
  scenario_name,
  combined_model_name,
  prior_id,
  sim_objects
) {

  prior_key <- paste0(
    "prior",
    prior_id
  )


  sim_obj <-
    sim_objects[[prior_key]][[scenario_name]]


  sim_data <-
    sim_obj$data


  # ----------------------------------------------------------
  # Prior label
  #
  # BOIN_Hill:
  #   prior 1 / prior 2 are actually compared.
  #
  # Other methods:
  #   BOIN-Hill prior is not used.
  # ----------------------------------------------------------

  if (
    combined_model_name
    == "BOIN_Hill"
  ) {

    prior_label <-
      get_prior_label(
        prior_id
      )

  } else {

    prior_label <-
      "Not used"
  }


  # ----------------------------------------------------------
  # Select model
  # ----------------------------------------------------------

  stage1_model <-
    combinations[[combined_model_name]]$stage1_model


  stage2_model <-
    combinations[[combined_model_name]]$stage2_model


  message(
    "Scenario: ",
    scenario_name,
    " | Model: ",
    combined_model_name,
    " | Prior: ",
    prior_label
  )


  # ----------------------------------------------------------
  # Run simulation
  # ----------------------------------------------------------

  result <- run_simulation_stage(
    para,
    sim_data,
    sim_obj,
    stage1_model,
    stage2_model
  )


  # ----------------------------------------------------------
  # Calculate operating characteristics
  # ----------------------------------------------------------

  output <- check_performance(
    para,
    sim_obj,
    result
  )


  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------

  result_row <- data.frame(

    Scenario =
      scenario_name,

    Model =
      combined_model_name,

    Prior =
      prior_label,

    PCS =
      output$PCS,

    TR_mean =
      output$TR_mean,

    TR_std =
      output$TR_std,

    early_stop =
      output$early_stop,

    select_OBD_stage1 =
      paste(
        output$stage1_counts,
        collapse = ", "
      ),

    select_OBD_stage2 =
      if (
        !is.null(
          output$stage2_counts
        )
      ) {

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
      if (
        !is.null(
          output$stage2_n
        )
      ) {

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

    stringsAsFactors =
      FALSE
  )


  return(
    result_row
  )
}



# ============================================================
# 11. Run all simulations
# ============================================================
#
# For every scenario:
#
#   3+3
#       run once
#
#   BF_BOIN
#       run once
#
#   UBOIN
#       run once
#
#   BOIN_Hill
#       run with prior 1
#       run with prior 2
#
# Therefore, each scenario produces five result rows.
# ============================================================

start_time <- Sys.time()


scenario_set <- names(
  scenarios
)


all_results <- data.frame()


result_filename <- file.path(
  results_dir,
  paste0(
    sim,
    "_simulation_results.csv"
  )
)


for (
  scenario_name
  in scenario_set
) {

  for (
    combined_model_name
    in names(combinations)
  ) {


    # --------------------------------------------------------
    # Determine which prior(s) to use
    # --------------------------------------------------------

    if (
      combined_model_name
      == "BOIN_Hill"
    ) {

      # BOIN_Hill is evaluated under both priors
      current_prior_ids <-
        boin_hill_prior_ids

    } else {

      # Other methods do not use the BOIN-Hill prior.
      # Use a single placeholder simulation object.
      current_prior_ids <-
        default_prior_id
    }


    # --------------------------------------------------------
    # Run required prior setting(s)
    # --------------------------------------------------------

    for (
      prior_id
      in current_prior_ids
    ) {

      result_row <-
        process_scenario(

          scenario_name =
            scenario_name,

          combined_model_name =
            combined_model_name,

          prior_id =
            prior_id,

          sim_objects =
            sim_objects
        )


      all_results <- rbind(
        all_results,
        result_row
      )


      # ------------------------------------------------------
      # Save intermediate results after every run
      #
      # This is useful in Code Ocean because partial progress
      # is preserved even if a later simulation fails.
      # ------------------------------------------------------

      write.csv(
        all_results,
        result_filename,
        row.names = FALSE
      )
    }
  }
}



# ============================================================
# 12. Finished
# ============================================================

end_time <- Sys.time()


cat(
  "\nSimulation completed.\n"
)


cat(
  "Number of scenarios:",
  length(scenarios),
  "\n"
)


cat(
  "Simulation replicates:",
  sim,
  "\n"
)


cat(
  "Results saved to:",
  result_filename,
  "\n"
)


cat(
  "Total execution time:",
  format(
    end_time -
      start_time
  ),
  "\n"
)