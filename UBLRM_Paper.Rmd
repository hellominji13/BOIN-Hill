---
title: "R Notebook"
output: html_notebook
---

Setting
Parameter 설정 및 실제 시나리오 (랜덤) 생성

# setting 1
```{r}
library(MASS)
library(shiny)
library(rhandsontable)
library(LaplacesDemon)
library(ggplot2)
library(readxl)
library(rjags)
library(future)
library(doParallel)
library(foreach)
library(parallel)
library(BayesLogit)

source("./Function_ver2/performance.R")
source("./Function_ver2/Scenario.R") 
source("./Function_ver2/Make_data.R") 
source("./Function_ver2/load_dependencies.R")
source("./Function_ver2/Plot.R")
source("./Function_ver2/UBLRM_Prior.R")
source("./Function_ver2/UBLRM_Model.R") 
source("./Function_ver2/Stage2_Run.R") 
source("./Function_ver2/Stage1_Run.R") 

library(VGAM)
library(BOIN)
library(UBCRM)
library(reshape2)
library(dplyr)
library(readr)

set.seed(42)

#앞으로 코딩에 필요한 parameter (사전설정해야함) 미리 정의
para <- list(
  seed=123,
  
  J = 5, 
  S1 = 12, # back fill 에서는 사용하지 않음 
  S2 = 24,
  N = 54,
  cohort = 3,
  
  window = 1,
  #12,9,30
  n_cap = 12, #back fill dose closed 기준 (전체데이터 기준)
  n_stop = 9, #stay + 현재 용량 유지 이면 용량 증가 끝 (전체데이터 기준)
  N_esc = 30, #용량 증가에 대한 환자 총 수 (=S1)

  DLT = 0.25, #0.25(UBOIN), 0.4
  T_admissible = 0.3, #0.3(UBOIN), #BOIN 에서 자동으로 + 0.05로 설정되어진다 , 0.45
  E_admissible = 0.2,
  CT = 0.95,                
  CE = 0.9
  )

```


```{r}

prior_table <- read_csv("파라미터_UBOIN.csv")

get_priors <- function(prior_table, scenario_name, candidate_num = 1) {
  # 컬럼 이름 구성
  mean_col <- paste0("Cand", candidate_num, "_mean(min,a)")
  sd_col  <- paste0("Cand", candidate_num, "_sd(max,b)")

  # 해당 시나리오 필터
  sub_df <- filter(prior_table, Scenario == scenario_name)

  # prior list
  E_prior <- list()
  T_prior <- list()

  for (i in seq_len(nrow(sub_df))) {
    param <- sub_df$Param[i]
    dist  <- sub_df$Dist[i]
    mean_val <- sub_df[[mean_col]][i]
    sd_val  <- sub_df[[sd_col]][i]

    ## Efficacy 파라미터 처리
    if (param %in% c("delta0", "delta1", "b0", "b1")) {
      if (dist == "Beta") {
        # Beta(shape1 = a, shape2 = b)
        # NOTE: 이 때는 mean = a/(a+b), var = ab/((a+b)^2*(a+b+1)) 이므로 직접 변환 안 함
        E_prior[[paste0(param, "_a")]] <- mean_val
        E_prior[[paste0(param, "_b")]] <- sd_val
      } else if (dist == "Normal") {
        E_prior[[paste0(param, "_mean")]] <- mean_val
        E_prior[[paste0(param, "_sd")]] <- sd_val
      }
    }

    ## Toxicity 파라미터 처리
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

  return(list(E_prior = E_prior, T_prior = T_prior))
}

```

#여러 시나리오에 대해서 
```{r}
source("./Function_ver2/Scenario.R")

file_path <- "시나리오_UBOIN.csv"
scenario_data <- read_csv(file_path)
score <- c(30, 100, 0, 30)
scenarios <- list()
score_matrix <- matrix(score, nrow = 2, byrow = TRUE, dimnames = list("Toxicity" = c("0", "1(occur)"),"Efficacy" = c("0", "1(occur)")))
cat("Score Table (2x2):\n")

for (i in seq(2, nrow(scenario_data), by = 2)) { 
  scenario_name <- as.character(scenario_data[i, 1])
  pi_T <- as.numeric(scenario_data[i, 3:ncol(scenario_data)])
  pi_E <- as.numeric(scenario_data[i+1, 3:ncol(scenario_data)])
  model <- "Gumbel"
  coff <- 1 
  
  scenarios[[scenario_name]] <- list(pi_T = pi_T, pi_E = pi_E, coff = coff, score=score)
}

# para 설정
para$J <- ncol(scenario_data) - 2
para$doses <- as.numeric(scenario_data[1, 3:ncol(scenario_data)])
para$ref_dose <- as.numeric(scenario_data[1, ncol(scenario_data)])

# T_prior_blrm <- Make_prior(para, "T", TRUE)
# E_prior_blrm <- Make_prior(para, "E", TRUE)
# save(T_prior_blrm, E_prior_blrm, file = "prior_blrm_data2.RData")

#pdf("시나리오.pdf", width = 9, height = 4)
par(mfrow = c(2, 4))
for (i in seq_along(scenarios)) {
  scenario_name <- names(scenarios)[i]
  scenario <- scenarios[[scenario_name]]
  prob <- calculate_joint_probabilities(scenario)
  score_prob <- t(apply(prob, 1, function(row) row * scenario$score))
  row_sums <- rowSums(score_prob)
  plot_scenario(para, scenario, name = scenario_name, utility = row_sums, show_legend = (i == 1))

}



```

```{r}
source("./Function_ver2/UBLRM_Model.R") 
source("./Function_ver2/Backfill_Stage1_Run.R") 
source("./Function_ver2/Stage1_Run.R") 
source("./Function_ver2/UBLRM_prior.R") 
source("./Function_ver2/Stage2_RUN.R") 
source("./Function_ver2/Performance.R") 
source("./Function_ver2/Plot.R") 


 
sim <- 10
prior_ids <- 1:2
# Generate all data upfront
sim_objects  <- setNames(vector("list", length(prior_ids)), paste0("prior", prior_ids))
scenario_data <- setNames(vector("list", length(prior_ids)), paste0("prior", prior_ids))


for (scenario_name in names(scenarios) ) { # 
  scenario <- scenarios[[scenario_name]]
  scenario_data[[scenario_name]] <- generate_toxicity_efficacy_time(scenario,sim * 100)}

for (p in prior_ids) {
  prior_key <- paste0("prior", p)
  sim_objects[[prior_key]]  <- setNames(vector("list", length(scenarios)), names(scenarios))
  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    priors <- get_priors(prior_table, scenario_name, candidate_num = p)
    T_prior_blrm <- priors$T_prior
    E_prior_blrm <- priors$E_prior

      sim_objects[[prior_key]][[scenario_name]] <- list(
      OBD = OBD_find(para, scenario),
      sim_num = sim,
      T_prior_blrm = T_prior_blrm,
      E_prior_blrm = E_prior_blrm,
      scenario = scenario,
      prior_id = p,
      # 이미 만들어둔 scenario_data 참조
      data = scenario_data[[scenario_name]]
    )
  }
}

.get_prior_label <- function(prior_id) {
  lab <- c("1" = "strong", "2" = "weak", "3" = "non")[as.character(prior_id)]
  if (is.na(lab)) lab <- paste0("prior", prior_id)
  lab
}


process_scenario <- function(scenario_name, prior_id, sim_objects) {

  prior_key   <- paste0("prior", prior_id)       # "prior1" / "prior2" / "prior3"
  prior_label <- .get_prior_label(prior_id)      # "strong inf" / "weak inf" / "non inf"

  sim_obj <- sim_objects[[prior_key]][[scenario_name]]
  data <- sim_obj$data


  sim <- sim_obj
  scenario_results <- data.frame()

  for (combined_model_name in names(combinations)) {
    stage1_model <- combinations[[combined_model_name]]$stage1_model
    stage2_model <- combinations[[combined_model_name]]$stage2_model

    message("시나리오 ", scenario_name,
            " | 모델 ", combined_model_name,
            " | prior ", prior_key, " (", prior_label, ") 실행 중...")

    result <- run_simulation_stage(
      para, data, sim_obj, stage1_model, stage2_model)
    output <- check_performance(para, sim_obj, result)

    result_row <- data.frame(
      Scenario = scenario_name,
      Model    = combined_model_name,
      Prior    = prior_label,
      PCS      = output$PCS,
      TR_mean  = output$TR_mean,
      TR_std   = output$TR_std,
      early_stop = output$early_stop,
      select_OBD_stage1 = paste(output$stage1_counts, collapse = ", "),
      select_OBD_stage2 = if (!is.null(output$stage2_counts))
                            paste(output$stage2_counts, collapse = ", ")
                          else "NULL",
      BF_num = output$BF_num,
      n_stage1 = output$stage1_n,
      n_stage2 = if (!is.null(output$stage2_n)) output$stage2_n else "NULL",
      TR_per_dose = paste(output$TR_per_dose, collapse = ", "),
      OBD = paste(sim$OBD, collapse = ", "),
      stringsAsFactors = FALSE
    )

    message(paste(result_row, collapse = " "), "\n")
    scenario_results <- rbind(scenario_results, result_row)
  }

  return(scenario_results)
}
```



```{r}
source("./Function_ver2/UBLRM_Model.R") 
source("./Function_ver2/Backfill_Stage1_Run.R") 
source("./Function_ver2/Stage1_Run.R") 
source("./Function_ver2/UBLRM_prior.R") 
source("./Function_ver2/Stage2_RUN.R") 
source("./Function_ver2/Performance.R") 
source("./Function_ver2/Plot.R") 

start_time <- Sys.time()

combinations <- list(
  "3+3" = list(stage1_model = dose3P3, stage2_model = NULL),
  "BF_BOIN" = list(stage1_model = Backfill.UBOIN.stage1, stage2_model = NULL),
  "BF_BOIN+UBOIN" = list(stage1_model = Backfill.UBOIN.stage1, stage2_model = UBOIN.stage2),
  "BF_BOIN+UBLRM" = list(stage1_model = Backfill.UBOIN.stage1, stage2_model = UBLRM.stage2)
  )


scenario_set <-names(scenarios)
prior_ids <- 2:1
start_time <-    Sys.time()

all_results <- data.frame()
result_filename <- paste0(sim, "_Test2.csv")

for (scenario_name in scenario_set) {
  for (prior_id in prior_ids) {#scenario_set,prior_ids
    scenario_results <- process_scenario(scenario_name, prior_id, sim_objects)
    all_results <- rbind(all_results, scenario_results)
    write.csv( all_results, result_filename, row.names = FALSE)
  }
}

end_time <- Sys.time() 
cat("Total execution time:", format(end_time - start_time), "\n")
```

