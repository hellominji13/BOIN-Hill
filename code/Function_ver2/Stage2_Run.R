CreData_stage2 <- function (ndose = 3, dosenames = paste("dose", 1:ndose, sep = " ")) 
{
  data <- data.frame(dose = 1:ndose, npt = rep(0, ndose), 
                     n00 = rep(0, ndose),n01 = rep(0, ndose),n10 = rep(0, ndose),n11 = rep(0, ndose), row.names = dosenames)
  return(data)
}

extract_responses <- function(trace, outputcom) {
  if (is.null(trace) || length(trace)==0) {
    results <- CreData_stage2(para$J)
  }
  else {
    results <- CreData_stage2(para$J)
    input_list <- strsplit(trace, " ")[[1]]
    for (input in input_list) {
      if (input == "") next
      
      # dose와 subject 값을 추출
      dose <- as.numeric(sub("D(\\d+)S.*", "\\1", input))  
      subject <- sub("D\\d+S(\\d+)", "\\1", input)
      
      # Toxicity와 Efficacy 응답값 추출
      toxicity_string <- paste0("Toxicity_Dose_", dose)
      efficacy_string <- paste0("Efficacy_Dose_", dose)
      subject_col <- paste0("subject_", subject)
      
      toxicity_response <- outputcom[toxicity_string, subject_col]
      efficacy_response <- outputcom[efficacy_string, subject_col]
      
      # 반응값 설정
      response_col <- paste0("n", toxicity_response,efficacy_response)
      # 결과값 업데이트
      results[dose, response_col] <- results[dose, response_col] + 1
      results[dose, "npt"] <- results[dose, "npt"] + 1  # npt 카운트 증가
    }
  }
  return(results)
}


############ UBOIN ################
admissible_dose_range_function <- function(data) {
  a_prior <- rep(1/4, 4)
  admissible_dose_range <- c()
  for (j in 1:para$J) {
    sample <- rdirichlet(1000, a_prior + as.numeric(data[j, 3:ncol(data)]))
    
    tox_prob <- sample[, 3] + sample[, 4]  # p10 + p11
    eff_prob <- sample[, 2] + sample[, 4]  # p01 + p11
    
    tox_condition <- mean(tox_prob > para$T_admissible)
    eff_condition <- mean(eff_prob < para$E_admissible)
    
    if (tox_condition <= para$CT && eff_condition <= para$CE) {
      admissible_dose_range <- c(admissible_dose_range, j)
    }
  }
  return(admissible_dose_range)
}


Dose_assign_function <- function(data,score,stage1_dose){
  if(is.na(stage1_dose)){
    return(NA)
  } 
  boundary <- get.boundary(target = para$DLT, ncohort = (para$N / para$cohort), cohortsize = para$cohort)
  a_prior <- rep(1/4, 4)
  #B1
  j = max(which(data$npt != 0))
  tox_prob <-  data$n10 + data$n11 #p10+p11
  if (j > 0 && para$J>j && ((tox_prob[j] / data$npt[j]) <= boundary$lambda_e)) {
    next_dose <- j+1
    return(next_dose)
  }
  #B2
  else {
    admissible_dose_range <- admissible_dose_range_function(data)
    expected_utilities <- rep(-10, para$J)
    if (length(admissible_dose_range) == 0){
      return(NA)
    }
    else{
      for (aj in admissible_dose_range) {
        expected_posterior_toxicity <- (a_prior + as.numeric(data[aj, 3:ncol(data)])) / (data$npt[aj] + sum(a_prior))
        expected_utilities[aj] <- sum(score * expected_posterior_toxicity)
      }
      next_dose <- which.max(expected_utilities) 
      return(next_dose)
    }
  }
}

UBOIN.stage2 <-function(para,outcome_table,stage1,T_prior,E_prior,score){
  stage1_data <- extract_responses(stage1$trace, outcome_table)
  index <- 1:para$cohort
  trace <- NULL
  stage2_data <- CreData_stage2(para$J)
  next_dose <- 1
  #B3  
  while (next_dose %in% 1:para$J) {
    data <- stage1_data
    data[, -1] <- stage1_data[, -1] + stage2_data[, -1]
    if ((max(data$npt) >= para$S2 )|| (sum(data$npt) >= para$N) ) {
      break
    }
    next_dose <- Dose_assign_function(data,score,stage1$obd)
    if (is.na(next_dose)) { break }
    temp <- paste(paste("D", next_dose, "S", index+stage1$n, sep = ""), collapse = " ")
    trace <- paste(trace, temp, collapse = " ")
    index <- index + para$cohort
    stage2_data <- extract_responses(trace, outcome_table)
  }
  
  observed_tox <- sapply(1:nrow(data), function(j) {
    if (data$npt[j] > 0) {
      return((data$n10[j] + data$n11[j]) / data$npt[j])
    } else {
      return(NA)
    }
  })
  observed_eff <- sapply(1:nrow(data), function(j) {
    if (data$npt[j] > 0) {
      return((data$n01[j] + data$n11[j]) / data$npt[j])
    } else {
      return(NA)
    }
  })
  fitting <- disjoint(observed_tox,observed_eff,0)
  return(list(obd = next_dose, n = min(index) - 1, trace = trace, data = stage2_data,
              fitting=fitting, prob_tox=observed_tox, prob_eff=observed_eff))
}




############ BLRM ################





UBLRM_Dose_assign_function <- function(data, T_prior, E_prior,current_dose,score) {
  #stage1, stage2 data 통합
  # print(data)
  prob <- UBLRM_model_fit_function(para,data, T_prior, E_prior)
  admissible <- Admissible_dose_range_function_BRLM(para, prob, current_dose)
  fitting <- disjoint(prob$tox,prob$eff,prob$coff)
  # print(prob$tox)
  # print(prob$eff)
  expected.utility <- colSums(fitting * score)  # apply 대신 벡터연산

  if (length(admissible) == 0) {
    return(list(
      next_dose = NA,
      obd = NA,
      T_prior_blrm = prob$T_prior_blrm,
      E_prior_blrm = prob$E_prior_blrm,
      fitting=fitting
    ))
  }
  else if (length(admissible) == 1) {
    return(list(
      next_dose = admissible,
      obd = admissible,
      T_prior_blrm = prob$T_prior_blrm,
      E_prior_blrm = prob$E_prior_blrm,
      fitting=fitting
    ))
  }
  else {
    obd <- which.max(replace(expected.utility, -admissible, -10))
    return(list(
      next_dose = obd,
      obd = obd,
      T_prior_blrm = prob$T_prior_blrm,
      E_prior_blrm = prob$E_prior_blrm,
      fitting =fitting
    ))
  }
}


UBLRM.stage2 <-function(para,outcome_table,stage1,T_prior,E_prior,score){
  stage1_data <- extract_responses(stage1$trace, outcome_table)
  index <- 1:para$cohort
  trace <- NULL
  dose_range <- stage1_data$dose
  stage2_data <- CreData_stage2(para$J)
  next_dose <- 1
  while (next_dose %in% 1:para$J) {
    data <- stage1_data
    data[, -1] <- stage1_data[, -1] + stage2_data[, -1]
    dose_assign <- UBLRM_Dose_assign_function(data,T_prior,E_prior,stage1$obd,score)
    next_dose <- dose_assign$next_dose
    # print("MINJI")
    # print(next_dose)
    if ((max(data$npt) >= para$S2) || (sum(data$npt) >= para$N)) {break}
    if (is.na(next_dose)) { break }
    temp <- paste(paste("D", next_dose, "S", index+stage1$n, sep = ""), collapse = " ")
    trace <- paste(trace, temp, collapse = " ")
    index <- index + para$cohort
    stage2_data <- extract_responses(trace, outcome_table)
  }
  prob <- UBLRM_model_fit_function(para, data, T_prior, E_prior)
  admissible <- Admissible_dose_range_function_BRLM(para, prob, stage1$obd)
  
  admissible_extra <- which(prob$tox <= para$DLT)  # 독성이 기준 이내
  admissible <- intersect(admissible, admissible_extra)
  fitting <- disjoint(prob$tox, prob$eff, prob$coff)
  expected.utility <- colSums(fitting * score)

  if (length(admissible) == 0 || all(is.na(admissible))) {
    # admissible이 NA이거나 빈 벡터일 때: 전체를 -10 처리
    expected.utility[] <- -10
    obd <- NA
  } else {
    expected.utility[-admissible] <- -10
    obd <- which.max(expected.utility)
  }
  # print("최종")
  # print(expected.utility)
  # print(obd)
  #dose_assign$obd
  
  return(list(obd = obd, n = min(index) - 1, trace = trace, data = stage2_data,fitting=fitting, 
              prob_tox=prob$tox, prob_eff=prob$eff, prob_coff=prob$coff))
}





