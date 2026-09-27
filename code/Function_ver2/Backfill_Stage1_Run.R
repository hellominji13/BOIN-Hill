


UBOIN_dose_assgin <- function(para, data, lastdose) {
  boundary <- get.boundary(target = para$DLT, ncohort = para$N*para$J, cohortsize = 1 ,cutoff.eli = para$CT)
  dose_range <- 1:para$J

  esc_label <- "Escalate if # of DLT <="
  de_label  <- "Deescalate if # of DLT >="
  eli_label <- "Eliminate if # of DLT >="
  n_label   <- "Number of patients treated"
  J <- para$J
  r <- lastdose
  m_r <- data$n10[r] + data$n11[r]
  n_r <- data$npt[r]
  
  if (m_r <= boundary$boundary_tab[esc_label, n_r]) {
    nextdose <- min(r + 1, J)
    
  } else if (m_r >= boundary$boundary_tab[de_label, n_r]) {
    
    nextdose <- max(r - 1, 1)
    
    eli_bd_r <- boundary$boundary_tab[eli_label, n_r]
    if (is.na(eli_bd_r)) {
      eli_bd_r <- boundary$boundary_tab[n_label, 1]
    }
    
    if (m_r >= eli_bd_r && n_r >= para$cohort) {
      # r 이상 dose 제거 → 그 아래만 사용
      allowed <- which(1:J < r)
      if (length(allowed) == 0) nextdose <- NA
    }
    
  } else {
    nextdose <- r
  }
  
  pooled_values <- numeric(length(data$n10))
  
  for (i in 1:length(data$n10)) {  
    pooled_values[i] <- if (i < length(data$n10)) {
      (data$n10[i] + data$n11[i] + data$n10[i+1] + data$n11[i+1]) / (data$npt[i] + data$npt[i+1])
    } else {
      (data$n10[i] + data$n11[i])/data$npt[i]
    }
  }
  
  eff_doses <- which((data$n01 + data$n11) > 0)
  if (length(eff_doses) == 0) {
    max_eff_dose <- 100   # 또는 1로 해도 됨 (원하는 구조에 맞춰 변경)
  } else {
    max_eff_dose <- max(eff_doses)
  }
  backfill_candidates <- which(
    data$npt <= para$n_cap &
      seq_along(para$doses) < nextdose &  # next_dose보다 작은 dose
      seq_along(para$doses) >= max_eff_dose &         # n01 또는 n11이 1 이상
      ((data$n10 + data$n11)/data$npt <= boundary$lambda_d |  
         pooled_values <= boundary$lambda_d))
  backfill_dose <- if (length(backfill_candidates) > 0) max(backfill_candidates) else NA
  
  
  get_decision <- function(d) {
    m <- data$n10[d] + data$n11[d]
    n <- data$npt[d]
    if (n <= 0) return(NA_character_)
    
    esc_bd <- boundary$boundary_tab[esc_label, n]
    de_bd  <- boundary$boundary_tab[de_label, n]
    eli_bd <- boundary$boundary_tab[eli_label, n]
    
    if (is.na(eli_bd)) {
      eli_bd <- boundary$boundary_tab[n_label, 1]
    }
    
    if (m <= esc_bd) {
      return("E")  # Escalate
    } else if (m >= de_bd) {
      return("D")  # De-escalation / Elimination
    } else {
      return("S")  # Stay
    }
  }
  
  is_conflict_pair <- function(r_dec, j_dec) {
    if (is.na(j_dec)) return(FALSE)
    # Table 3:
    # Regular Stay vs Backfill Escalation -> conflict
    if (j_dec == "S" && r_dec == "E") return(TRUE)
    # Regular De-escalation/Elimination vs (E, S, D) -> 모두 conflict
    if (j_dec == "D") return(TRUE)
    return(FALSE)
  }
  
  reg_decision <- get_decision(r)
  all_decisions <- sapply(1:J, get_decision)
  candidate_js  <- which(1:J < r)
  conflict_js   <- candidate_js[
    vapply(candidate_js,
           function(j) is_conflict_pair(reg_decision, all_decisions[j]),
           logical(1))
  ]

  # conflict가 없으면 여기서 종료
  if (length(conflict_js) == 0 ) {
    return(list(next_dose    = nextdose,
                backfill_dose = backfill_dose))
  }
  
  # 가장 높은 j* (conflict를 일으키는 backfill dose 중 최고)
  j_star <- max(conflict_js)
  y_vec <- data$n10 + data$n11
  n_vec <- data$npt
  
  pooled_yn <- function(a, b) {
    idx <- a:b
    c(Y = sum(y_vec[idx]), N = sum(n_vec[idx]))
  }
  
  # [j*, r] 구간의 Y, N
  pooled_jr <- pooled_yn(j_star, r)
  Y_jr <- pooled_jr["Y"]
  N_jr <- pooled_jr["N"]
  
  if (N_jr > 0) {
    esc_bd_jr <- boundary$boundary_tab[esc_label, N_jr]
    de_bd_jr  <- boundary$boundary_tab[de_label, N_jr]
    
    if (Y_jr <= esc_bd_jr) {
      ## Escalate: r + 1
      nextdose <- min(r + 1, J)
      
    } else if (Y_jr >= de_bd_jr) {
      ## De-escalate:
      ## [j*, r-1] 중에서 [k, r] pooled가 "안전"한 (Y_kr <= de_bd_kr) 가장 높은 k
      k_candidates <- j_star:(r - 1)
      k_valid <- NA_integer_
      for (k in k_candidates) {
        pooled_kr <- pooled_yn(k, r)
        Y_kr <- pooled_kr["Y"]
        N_kr <- pooled_kr["N"]
        if (N_kr == 0) next
        if (Y_kr <= boundary$boundary_tab[de_label, N_kr]) {
          k_valid <- k
        }
      }
      if (!is.na(k_valid)) {
        nextdose <- k_valid
      } else {
        # 그런 k가 없으면 j* - 1으로 이동
        nextdose <- max(j_star - 1, 1)
      }
      
    } else {
      ## Stay: r 유지
      nextdose <- r
    }
  }
  return(list(next_dose=nextdose, backfill_dose=backfill_dose))
}


Backfill.UBOIN.stage1 <- function(para,outcome_table,time_table,T_prior){
  index <- 0
  esclation_trace <- NULL
  valid_backfill_trace <- NULL
  backfill_trace <- NULL
  all_trace <- NULL
  all_data <- CreData_stage2(para$J)
  data <- CreData_stage2(para$J)
  ovserved_data <- CreData_stage2(para$J) 
  next_dose <- 1
  enter_times <- NULL 
  response_tox_times <- NULL
  response_eff_times <- NULL
  
  current_time <- 0 
  window_time <- 0 
  backfill_dose <- NA
  BF_num <- 0 
  while (next_dose %in% 1:para$J) {
    
    #환자 할당
    current_time <- current_time 
    enter_times <- c(enter_times, rep(current_time, para$cohort))
    index <-  tail(index, 1) + (1:para$cohort)
    temp <- paste(paste("D", next_dose, "S", index, sep = ""), collapse = " ")
    esclation_trace <- paste(esclation_trace, temp, collapse = " ")
    all_trace <- paste(all_trace, temp, collapse = " ")
    times_tox <- time_table[paste("Toxicity_Dose", next_dose, sep = "_"), index, drop = FALSE]
    times_eff <- time_table[paste("Efficacy_Dose", next_dose, sep = "_"), index, drop = FALSE]
    response_tox_times <- c(response_tox_times,times_tox)
    response_eff_times <- c(response_eff_times,times_eff)
    
    current_time <- current_time + para$window
    # Backfill
    # 종료 3
    if (!is.na(backfill_dose)){
      while (TRUE) {
        enter <- tail(enter_times, 1) +rexp(1, rate = 3 / para$window) #한달에 3명 푸아송 
        if (enter >= current_time) {break}
        index <- tail(index, 1)  + 1
        BF_num <- BF_num+1
        enter_times <- c(enter_times, enter)
        temp <- paste(paste("D",backfill_dose, "S", index, sep = ""), collapse = " ")
        backfill_trace <- paste(backfill_trace, temp, collapse = " ")
        all_trace <- paste(all_trace, temp, collapse = " ")
        times_tox <- time_table[paste("Toxicity_Dose", backfill_dose, sep = "_"), index, drop = FALSE]
        times_eff <- time_table[paste("Efficacy_Dose", backfill_dose, sep = "_"), index, drop = FALSE]
        response_tox_times <- c(response_tox_times,times_tox)
        response_eff_times <- c(response_eff_times,times_eff)
      }
    }
    # valid 한 backfill 데이터만 update 하기 
    valid_backfill_trace <- valid_data_trace_stage1(time_table,backfill_trace,enter_times,current_time)
    esclation_data <- extract_responses(esclation_trace, outcome_table)
    vaild_backfill_data <- extract_responses(valid_backfill_trace, outcome_table)
    
    # 아직 안 나온 데이터 예측하기
    if (!is.null(backfill_trace)) {
      pred_backfill_trace <- if (is.null(valid_backfill_trace)) backfill_trace else sub(valid_backfill_trace, "", backfill_trace, fixed = TRUE)
    }
    else{pred_backfill_trace<-NULL}
    ovserved_data[, -1] <- esclation_data[, -1] + vaild_backfill_data[, -1]
    #정확한 값을 출력하지는 않고 예상값 (n10 에 다 몰아서) 를 출력 
    pred_backfill_data <- pred_data(para,time_table,pred_backfill_trace,ovserved_data,enter_times,current_time)
    data[, -1] <- esclation_data[, -1] + vaild_backfill_data[, -1] + pred_backfill_data[,-1]

    #종료 1
    if ((sum(esclation_data$npt) >= para$N_esc)){break}
    
    dose <- UBOIN_dose_assgin(para, data, next_dose)
    backfill_dose <- dose$backfill_dose
    if (is.na(dose$next_dose)){ break}
    if (next_dose==dose$next_dose && data$npt[dose$next_dose] >= para$n_stop){ break}
    next_dose <- dose$next_dose
  }
  # 종료 후 마지막 결과 정리 
  all_data <- esclation_data[, -1] + extract_responses(backfill_trace, outcome_table)[, -1] 
  result <- select.mtd(target = para$DLT, npts = all_data$npt, ntox = (all_data$n10+all_data$n11),cutoff.eli=para$CT )
  mtd <- result$MTD
  

  #early stop 경우 정리
  if (mtd==99){
    print("early stop")
    mtd <- NA 
  }
  return(return(list(obd = mtd,BF_num=BF_num, n = max(index), trace = all_trace , esc_trace = esclation_trace, backfill_trace=backfill_trace, data = all_data, time = current_time +1 , enter_times=enter_times, response_tox_times=response_tox_times, response_eff_times=response_eff_times)))
}




valid_data_trace_stage1 <- function(time_table,trace,enter_times,current_time){
  
  valid_backfill_trace <- NULL
  if (!is.null(trace)) {
    trace_elements <- unlist(strsplit(trace, " "))
    trace_elements <- trace_elements[trace_elements != ""]
    for (trace_element in trace_elements) {
      dose_index <- as.numeric(sub("D(\\d+)S.*", "\\1", trace_element))  
      subject_index <- as.numeric(sub("D\\d+S(\\d+)", "\\1", trace_element))
      backfill_response_time <- time_table[paste("Toxicity_Dose", dose_index, sep = "_"), subject_index, drop = FALSE]
      adjusted_times <- enter_times[subject_index] + backfill_response_time  
      # 반응이 나왔으면 update 
      if (adjusted_times <= current_time){
        valid_backfill_trace <- paste(valid_backfill_trace, trace_element, collapse = " ")
      }
    }
  }
  return(valid_backfill_trace)
}




pred_data <- function(para, time_table, trace, data, enter_times, current_time){
  
  pred_backfill_data <- CreData_stage2(para$J)
  if (is.null(trace) || trace == "") {
    return(pred_backfill_data)  # trace가 없으면 원래 데이터 반환
  }
  
  trace_elements <- unlist(strsplit(trace, " "))
  trace_elements <- trace_elements[trace_elements != ""]
  
  dose_subjects <- vector("list", para$J)
  for (trace_element in trace_elements) {
    dose_index <- as.numeric(sub("D(\\d+)S.*", "\\1", trace_element))
    subject_index <- as.numeric(sub("D\\d+S(\\d+)", "\\1", trace_element))
    
    dose_subjects[[dose_index]] <- c(dose_subjects[[dose_index]], subject_index)
  }
  
  # 베이지안 사후 평균 계산
  s <- data$n10 + data$n11  # 독성을 경험한 환자 수
  r <- data$npt             # 총 평가가 완료된 환자 수
  alpha <- 0 # 0.5 * para$DLT
  beta <- 0 #1 - para$DLT
  p <- (s + alpha) / (r + alpha + beta)
  
  # 각 용량 수준별로 계산
  for (dose_idx in 1:nrow(pred_backfill_data)) {
    if (!is.null(dose_subjects[[dose_idx]])) {
      # 해당 dose의 subjects 가져오기
      subjects <- dose_subjects[[dose_idx]]
      STFT <- 0 
      m <- length(subjects)  # 미완료 환자 수
      
      for (sbj in subjects) {
        backfill_response_time <- time_table[paste("Toxicity_Dose", dose_idx, sep = "_"), sbj, drop = TRUE]
        adjusted_time <- enter_times[sbj] + backfill_response_time
        # STFT 계산 - 표준화된 추적 시간
        STFT <- STFT + ((adjusted_time- current_time) / para$window)
      }
      
      # TITE-BOIN 공식에 따른 예상 독성 수 계산
      tox <-(p[dose_idx] / (1 - p[dose_idx])) * (m - STFT)
      
      # 결과 저장
      pred_backfill_data[dose_idx, "npt"] <- m
      pred_backfill_data[dose_idx, "n10"] <- min(m,tox)  # 예상 독성 수
    }
  }
  
  return(pred_backfill_data)
}



