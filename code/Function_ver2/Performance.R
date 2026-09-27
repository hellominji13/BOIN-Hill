
run_simulation_stage <- function(para, data, sim, stage1, stage2 = NULL) {
  total_n <- numeric(sim$sim_num)
  save_OBD_stage1 <- numeric(sim$sim_num)
  save_OBD_stage2 <- numeric(sim$sim_num)
  dose_treatments <- matrix(0, nrow = para$J, ncol =sim$sim_num)
  total_n_stage1 <- numeric(sim$sim_num)
  total_n_stage2 <- numeric(sim$sim_num)
  BF_number <- numeric(sim$sim_num)
  
  all_pi_T <- matrix(0, nrow = sim$sim_num, ncol = para$J)
  all_pi_E <- matrix(0, nrow = sim$sim_num, ncol = para$J)
  
  all_out_come <- data$outcome_table
  all_out_come_time <- data$time_table
  outcome_data <- lapply(1:sim$sim_num, function(i) {
    start_col <- ((i-1) * 100) + 1
    end_col <- i * 100
    df <- all_out_come[, start_col:end_col]
    colnames(df) <- colnames(all_out_come)[1:100]  # 열 이름 유지
    return(df)
  })
  time_data <- lapply(1:sim$sim_num, function(i) {
    start_col <- ((i-1) * 100) + 1
    end_col <- i * 100
    df <- all_out_come_time[, start_col:end_col]
    colnames(df) <- colnames(all_out_come_time)[1:100]  # 열 이름 유지
    return(df)
  })
  for (i in 1:sim$sim_num) {
    outcome_table <- outcome_data[[i]]
    time_table <- time_data[[i]]
    result <- stage1(para, outcome_table,time_table,sim$T_prior_blrm)
    dose_treatments[,i] <- result$data["npt"][,]
    total_n_stage1[i] <- result$n
    save_OBD_stage1[i] <- result$obd
    BF_number[i] <- if (!is.null(result$BF_num) && !is.na(result$BF_num)) result$BF_num else 0
    
    # Stage2
    if (!is.null(stage2)) {
      result2 <- stage2(para,outcome_table, result, sim$T_prior_blrm,sim$E_prior_blrm, sim$scenario$score)
      pi_T <- result2$fitting[3,] + result2$fitting[4,]
      pi_E <- result2$fitting[2,] + result2$fitting[4,]
      all_pi_T[i,] <- pi_T
      all_pi_E[i,] <- pi_E 
      dose_treatments[,i] <- dose_treatments[,i]+result2$data["npt"][,]
      total_n_stage2[i] <- result2$n # stage1과 합침
      save_OBD_stage2[i] <- result2$obd #stage1 결과 덮기
    }
  }
  return(list(
    dose_treatments = dose_treatments,
    total_n_stage1 = total_n_stage1,
    total_n_stage2 = if(!is.null(stage2)) total_n_stage2 else NULL,
    OBD_stage1 = save_OBD_stage1,
    OBD_stage2 = if(!is.null(stage2)) save_OBD_stage2 else NULL,
    plot = if(!is.null(stage2)) single_plot_curve2(all_pi_T, all_pi_E, sim$scenario,"none")else NULL,
    estimate_pi_T=all_pi_T, estimate_pi_E=all_pi_E,
    BF_number =BF_number
    
  ))
}




check_performance <- function(para,sim,output) {
  # Stage 1 early stopping check
  OBD_count <- 0
  early_stop <- 0
  OBD_treatd <- numeric(sim$sim_num)
  tot_n <- numeric(sim$sim_num)
  for (i in 1:sim$sim_num) {
    #PCS
    if (is.null(output$OBD_stage2[i])){
      OBD_count <- OBD_count  + (output$OBD_stage1[i] %in% sim$OBD)
      
      tot_n[i] <- output$total_n_stage1[i]
      early_stop <- early_stop + is.na(output$OBD_stage1[i])
    }
    else{
      OBD_count <- OBD_count  + (output$OBD_stage2[i] %in% sim$OBD)
      tot_n[i] <- output$total_n_stage1[i] + output$total_n_stage2[i]
      early_stop <- early_stop + is.na(output$OBD_stage2[i])
    }
    
    #TR
    OBD_treatd[i] <- if(is.na(sim$OBD[1])) 0 else sum(output$dose_treatments[sim$OBD,i])
    
  }
  
  return(list(PCS = OBD_count*100/sim$sim_num,
              TR_std = sd(OBD_treatd*100/tot_n),TR_mean = mean(OBD_treatd*100/tot_n),
              early_stop = early_stop, TR_per_dose = rowMeans(output$dose_treatments, na.rm = TRUE),
              BF_num = mean(output$BF_number),
              stage1_n = mean(output$total_n_stage1), 
              stage2_n = if(!is.null(output$total_n_stage2)) mean(output$total_n_stage2) else NULL,
              stage1_counts = as.numeric((table(factor(output$OBD_stage1, levels = 1:para$J)))),
              stage2_counts =  if(!is.null(output$OBD_stage2)) as.numeric(table(factor(output$OBD_stage2, levels = 1:para$J))) else NULL))
}

# 이름 넣으면 plot하게 만드는 함수 
plot_performance <- function(model_names) {
  # 모든 모델의 성능 지표를 저장할 데이터 프레임 초기화
  results <- data.frame()
  
  # 각 모델의 성능 지표를 데이터 프레임에 추가
  for (model_name in model_names) {
    # 모델 객체를 동적으로 참조
    model_data <- get(model_name)
    
    model_results <- data.frame(
      Model = model_name,
      PCS = model_data$PCS,
      TR_mean = model_data$TR_mean,
      TR_std = model_data$TR_std,
      early_stop = model_data$early_stop
    )
    
    results <- rbind(results, model_results)
  }
  
  # 소수점 둘째 자리까지 반올림
  results <- results %>% mutate(across(where(is.numeric), ~ round(., 2)))
  
  # 데이터 변환: PCS, TR_mean, early_stop을 길게 변환
  results_long <- reshape(
    results,
    varying = c("PCS", "TR_mean","early_stop"),
    v.names = "Value",
    times = c("PCS", "TR_mean","early_stop"),
    timevar = "Metric",
    direction = "long"
  )
  
  # 막대 그래프 생성
  ggplot(results_long, aes(x = Model, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    scale_fill_manual(values = c("PCS" = "blue", "TR_mean" = "red", "early_stop"="green")) +
    labs(title = "Model Performance Metrics",
         x = "Model",
         y = "Value") +
    theme_minimal() +
    theme(legend.title = element_blank())
}

