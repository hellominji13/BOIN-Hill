plot_clinical_timeline <- function(result, result2 = NULL, time_table, para) {
  library(ggplot2)
  
  # 반응 시간 데이터 추출
  get_response_times <- function(trace, enter_times, is_backfill = FALSE) {
    if(is.null(trace) || trace == "") return(data.frame())
    
    response_times <- list()
    
    trace_elements <- unlist(strsplit(trace, " "))
    trace_elements <- trace_elements[trace_elements != ""]
    
    for(trace_element in trace_elements) {
      dose_index <- as.numeric(sub("D(\\d+)S.*", "\\1", trace_element))  
      subject_index <- as.numeric(sub("D\\d+S(\\d+)", "\\1", trace_element))
      
      # Toxicity 체크
      tox_response <- time_table[paste("Toxicity_Dose", dose_index, sep = "_"), subject_index]
      if(tox_response < 1) {
        response_times[[length(response_times) + 1]] <- data.frame(
          subject_index = subject_index,
          dose = dose_index,
          enter_time = enter_times[subject_index],
          response_time = enter_times[subject_index] + tox_response,
          type = "t",
          is_backfill = is_backfill
        )
      }
      
      # Efficacy 체크
      eff_response <- time_table[paste("Efficacy_Dose", dose_index, sep = "_"), subject_index]
      if(eff_response < 1) {
        response_times[[length(response_times) + 1]] <- data.frame(
          subject_index = subject_index,
          dose = dose_index,
          enter_time = enter_times[subject_index],
          response_time = enter_times[subject_index] + eff_response,
          type = "e",
          is_backfill = is_backfill
        )
      }
    }
    
    if(length(response_times) > 0) {
      return(do.call(rbind, response_times))
    } else {
      return(data.frame())
    }
  }
  
  # 메인 데이터와 backfill 데이터 준비
  get_subject_data <- function(trace, is_backfill = FALSE) {
    if(is.null(trace) || trace == "") return(data.frame(
      subject_index = numeric(0),
      dose = numeric(0),
      is_backfill = logical(0),
      time = numeric(0)
    ))
    
    trace_elements <- unlist(strsplit(trace, " "))
    trace_elements <- trace_elements[trace_elements != ""]
    
    data.frame(
      subject_index = as.numeric(sub("D\\d+S(\\d+)", "\\1", trace_elements)),
      dose = as.numeric(sub("D(\\d+)S.*", "\\1", trace_elements)),
      is_backfill = is_backfill
    )
  }
  
  # 데이터 준비 (이전과 동일)
  main_subjects <- get_subject_data(result$esc_trace, FALSE)
  backfill_subjects <- get_subject_data(result$backfill_trace, TRUE)
  
  times_df <- data.frame(
    time = result$enter_times,
    subject_index = seq_along(result$enter_times)
  )
  
  # 겹치는 enter times 처리
  adjust_overlapping_times <- function(plot_data) {
    if(nrow(plot_data) == 0) return(plot_data)
    
    # 시간 순서로 정렬하되, 같은 시간 내에서는 subject_index 순서로 정렬
    plot_data <- plot_data[order(plot_data$time, plot_data$subject_index),]
    result <- plot_data
    
    # 같은 시간에 있는 점들을 찾아서 약간씩 이동
    for(t in unique(plot_data$time)) {
      same_time <- which(plot_data$time == t)
      if(length(same_time) > 1) {
        # offset 크기를 더 작게 조정하고 subject_index 순서대로 적용
        offset <- seq(-0.05, 0.05, length.out = length(same_time))
        result$time[same_time] <- plot_data$time[same_time] + offset[order(plot_data$subject_index[same_time])]
      }
    }
    return(result)
  }
  calculate_y_offsets <- function(responses) {
    if(nrow(responses) == 0) return(responses)
    
    responses <- responses[order(responses$adjusted_enter_time),]
    responses$y_offset <- 0
    
    for(d in unique(responses$dose)) {
      dose_responses <- which(responses$dose == d)
      if(length(dose_responses) > 1) {
        used_offsets <- c()
        for(i in dose_responses) {
          current_range <- c(responses$adjusted_enter_time[i], responses$response_time[i])
          
          # 현재 response와 시간이 겹치는 다른 response 찾기
          overlapping <- sapply(dose_responses[dose_responses != i], function(j) {
            other_range <- c(responses$adjusted_enter_time[j], responses$response_time[j])
            max(current_range[1], other_range[1]) <= min(current_range[2], other_range[2])
          })
          
          if(any(overlapping)) {
            # 사용 중인 offset 찾기
            used <- unique(responses$y_offset[dose_responses[overlapping]])
            # 사용되지 않은 가장 작은 offset 찾기
            available_offsets <- setdiff(seq(0.2, 1, by = 0.2), used)
            if(length(available_offsets) > 0) {
              responses$y_offset[i] <- min(available_offsets)
            } else {
              responses$y_offset[i] <- max(used) + 0.2
            }
          }
        }
      }
    }
    return(responses)
  }
  
  if (!is.null(result2)) {
    # result2의 시간을 result의 마지막 시간부터 시작하도록 조정
    time_shift <- max(result$enter_times)
    result2$enter_times <- result2$enter_times + time_shift
  }
  
  # result 데이터 준비
  main_subjects1 <- get_subject_data(result$esc_trace, FALSE)
  backfill_subjects1 <- get_subject_data(result$backfill_trace, TRUE)
  
  times_df1 <- data.frame(
    time = result$enter_times,
    subject_index = seq_along(result$enter_times)
  )
  
  # result2 데이터 준비
  if (!is.null(result2)) {
    main_subjects2 <- get_subject_data(result2$esc_trace, FALSE)
    backfill_subjects2 <- get_subject_data(result2$backfill_trace, TRUE)
    
    times_df2 <- data.frame(
      time = result2$enter_times,
      subject_index = seq_along(result2$enter_times)
    )
  }
  
  # 데이터 처리
  main_plot_data1 <- if(nrow(main_subjects1) > 0) {
    adjust_overlapping_times(merge(main_subjects1, times_df1))
  } else data.frame()
  
  backfill_plot_data1 <- if(nrow(backfill_subjects1) > 0) {
    adjust_overlapping_times(merge(backfill_subjects1, times_df1))
  } else data.frame()
  
  if (!is.null(result2)) {
    main_plot_data2 <- if(nrow(main_subjects2) > 0) {
      adjust_overlapping_times(merge(main_subjects2, times_df2))
    } else data.frame()
    
    backfill_plot_data2 <- if(nrow(backfill_subjects2) > 0) {
      adjust_overlapping_times(merge(backfill_subjects2, times_df2))
    } else data.frame()
  }
  
  # 반응 데이터 준비
  all_responses1 <- rbind(
    get_response_times(result$esc_trace, result$enter_times, FALSE),
    get_response_times(result$backfill_trace, result$enter_times, TRUE)
  )
  
  if (!is.null(result2)) {
    all_responses2 <- rbind(
      get_response_times(result2$esc_trace, result2$enter_times, FALSE),
      get_response_times(result2$backfill_trace, result2$enter_times, TRUE)
    )
  }
  
  # result1 반응 데이터 처리
  if(nrow(all_responses1) > 0) {
    all_subjects_adjusted1 <- rbind(
      if(nrow(main_plot_data1) > 0) main_plot_data1 else NULL,
      if(nrow(backfill_plot_data1) > 0) backfill_plot_data1 else NULL
    )
    
    all_responses1 <- merge(
      all_responses1,
      all_subjects_adjusted1[, c("subject_index", "time")],
      by = "subject_index"
    )
    all_responses1$adjusted_enter_time <- all_responses1$time
    all_responses1 <- calculate_y_offsets(all_responses1)
  }
  
  # result2 반응 데이터 처리
  if (!is.null(result2) && nrow(all_responses2) > 0) {
    all_subjects_adjusted2 <- rbind(
      if(nrow(main_plot_data2) > 0) main_plot_data2 else NULL,
      if(nrow(backfill_plot_data2) > 0) backfill_plot_data2 else NULL
    )
    
    all_responses2 <- merge(
      all_responses2,
      all_subjects_adjusted2[, c("subject_index", "time")],
      by = "subject_index"
    )
    all_responses2$adjusted_enter_time <- all_responses2$time
    all_responses2 <- calculate_y_offsets(all_responses2)
  }
  
  # Plot 생성
  p <- ggplot() +
    scale_y_continuous(breaks = 1:para$J, limits = c(0.5, para$J+1)) +
    
    # 구분선 추가
    {if (!is.null(result2))
      geom_vline(xintercept = max(result$enter_times), 
                 linetype = "dashed", color = "gray50", linewidth = 0.3)} +
    
    # Result1 데이터 플로팅
    {if(nrow(main_plot_data1) > 0)
      geom_point(data = main_plot_data1, 
                 aes(x = time, y = dose), 
                 shape = 1, size = 3)} +
    {if(nrow(backfill_plot_data1) > 0)
      geom_point(data = backfill_plot_data1, 
                 aes(x = time, y = dose), 
                 shape = 2, size = 3)} +
    
    # Result2 데이터 플로팅
    {if(!is.null(result2) && nrow(main_plot_data2) > 0)
      geom_point(data = main_plot_data2, 
                 aes(x = time, y = dose), 
                 shape = 1, size = 3)} +
    {if(!is.null(result2) && nrow(backfill_plot_data2) > 0)
      geom_point(data = backfill_plot_data2, 
                 aes(x = time, y = dose), 
                 shape = 2, size = 3)} +
    
    # Subject 번호 표시
    {if(nrow(main_plot_data1) > 0)
      geom_text(data = main_plot_data1,
                aes(x = time, y = dose, label = subject_index),
                vjust = 2, size = 3)} +
    {if(nrow(backfill_plot_data1) > 0)
      geom_text(data = backfill_plot_data1,
                aes(x = time, y = dose, label = subject_index),
                vjust = 2, size = 3)} +
    
    {if(!is.null(result2) && nrow(main_plot_data2) > 0)
      geom_text(data = main_plot_data2,
                aes(x = time, y = dose, label = subject_index),
                vjust = 2, size = 3)} +
    {if(!is.null(result2) && nrow(backfill_plot_data2) > 0)
      geom_text(data = backfill_plot_data2,
                aes(x = time, y = dose, label = subject_index),
                vjust = 2, size = 3)} +
    
    # Result1 화살표와 라벨
    {if(nrow(all_responses1) > 0) 
      geom_segment(data = all_responses1,
                   aes(x = adjusted_enter_time, xend = adjusted_enter_time,
                       y = dose, yend = dose + y_offset,
                       linetype = type, color = type),
                   linewidth = 0.2)} +
    {if(nrow(all_responses1) > 0) 
      geom_segment(data = all_responses1,
                   aes(x = adjusted_enter_time, xend = response_time,
                       y = dose + y_offset, yend = dose + y_offset,
                       linetype = type, color = type),
                   arrow = arrow(length = unit(0.1, "cm"), type = "open", angle = 15),
                   linewidth = 0.2)} +
    
    # Result2 화살표와 라벨
    {if(!is.null(result2) && nrow(all_responses2) > 0) 
      geom_segment(data = all_responses2,
                   aes(x = adjusted_enter_time, xend = adjusted_enter_time,
                       y = dose, yend = dose + y_offset,
                       linetype = type, color = type),
                   linewidth = 0.2)} +
    {if(!is.null(result2) && nrow(all_responses2) > 0) 
      geom_segment(data = all_responses2,
                   aes(x = adjusted_enter_time, xend = response_time,
                       y = dose + y_offset, yend = dose + y_offset,
                       linetype = type, color = type),
                   arrow = arrow(length = unit(0.1, "cm"), type = "open", angle = 15),
                   linewidth = 0.2)} +
    
    # Response labels
    {if(nrow(all_responses1) > 0)
      geom_text(data = all_responses1,
                aes(x = response_time, y = dose + y_offset, label = type,
                    color = type),
                hjust = -0.5)} +
    {if(!is.null(result2) && nrow(all_responses2) > 0)
      geom_text(data = all_responses2,
                aes(x = response_time, y = dose + y_offset, label = type,
                    color = type),
                hjust = -0.5)} +
    
    scale_linetype_manual(values = c("e" = "solid", "t" = "dashed")) +
    scale_color_manual(values = c("e" = "blue", "t" = "red")) +
    
    labs(x = "Time(Months)", y = "Dose Level", 
         title = "Clinical Trial Timeline") +
    
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  return(p)
}



single_plot_curve <- function(tox_sample, eff_sample, scenario) {
  # Toxicity data
  tox_means <- colMeans(as.matrix(tox_sample))
  tox_quantiles <- apply(tox_sample, 2, function(x) quantile(x, c(0.05, 0.95)))
  
  # Efficacy data
  eff_means <- colMeans(as.matrix(eff_sample))
  eff_quantiles <- apply(eff_sample, 2, function(x) quantile(x, c(0.05, 0.95)))
  
  plot_data <- data.frame(
    dose = rep(para$doses, 2),
    Mean = c(tox_means, eff_means),
    Lower = c(tox_quantiles[1,], eff_quantiles[1,]),
    Upper = c(tox_quantiles[2,], eff_quantiles[2,]),
    Type = rep(c("Toxicity", "Efficacy"), each = length(para$doses)),
    True = c(scenario$pi_T, scenario$pi_E)  # 실제 값 추가
  )
  
  p <- ggplot(plot_data, aes(x = dose, color = Type, fill = Type)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2) +
    geom_line(aes(y = Mean), size = 1) +
    geom_point(aes(y = True), size = 3) +  # 실제 값 점으로 표시
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("blue", "red")) +
    scale_fill_manual(values = c("blue", "red")) +
    labs(title = "Toxicity and Efficacy Probability by Dose",
         x = "Dose",
         y = "Probability") +
    theme_minimal()
  
  return(p)
}


single_plot_curve2 <- function(tox_sample, eff_sample, scenario, model_name) {
  # NA 무시하고 평균/표준편차 계산
  tox_means <- colMeans(tox_sample, na.rm = TRUE)
  tox_sd <- apply(tox_sample, 2, sd, na.rm = TRUE)
  eff_means <- colMeans(eff_sample, na.rm = TRUE)
  eff_sd <- apply(eff_sample, 2, sd, na.rm = TRUE)
  
  # 모든 값이 NA인 열 찾아내기
  tox_all_na <- apply(tox_sample, 2, function(x) all(is.na(x)))
  eff_all_na <- apply(eff_sample, 2, function(x) all(is.na(x)))
  
  # 유효한 값만 필터링
  valid_tox_idx <- which(!tox_all_na)
  valid_eff_idx <- which(!eff_all_na)
  
  # 하한/상한 계산
  tox_lower <- pmax(tox_means - tox_sd, 0)
  tox_upper <- pmin(tox_means + tox_sd, 1)
  eff_lower <- pmax(eff_means - eff_sd, 0)
  eff_upper <- pmin(eff_means + eff_sd, 1)
  
  # plot용 데이터 생성
  df_tox <- data.frame(
    dose  = para$doses[valid_tox_idx],
    Mean  = tox_means[valid_tox_idx],
    Lower = tox_lower[valid_tox_idx],
    Upper = tox_upper[valid_tox_idx],
    Type  = "Toxicity",
    True  = scenario$pi_T[valid_tox_idx],
    stringsAsFactors = FALSE
  )
  
  df_eff <- data.frame(
    dose  = para$doses[valid_eff_idx],
    Mean  = eff_means[valid_eff_idx],
    Lower = eff_lower[valid_eff_idx],
    Upper = eff_upper[valid_eff_idx],
    Type  = "Efficacy",
    True  = scenario$pi_E[valid_eff_idx],
    stringsAsFactors = FALSE
  )
  
  plot_data <- rbind(df_tox, df_eff)
  # 그래프
  p <- ggplot(plot_data, aes(x = dose, color = Type, fill = Type)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2) +
    geom_line(aes(y = Mean), size = 1) +
    geom_point(aes(y = Mean), size = 3) +
    geom_point(aes(y = True), shape = 4, size = 3, stroke = 1.5) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("Efficacy" = "green4", "Toxicity" = "red3")) +
    scale_fill_manual(values = c("Efficacy" = "green4", "Toxicity" = "red3")) +
    labs(
      title = paste("Estimated Toxicity and Efficacy Probability by Dose:", model_name),
      subtitle = "Points = estimated means, 'X' = true values, bands = ±1 SD",
      x = "Dose",
      y = "Probability"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80")
    )
  
  return(p)
}


save_plot <- function(plot, scenario_name, stage2_model) {
  # Create directory for plots if it doesn't exist
  plot_dir <- "simulation_plots"
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  # Extract stage2 model name or use "None" if NULL
  if (identical(deparse(stage2_model), deparse(UBLRM.stage2)))  {
    stage2_name <- "UBLRM"
  } else if (identical(deparse(stage2_model), deparse(UBOIN.stage2))) {
    stage2_name <- "UBOIN"
  } else {
    stage2_name <- "None"
  }
  
  # Create filename
  filename <- file.path(plot_dir, paste0(scenario_name, "_", stage2_name, ".png"))
  
  # Save the plot
  ggsave(filename, plot, width = 10, height = 6, dpi = 300)
  message("Plot saved to: ", filename)
  
  return(filename)
}