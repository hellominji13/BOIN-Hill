

beta_quantiles <- function(q_p,p,c) { #Pr( < q_p) = p
  if (q_p >= p) {
    quantiles <- qbeta(c, (log(p) / log(q_p)), 1)
  } else {
    quantiles <- qbeta(c, 1, (log(1 - p) / log(1 - q_p)))}
  return(quantiles)
}
  
MVN_quantiles<- function(prob_list) {
  quantile_matrices <- lapply(prob_list, function(probs) {
    quant_values <- quantile(probs, probs = para$c)
    # 분위수 값들을 행렬로 변환
    unlist(quant_values)  # 벡터 형태로 변환
  })
  # 리스트를 행렬로 변환
  do.call(rbind, quantile_matrices)
}

#Assume 
#Pr(min_dose > T_max(or E_max)) = 0.025, Pr(max_dose < T_min) = 0.975
Make_Beta_q_j <- function(para,min_dose,max_dose,kind) {
  Beta_q_j <- matrix(0, nrow = para$J, ncol = length(para$c))
  x_dose <- seq(0.5, 500, length.out = (500 %/% 0.5))
  mean_vale <- rep(0, length(x_dose))
  
  # plot 을 위해서 필요 
  whole_q_j <- matrix(0, nrow = length(x_dose), ncol = length(para$c))
  max_index <- para$doses[max_dose] %/% 0.5
  min_index <- para$doses[min_dose] %/% 0.5
  
  # 각 dose 의 mean 값 계산하기 
  # E 타입일 경우
  if (kind == "E") {
    mean_vale[min_index] <- beta_quantiles(para$E_max, 0.95, para$c)[2] 
    mean_vale[max_index] <- beta_quantiles(para$E_min, 0.05, para$c)[2] 

    # max_dose 이하에 대해서는 선형 (로지스틱 선형 형태)
    for (i in 1:length(x_dose)) {
      if (x_dose[i] <= para$doses[max_dose]) {
        prob <- ((logit(mean_vale[min_index])-logit(mean_vale[max_index]))/(log(para$doses[min_dose] / para$doses[max_dose]))) * (log(x_dose[i]/para$doses[max_dose])) + logit(mean_vale[max_index])
        mean_vale[i] <- inv_logit(prob)
        #mean_vale[i] <- para$E_admissible
      } else {
        # max_dose 이후에는 감소 -->유지 
        prob <- ((logit(mean_vale[max_index]) - logit(mean_vale[min_index])) / (log(para$doses[max_dose] / para$doses[para$J]))) * (log(x_dose[i] /para$doses[max_dose])) + logit(mean_vale[max_index])
        mean_vale[i] <- inv_logit(prob)
        # mean_vale[i] <- mean_vale[max_index]
      }
    }
  }
  #logistic linear 하다가고 가정
  if (kind=="T"){
    mean_vale[min_index]  <- beta_quantiles(para$T_max, 0.95, para$c)[2] 
    mean_vale[max_index] <- beta_quantiles(para$T_min, 0.05, para$c)[2] 
    prob <- ((logit(mean_vale[min_index])-logit(mean_vale[max_index]))/(log(para$doses[min_dose] / para$doses[max_dose]))) * (log(x_dose/para$doses[max_dose])) + logit(mean_vale[max_index])
    mean_vale <- inv_logit(prob)
  }
  
  #beta distribution 구한 후 quarantines 구하기 
  for (index in 1:para$J) {
    Beta_q_j[index, ] <- beta_quantiles(mean_vale[(para$doses[index] %/% 0.5)], 0.5,para$c)}
  for (index in 1:length(x_dose)) {
    whole_q_j[index, ] <- beta_quantiles(mean_vale[index], 0.5, para$c)
  }
  return (list(Beta_q_j=Beta_q_j,whole_q_j = whole_q_j))
}

Make_MVN_q_j <- function(para, mean_theta, sd_theta, corr, kind) {
  prob_list <- vector("list", length(para$doses)) #여기에 Toxicity 와 Efficacy 의 확률값을 담아야함 
  
  #Parameter의 사전 분포를 넣으면, 확률값의 분위값이 나오도록 계산하는 코드 
  if (kind == "T") {
    cov_matrix <- matrix(c(sd_theta[1]^2, corr * sd_theta[1] * sd_theta[2],
                           corr * sd_theta[1] * sd_theta[2], sd_theta[2]^2), ncol=2, byrow=T)
    samples <- mvrnorm(n = 500, mu = mean_theta, Sigma = cov_matrix)
    
    for (k in seq_along(para$doses)) {
      dose_probs <- numeric(500)  
      for (j in 1:500) {
        theta <- samples[j, ]  
        dose_probs[j] <- BLRM_model_T(para, k, theta)
      }
      prob_list[[k]] <- dose_probs
    }
  }
  if (kind=="E") {
    # 대각행렬: 비대각 요소는 모두 0
    cov_matrix <- matrix(0, ncol=length(mean_theta), nrow=length(mean_theta))
    diag(cov_matrix) <- sd_theta^2  # 대각 요소만 분산값 설정
    samples <- mvrnorm(n = 500, mu = mean_theta, Sigma = cov_matrix)
    for (k in seq_along(para$doses)) {
      dose_probs <- numeric(500)  
      for (j in 1:500) {
        theta <- samples[j, ]  
        dose_probs[j] <- BLRM_model_E(para, k, theta)
      }
      prob_list[[k]] <- dose_probs
    }
  }
  MVN_q_j <- MVN_quantiles(prob_list)
  return(list(MVN_q_j = MVN_q_j))
}


Make_prior <- function(para, kind, plot = FALSE, max_iterations = 3000, tolerance = 0.1) {
  # 초기 설정
  if (kind == "T") {
    current_mean_theta <- c(-1, 1)
    current_sd_theta <- c(0.5, 0.5)
    current_cor <- 0
    Beta_q_j <- Make_Beta_q_j(para, 1, para$J, kind)
  } else if (kind == "E") {
    current_mean_theta <- c(2,2,2,2)
    current_sd_theta <- c(1, 1, 1,1)
    current_cor <- 0
    Beta_q_j <- Make_Beta_q_j(para, 1, para$J %/% 2 +1, kind)
  } else {
    stop("Invalid 'kind' parameter. Must be 'T' or 'E'.")
  }
  
  # 초기화
  best_mean_theta <- current_mean_theta
  best_sd_theta <- current_sd_theta
  best_cor <- current_cor
  
  # 초기 차이 계산
  initial_MVN <- Make_MVN_q_j(para, current_mean_theta, current_sd_theta, current_cor, kind)
  best_diff <- max(abs(initial_MVN$MVN_q_j - Beta_q_j$Beta_q_j))  # MSE 사용
  if (plot) {
    print("MINJI")
    if (kind == "T") {
      print(plot_curves(para, initial_MVN$MVN_q_j, NULL, Beta_q_j, NULL)$toxicity)
    } else {
      print(plot_curves(para, NULL, initial_MVN$MVN_q_j, NULL, Beta_q_j)$efficacy)
    }
  }
  # 반복 최적화
  for (iter in 1:max_iterations) {
    # 스텝 크기 조절
    step_size <- max(0.1, 1 - iter/max_iterations)
    
    # 여러 후보 생성 및 평가
    n_candidates <- 5
    best_candidate_diff <- best_diff
    
    for(candidate in 1:n_candidates) {
      # 파라미터 변화량 생성
      mean_changes <- rnorm(length(current_mean_theta), 0, 1 * step_size)
      sd_changes <- rnorm(length(current_sd_theta), 0, 1 * step_size)
      cor_change <- rnorm(1, 0, 0.1 * step_size)
      
      # 새로운 파라미터 제안
      proposed_mean_theta <- best_mean_theta + mean_changes
      proposed_sd_theta <- pmax(0.1, best_sd_theta + sd_changes)  # 최소값 제한
      proposed_cor <- pmax(-0.99, pmin(0.99, best_cor + cor_change))  # 상관계수 범위 제한
      
      # 제안된 분포 계산
      MVN_q_j_proposed <- Make_MVN_q_j(para, proposed_mean_theta, proposed_sd_theta, proposed_cor, kind)
      
      # MSE 계산
      proposed_diff <- max(abs(MVN_q_j_proposed$MVN_q_j - Beta_q_j$Beta_q_j))
      
      # 더 나은 결과를 찾은 경우 업데이트
      if (proposed_diff < best_candidate_diff) {
        best_candidate_diff <- proposed_diff
        best_candidate_mean <- proposed_mean_theta
        best_candidate_sd <- proposed_sd_theta
        best_candidate_cor <- proposed_cor
      }
    }
    
    # 개선된 결과가 있으면 업데이트
    if (best_candidate_diff < best_diff) {
      best_diff <- best_candidate_diff
      best_mean_theta <- best_candidate_mean
      best_sd_theta <- best_candidate_sd
      best_cor <- best_candidate_cor
      
      # 플롯 업데이트
      if (plot) {
        MVN_q_j <- Make_MVN_q_j(para, best_mean_theta, best_sd_theta, best_cor, kind)
        if (kind == "T") {
          print(plot_curves(para, MVN_q_j$MVN_q_j, NULL, Beta_q_j, NULL)$toxicity)
        } else {
          print(plot_curves(para, NULL, MVN_q_j$MVN_q_j, NULL, Beta_q_j)$efficacy)
        }
      }
    }
    
    # 수렴 체크
    if (best_diff < tolerance) {
      break
    }
  }
  
  # 결과 반환
  return(list(
    mean = best_mean_theta,
    se = best_sd_theta,
    corr = best_cor,
    final_diff = best_diff,
    iterations = iter
  ))
}



plot_curves <- function(para, MVN_q_j_T, MVN_q_j_E, Beta_q_j_T, Beta_q_j_E) {
  # Toxicity plot 데이터 준비
  if(!is.null(MVN_q_j_T) && !is.null(Beta_q_j_T)) {
    plot_data_T_beta <- data.frame(
      dose = para$doses,
      Mean = Beta_q_j_T$Beta_q_j[,2],
      Lower = Beta_q_j_T$Beta_q_j[,1],
      Upper = Beta_q_j_T$Beta_q_j[,3],
      Method = "Beta Prior"
    )
    
    plot_data_T_MVN <- data.frame(
      dose = para$doses,
      Mean = MVN_q_j_T[,2],
      Lower = MVN_q_j_T[,1],
      Upper = MVN_q_j_T[,3],
      Method = "MVN Fitting"
    )
    
    plot_data_T <- rbind(plot_data_T_beta, plot_data_T_MVN)
    
    plot_Toxicity <- ggplot(plot_data_T, aes(x = dose)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = 0.2) +
      geom_line(aes(y = Mean, color = Method), size = 1) +
      scale_fill_manual(values = c("Beta Prior" = "black", "MVN Fitting" = "red")) +
      scale_color_manual(values = c("Beta Prior" = "black", "MVN Fitting" = "red")) +
      labs(title = "Toxicity Probability by Dose",
           x = "Dose",
           y = "Probability",
           fill = "Method",
           color = "Method") +
      theme_minimal()
  } else {
    plot_Toxicity <- NULL
  }
  
  # Efficacy plot 데이터 준비
  if(!is.null(MVN_q_j_E) && !is.null(Beta_q_j_E)) {
    plot_data_E_beta <- data.frame(
      dose = para$doses,
      Mean = Beta_q_j_E$Beta_q_j[,2],
      Lower = Beta_q_j_E$Beta_q_j[,1],
      Upper = Beta_q_j_E$Beta_q_j[,3],
      Method = "Beta Prior"
    )
    
    plot_data_E_MVN <- data.frame(
      dose = para$doses, 
      Mean = MVN_q_j_E[,2],
      Lower = MVN_q_j_E[,1],
      Upper = MVN_q_j_E[,3],
      Method = "MVN Fitting"
    )
    
    plot_data_E <- rbind(plot_data_E_beta, plot_data_E_MVN)
    
    plot_Efficacy <- ggplot(plot_data_E, aes(x = dose)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = 0.2) +
      geom_line(aes(y = Mean, color = Method), size = 1) +
      scale_fill_manual(values = c("Beta Prior" = "black", "MVN Fitting" = "blue")) +
      scale_color_manual(values = c("Beta Prior" = "black", "MVN Fitting" = "blue")) +
      labs(title = "Efficacy Probability by Dose",
           x = "Dose",
           y = "Probability",
           fill = "Method",
           color = "Method") +
      theme_minimal()
  } else {
    plot_Efficacy <- NULL
  }
  
  return(list(toxicity = plot_Toxicity, efficacy = plot_Efficacy))
}