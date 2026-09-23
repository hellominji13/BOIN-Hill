# Backfill을 하지 않을 때의 주어진 prob에 대한 시나리오 생성 

# generate_toxicity_efficacy <- function(scenario, subject_num) {
#   J<-para$J
#   #output : Jx4 table (outcome(p_te) p00, p01,p10,p11 )
#   joint_probs <- calculate_joint_probabilities(scenario) 
#   response_matrix <- matrix(0, nrow = J * 2, ncol = subject_num)
#   for (j in 1:J) {
#     p00 <- joint_probs$p00[j]
#     p01 <- joint_probs$p01[j]
#     p10 <- joint_probs$p10[j]
#     p11 <- joint_probs$p11[j]
#     
#     for (s in 1:subject_num) {
#       response <- sample(c("00", "01", "10", "11"), 1, prob = c(p00, p01, p10, p11))
#       
#       if (response == "00") {
#         response_matrix[j, s] <- 0  # No toxicity
#         response_matrix[j + J, s] <- 0  # No efficacy
#       } 
#       else if (response == "01") {
#         response_matrix[j, s] <- 0  # No toxicity
#         response_matrix[j + J, s] <- 1  # Efficacy
#       } 
#       else if (response == "10") {
#         response_matrix[j, s] <- 1  # Toxicity
#         response_matrix[j + J, s] <- 0  # No efficacy
#       } 
#       else if (response == "11") {
#         response_matrix[j, s] <- 1  # Toxicity
#         response_matrix[j + J, s] <- 1  # Efficacy
#       }
#       
#     }
#   }
#   outcome_table <- as.data.frame(response_matrix)
#   colnames(outcome_table) <- c(paste("subject", 1:subject_num, sep="_"))
#   row_names <- c(paste("Toxicity_Dose", 1:J, sep="_"), paste("Efficacy_Dose", 1:J, sep="_"))
#   rownames(outcome_table) <- row_names
#   return(outcome_table)
# }




# Time 분포를 균등분포(Uniform) 또는 와이블분포(Weibull)로 고려한 시나리오 생성 
generate_toxicity_efficacy_time <- function(scenario, subject_num, distribution = "uniform") {
  mv.ncop <- list()
  rho <- 1
  
  # Uniform 분포 매개변수 계산 함수
  calculate_uniform_params <- function(p) {
    # 발생 확률이 p가 되도록 최대값 계산
    # Uniform(0, max)에서 window 내에 발생할 확률이 p가 되기 위한 max값
    max_time <- para$window / p
    return(list(min = 0, max = max_time))
  }
  
  # Weibull 파라미터 계산 함수
  calculate_weibull_params <- function(p) {
    # 1개월(30일)에서 p, 0.5개월(15일)에서 p/2가 되도록 파라미터 계산
    k <- log(log(1-p/2)/log(1-p))/log(0.5)
    lambda <- 1/((-log(1-p))^(1/k))
    return(list(shape = k, scale = lambda))
  }
  
  # 분포 선택 확인
  if(!(distribution %in% c("uniform", "weibull"))) {
    stop("분포는 'uniform' 또는 'weibull'이어야 합니다.")
  }
  
  # J는 시나리오에서 정의된 독성과 효능에 대한 수치
  for (i in 1:para$J) {
    psi.T <- scenario$pi_T[i]
    psi.E <- scenario$pi_E[i]
    if(distribution == "uniform") {
      # Uniform 파라미터 계산
      uniform_T <- calculate_uniform_params(psi.T)
      uniform_E <- calculate_uniform_params(psi.E)
      
      # Frank copula로 상관관계 유지
      mv.ncop <- append(mv.ncop, copula::mvdc(
        copula = copula::frankCopula(param = rho, dim = 2),
        margins = c("unif", "unif"),
        paramMargins = list(
          list(min = uniform_T$min, max = uniform_T$max),
          list(min = uniform_E$min, max = uniform_E$max)
        )
      ))
    } else { # weibull
      # Weibull 파라미터 계산
      weibull_T <- calculate_weibull_params(psi.T)
      weibull_E <- calculate_weibull_params(psi.E)
      
      # Frank copula로 상관관계 유지
      mv.ncop <- append(mv.ncop, copula::mvdc(
        copula = copula::frankCopula(param = rho, dim = 2),
        margins = c("weibull", "weibull"),
        paramMargins = list(
          list(shape = weibull_T$shape, scale = weibull_T$scale),
          list(shape = weibull_E$shape, scale = weibull_E$scale)
        )
      ))
    }
  }
  
  # 결과 저장을 위한 매트릭스 초기화
  response_matrix <- matrix(0, nrow = para$J * 2, ncol = subject_num)
  time_matrix <- matrix(0, nrow = para$J * 2, ncol = subject_num)
  # 각 대상에 대해 독성과 효능 시간 및 반응 생성
  for (j in 1:para$J) {
    for (s in 1:subject_num) {
      
      time.te <- copula::rMvdc(1, mv.ncop[[j]])  # 1개의 샘플 생성
      # 시간이 window를 초과하면 window 값으로 자름
      time_matrix[j, s] <- min(time.te[1], para$window)  # 독성 시간
      time_matrix[j + para$J, s] <- min(time.te[2], para$window)  # 효능 시간
      
      # window 이내에 발생 여부 확인
      response_matrix[j, s] <- as.numeric(time.te[1] <= para$window)
      response_matrix[j + para$J, s] <- as.numeric(time.te[2] <= para$window)
    }
  }
  
  # 결과를 데이터프레임으로 변환
  outcome_table <- as.data.frame(response_matrix)
  colnames(outcome_table) <- paste("subject", 1:subject_num, sep = "_")
  row_names <- c(
    paste("Toxicity_Dose", 1:para$J, sep = "_"), 
    paste("Efficacy_Dose", 1:para$J, sep = "_"))
  rownames(outcome_table) <- row_names
  
  time_table <- as.data.frame(time_matrix)
  colnames(time_table) <- paste("subject", 1:subject_num, sep = "_")
  rownames(time_table) <- row_names
  
  # 파라미터 저장
  params <- list()
  if(distribution == "uniform") {
    params$toxicity <- lapply(1:para$J, function(i) calculate_uniform_params(scenario$pi_T[i]))
    params$efficacy <- lapply(1:para$J, function(i) calculate_uniform_params(scenario$pi_E[i]))
  } else { # weibull
    params$toxicity <- lapply(1:para$J, function(i) calculate_weibull_params(scenario$pi_T[i]))
    params$efficacy <- lapply(1:para$J, function(i) calculate_weibull_params(scenario$pi_E[i]))
  }
  
  return(list(outcome_table = outcome_table, 
              time_table = time_table,
              distribution = distribution,
              params = params))
}
