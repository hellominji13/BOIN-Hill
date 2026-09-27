generate_random_scenario <- function(para) {
  pi_T <- sort(runif(para$J, min = 0, max = 1))
  pi_E <- c(sort(runif(3, min = 0, max = 1)),sort(runif(para$J-3, min = 0, max = 1),decreasing = TRUE))
  return(list(pi_T = pi_T, pi_E = pi_E))
}

plot_scenario <- function(para, scenario, sim, name = NULL, utility = NULL, show_legend = FALSE) {  
  par(mar = c(3, 2, 3, 1)) 
  
  # Round utility scores to two decimal places and create labels for doses with utility
  if (!is.null(utility)) {
    dose_labels <- paste(para$doses, "(", round(utility, 0), ")", sep = "")
  }
  else {
    dose_labels <- para$doses
  }
  # Plot 1: P_Tj(1) and P_Ej(1) vs Dose j
  plot(para$doses, scenario$pi_T, type = "l", col = "blue", ylim = c(0, 1), 
       xlab = "Dose (utility score)", ylab = "Probability", axes = FALSE)
  axis(1, at = para$doses, labels = FALSE)  # Remove las parameter
  axis(2, ylim = c(0, 1), col = "black")
  
  # Rotate x-axis labels by 45 degrees
  text(x = para$doses, y = par("usr")[3] - 0.02, labels = dose_labels, srt = 45, 
       adj = 1, xpd = TRUE)
  
  lines(para$doses, scenario$pi_E, type = "l", col = "red")
  
  # Add DLT threshold as a dashed line
  abline(h = para$DLT, col = "black", lty = 2)  # Adds dashed horizontal line at DLT threshold
  
  # Add DLT label near the dashed line
  text(x = para$doses[2], y = para$DLT, labels = "DLT", pos = 2, col = "black", offset = 0.5) # Move label away from axis
  
  # Find OBD
  OBD <- OBD_find(para, scenario)
  
  # Add OBD label on the y-axis
  if (!all(is.na(OBD))){
    text(x = para$doses[OBD], y = rep(0.05, length(OBD)), labels = "OBD", 
       pos = 2, col = "green", offset = 0.5)  # Adjust y position as needed
  # Add marker (circle) at OBD position on x-axis
    points(para$doses[OBD], rep(0, length(OBD)), pch = 16, col = "green", cex = 1.5)
  }
  # Add name label if provided
  if (!is.null(name)) {
    mtext(text = name, side = 3, line = 0.5, at = mean(para$doses), col = "black", cex = 1)
  }
  
  # Add legend
  if (show_legend) {
    legend("topright", legend = c("P_Tj(1)", "P_Ej(1)"), col = c("blue", "red"), lty = 1)
  }
}

# 
# # OR pT PE 시 각각 확률 retrun
# calculate_joint_probabilities <- function(scenario) {
#   pT <- scenario$pi_T  # Toxicity
#   pE <- scenario$pi_E  # Efficacy
# 
# 
#   # p_te (toxicity, efficacy) if 1 : occur, 0 : no response
#   p11 <- numeric(length(pT))
#   p10 <- numeric(length(pT)) # Toxicity but not efficacy
#   p01 <- numeric(length(pT))
#   p00 <- numeric(length(pT))
# 
#   if (scenario$model=="OR") {
#     OR <- scenario$coff
#     for (i in 1:length(pT)) {
#       if (OR[i] == 1) {
#         p11[i] <- pT[i]*pE[i]
#       }
#       else {
#         a <- 1 + (pT[i] + pE[i]) * (OR[i] - 1)
#         b <- -4 * OR[i] * (OR[i] - 1) * pT[i] * pE[i]
# 
#         # Calculate joint probability p11
#         p11[i] <- (1 / (2 * (OR[i] - 1))) * (a - sqrt(a^2 + b))
#       }
# 
#       # Calculate other joint probabilities
#       p10[i] <- pT[i] - p11[i]
#       p01[i] <- pE[i] - p11[i]
#       p00[i] <- 1 - pT[i] - pE[i] + p11[i]
# 
# 
#       if (p11[i] < 0 || p10[i] < 0 || p01[i] < 0 || p00[i] < 0) {
#         stop(paste("Error: Negative joint probability found at index", i,
#                    "p00:", p00[i], "p01:", p01[i], "p10:", p10[i], "p11:", p11[i]))
#       }
#     }
#   }
# 
#   else if (scenario$model == "Gumbel") {
#     coff <- scenario$coff
#     gumbel <- function(t, e, pT, pE,coff) {
#       (pT)^t * (1 - pT)^(1 - t) * (pE)^e * (1 - pE)^(1 - e) +
#         (pE) * (1 - pE) * (pT) * (1 - pT) * (-1)^(t + e) * ((e^coff - 1) / (e^coff + 1))
#     }
#     for (i in 1:length(pT)) {
#       p00[i] <- gumbel(0, 0, pT[i], pE[i],coff[i])
#       p01[i] <- gumbel(0, 1, pT[i], pE[i],coff[i])
#       p10[i] <- gumbel(1, 0, pT[i], pE[i],coff[i])
#       p11[i] <- gumbel(1, 1, pT[i], pE[i],coff[i])
#     }
#   }
# 
#   else {
#     stop("Error: model can be OR or Gumbel")
#   }
# 
#   joint_probabilities <- data.frame(
#     p00 = p00,
#     p01 = p01,
#     p10 = p10,
#     p11 = p11
#   )
# 
#   return(joint_probabilities)
# }


calculate_joint_probabilities <- function(scenario, rho = scenario$coff) {
  # Weibull 파라미터 계산 함수
  calculate_weibull_params <- function(p) {
    # window=1에서 p, window=0.5에서 p/2가 되도록 파라미터 계산
    k <- log(log(1-p/2)/log(1-p))/log(0.5)
    lambda <- 1/((-log(1-p))^(1/k))
    return(list(shape = k, scale = lambda))
  }

  pT <- scenario$pi_T  # Toxicity probabilities
  pE <- scenario$pi_E  # Efficacy probabilities

  # 결과 저장을 위한 벡터 초기화
  p00 <- numeric(length(pT))
  p01 <- numeric(length(pT))
  p10 <- numeric(length(pT))
  p11 <- numeric(length(pT))

  # 각 용량 수준에 대해 계산
  for (i in 1:length(pT)) {
    # Weibull 파라미터 계산
    weibull_T <- calculate_weibull_params(pT[i])
    weibull_E <- calculate_weibull_params(pE[i])

    # Copula 객체 생성
    mv.cop <- copula::mvdc(
      copula = copula::frankCopula(param = rho, dim = 2),
      margins = c("weibull", "weibull"),
      paramMargins = list(
        list(shape = weibull_T$shape, scale = weibull_T$scale),
        list(shape = weibull_E$shape, scale = weibull_E$scale)
      )
    )

    # 큰 수의 샘플을 생성하여 확률 추정
    n_samples <- 100000
    samples <- copula::rMvdc(n_samples, mv.cop)

    # window=1 이내 발생 여부로 확률 계산
    T_occur <- samples[,1] <= 1
    E_occur <- samples[,2] <= 1

    # 결합 확률 계산
    p11[i] <- mean(T_occur & E_occur)    # 둘 다 발생
    p10[i] <- mean(T_occur & !E_occur)   # Toxicity만 발생
    p01[i] <- mean(!T_occur & E_occur)   # Efficacy만 발생
    p00[i] <- mean(!T_occur & !E_occur)  # 둘 다 미발생
  }

  # 결과를 데이터프레임으로 반환
  joint_probabilities <- data.frame(
    p00 = p00,
    p01 = p01,
    p10 = p10,
    p11 = p11
  )

  # 수치적 오차 검증
  total_prob <- rowSums(joint_probabilities)
  if (any(abs(total_prob - 1) > 0.01)) {
    warning("Warning: Some joint probabilities do not sum to 1 (within 0.01 tolerance)")
  }

  return(joint_probabilities)
}

# OBD 를 정의 하는 함수 
OBD_find <- function(para, scenario){
  
  # score table plot 
  score_matrix <- matrix(scenario$score, nrow = 2, byrow = TRUE, dimnames = list("Toxicity" = c("0", "1(occur)"),"Efficacy" = c("0", "1(occur)")))

  # Find OBD
  prob <- calculate_joint_probabilities(scenario)
  score_porb <- t(apply(prob, 1, function(row) row * scenario$score))
  row_sums <- rowSums(score_porb)
  
  dosen <- 1:para$J
  
  admissible_doses <- dosen[scenario$pi_T[dosen] < para$T_admissible]
  admissible_doses <- admissible_doses[scenario$pi_E[admissible_doses] > para$E_admissible]

  if (length(admissible_doses) > 0) {
    max_row_sum <- max(row_sums[admissible_doses])
    admissible_doses <- admissible_doses[row_sums[admissible_doses] >= (max_row_sum - 5)]
    OBD = admissible_doses
  }
  else {
    OBD <- NA
  }
  return(OBD)
}