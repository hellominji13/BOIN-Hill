
logit <- function(p) log(p / (1 - p))

inv_logit <- function(x) (1 / (1 + exp(-x)))




Bayesian.Hill <- function(Y, x,
                          burn.in = 3000, B = 1000, seed = 123,
                          delta0_initial = 0.2, delta1_initial = 0.6,
                          E_prior,
                          unif = TRUE, chol = TRUE,
                          mat.calc = FALSE) {
  set.seed(seed)
  from.prior.only <- is.null(Y) || length(Y) == 0
  I <- B + burn.in
  N <- if (from.prior.only) 0 else length(Y)
  
  delta0 <- delta1 <- b0 <- b1 <- numeric()
  
  if (from.prior.only) {
    sim <- 1000000
    delta0 <- rbeta(sim, E_prior$delta0_a, E_prior$delta0_b)
    delta1 <- rbeta(sim, E_prior$delta1_a, E_prior$delta1_b)
    b0     <- rnorm(sim, E_prior$b0_mean, E_prior$b0_sd)
    b1     <- rnorm(sim, E_prior$b1_mean, E_prior$b1_sd)
    return(data.frame(delta0, delta1, b0, b1))
  }
  
  # Prior from E_prior list
  a00 <- E_prior$delta0_a; a01 <- E_prior$delta0_b
  a10 <- E_prior$delta1_a; a11 <- E_prior$delta1_b
  mu0 <- E_prior$b0_mean;  s2.0 <- (E_prior$b0_sd)^2
  mu1 <- E_prior$b1_mean;  s2.1 <- (E_prior$b1_sd)^2
  
  # 초기값과 저장 벡터
  delta0 <- delta1 <- b0 <- b1 <- numeric()
  delta0_current <- delta0_initial
  delta1_current <- delta1_initial
  b0_current <-  E_prior$b0_mean
  b1_current <-  E_prior$b1_mean
  
  # 유틸리티 벡터
  if (unif) U <- matrix(runif(N * I), N, I)
  if (chol) G <- matrix(rnorm(2 * I), 2, I)
  if (mat.calc) X <- cbind(1, x)
  sY <- sum(Y); sx <- sum(x)
  
  for (i in 1:I) {
    # Step 1: sample z | rest
    BX <- b0_current + b1_current * x
    temp <- exp(BX)
    temp2 <- ifelse(Y == 1,
                    delta0_current / delta1_current,
                    (1 - delta0_current) / (1 - delta1_current))
    eta <- temp / (temp + temp2)
    z_next <- if (unif) (U[, i] < eta) else rbinom(n = N, size = 1, prob = eta)
    
    # Step 2: sample delta0, delta1 | rest
    sYz <- sum(Y * z_next)
    sz <- sum(z_next)
    delta0_current <- rbeta(1, a00 + sY - sYz, a01 + N - sz - sY + sYz)
    delta1_current <- rbeta(1, a10 + sYz, a11 + sz - sYz)
    
    # Step 3: sample b0, b1 | rest via Polya-Gamma
    W <- BayesLogit::rpg(N, 1, BX)
    if (mat.calc) {
      SIGMA <- solve(t(X) %*% diag(W) %*% X + diag(c(1/s2.0, 1/s2.1)))
      m <- SIGMA %*% (t(X) %*% (z_next - 0.5) + diag(c(1/s2.0, 1/s2.1)) %*% c(mu0, mu1))
      m1 <- m[1]; m2 <- m[2]
    } else {
      sW <- sum(W); sWX <- sum(W * x); sWX2 <- sum(W * x^2)
      a <- sW + 1/s2.0; b <- sWX; d <- sWX2 + 1/s2.1
      D <- a * d - b^2
      SIGMA <- (1/D) * matrix(c(d, -b, -b, a), 2, 2)
      k1 <- sum(z_next) - N/2 + mu0/s2.0
      k2 <- sum(x * (z_next - 0.5)) + mu1/s2.1
      m1 <- (d * k1 - b * k2) / D
      m2 <- (-b * k1 + a * k2) / D
    }
    
    b_next <- if (chol) c(m1, m2) + chol(SIGMA) %*% G[, i]
    else MASS::mvrnorm(1, mu = c(m1, m2), Sigma = SIGMA)
    
    b0_current <- b_next[1]; b1_current <- b_next[2]
    
    delta0 <- c(delta0, delta0_current)
    delta1 <- c(delta1, delta1_current)
    b0 <- c(b0, b0_current)
    b1 <- c(b1, b1_current)
  }
  
  return(data.frame(delta0, delta1, b0, b1)[-(1:burn.in), ])
}

Bayesian.Toxicity <- function(Y, x,
                              burn.in = 2000, B = 1000, seed = 123,
                              T_prior,
                              chol = TRUE) {
  set.seed(seed)
  from.prior.only <- is.null(Y) || length(Y) == 0
  I <- burn.in + B
  N <- if (from.prior.only) 0 else length(Y)
  alpha0 <- alpha1 <- numeric()
  
  if (from.prior.only) {
    sim <- 1000
    alpha0 <- rnorm(sim, T_prior$a0_mean, T_prior$a0_sd)
    alpha1 <- rnorm(sim, T_prior$a1_mean, T_prior$a1_sd)
    return(data.frame(alpha0 = alpha0, alpha1 = alpha1)[-(1:burn.in), ])
  }
  
  # prior parameters
  mu0 <- T_prior$a0_mean
  s2.0 <- (T_prior$a0_sd)^2
  mu1 <- T_prior$a1_mean
  s2.1 <- (T_prior$a1_sd)^2
  
  prec0 <- 1 / s2.0
  prec1 <- 1 / s2.1
  
  # initial values
  a0_current <- mu0
  a1_current <- mu1
  a0 <- a1 <- numeric()
  
  # utility
  G <- if (chol) matrix(rnorm(2 * I), 2, I)
  X <- cbind(1, x)
  V0inv <- diag(c(prec0, prec1))
  mu_vec <- c(mu0, mu1)
  
  for (i in 1:I) {
    # Step 1: sample w | a0, a1
    psi <- a0_current + a1_current * x
    w <- BayesLogit::rpg(N, 1, psi)
    
    # Step 2: sample a0, a1 | w, Y
    SIGMA <- solve(t(X) %*% diag(w) %*% X + diag(c(1/s2.0, 1/s2.1)))
    m <- SIGMA %*% (t(X) %*% (Y - 0.5) + diag(c(1/s2.0, 1/s2.1)) %*% c(mu0, mu1))
    
    a_next <- if (chol) c(m[1], m[2]) + chol(SIGMA) %*% G[, i]
    else MASS::mvrnorm(1, mu = c(m[1], m[2]), Sigma = SIGMA)
    
    a0_current <- a_next[1]; a1_current <- a_next[2]
    
    a0 <- c(a0, a0_current)
    a1 <- c(a1, a1_current)
  }
  
  return(data.frame(alpha0 = a0, alpha1 = a1)[-(1:burn.in), ])
}

UBLRM_model_fit_function <- function(para, data,
                                     T_prior, E_prior, thin = 2,
                                     seed = 123) {
  set.seed(seed)
  
  # log-scaled dose
  D <- log(para$doses / para$ref_dose)
  nd <- length(D)
  
  ## --------- DATA-BASED MCMC MODE ---------
  dose.level <- c()
  y.tox <- c()
  y.eff <- c()
  
  for (j in seq_len(nrow(data))) {
    dose_j <- data$dose[j]
    
    dose.level <- c(
      dose.level,
      rep(dose_j, data$n00[j]),
      rep(dose_j, data$n01[j]),
      rep(dose_j, data$n10[j]),
      rep(dose_j, data$n11[j])
    )
    
    y.tox <- c(
      y.tox,
      rep(0, data$n00[j]),
      rep(0, data$n01[j]),
      rep(1, data$n10[j]),
      rep(1, data$n11[j])
    )
    
    y.eff <- c(
      y.eff,
      rep(0, data$n00[j]),
      rep(1, data$n01[j]),
      rep(0, data$n10[j]),
      rep(1, data$n11[j])
    )
  }
  
  # log-scaled X
  X <- log(para$doses[dose.level] / para$ref_dose)
  
  hill.draws <- Bayesian.Hill(Y = y.eff, x = X,
                              E_prior = E_prior)
  
  eff_par_mean <- colMeans(hill.draws)
  eff_sample <- do.call(cbind, lapply(1:nd, function(j) {
    delta0 <- hill.draws[, "delta0"]
    delta1 <- hill.draws[, "delta1"]
    b0     <- hill.draws[, "b0"]
    b1     <- hill.draws[, "b1"]
    delta0 + (delta1) * plogis(b0 + b1 * D[j])
  }))
  pE_pred <- colMeans(eff_sample)
  
  tox.draws <- Bayesian.Toxicity(Y = y.tox, x = X,
                                 T_prior = T_prior)
  tox_sample <- do.call(cbind, lapply(1:nd, function(j) {
    plogis(tox.draws$alpha0 + tox.draws$alpha1 * D[j])
  }))
  pT_pred <- colMeans(tox_sample)
  
  coff_mean <- 0
  joint_mat <- disjoint(tox = pT_pred, eff = pE_pred, coff = coff_mean)
  rownames(joint_mat) <- c("00", "01", "10", "11")
  
  list(
    joint         = joint_mat,
    eff           = pE_pred,
    tox           = pT_pred,
    eff_sample    = eff_sample,
    tox_sample    = tox_sample,
    coff          = coff_mean,
    efficacy_params = as.list(eff_par_mean),
    tox_params = list(
      alpha0 = mean(tox.draws$alpha0),
      alpha1 = mean(tox.draws$alpha1)
    ),
    source = "posterior"
  )
}




disjoint <- function(tox, eff, coff){
  joint <- matrix(NA, nrow = 4, ncol = length(tox))  # 초기값을 NA로 설정
  
  for (j in 1:length(tox)) {
    for (t in 0:1) {
      for (e in 0:1) {
        if (is.na(tox[j]) || is.na(eff[j])) {
          joint[1 + 2 * t + e, j] <- NA  # 명시적으로 NA 할당
        } else {
          joint[1 + 2 * t + e, j] <- (tox[j]^t * (1 - tox[j])^(1 - t) * eff[j]^e * (1 - eff[j])^(1 - e)) +
            (eff[j] * (1 - eff[j]) * tox[j] * (1 - tox[j]) * (-1)^(t + e) * ((exp(coff) - 1) / (exp(coff) + 1)))
        }
      }
    }
  }
  
  return(joint)
}
Admissible_dose_range_function_BRLM <- function(para, prob, current_dose) {
  Admssible_dose_range <- c()
  
  if(!is.na(current_dose)) {  # current_dose가 NA가 아닐 때만 실행
    for (j in 1:para$J) {
      tox_posterior <- prob$tox_sample[, j]  
      eff_posterior <- prob$eff_sample[, j]
      tox_condition <- mean(tox_posterior > para$T_admissible)
      eff_condition <- mean(eff_posterior < (para$E_admissible))
      if (tox_condition <= para$CT && eff_condition <= para$CE) {
        Admssible_dose_range <- c(Admssible_dose_range, j)
      }
    }
  }
  # admissible_extra <- which(prob$tox <= para$T_admissible)  # 독성이 기준 이내
  # Admssible_dose_range <- intersect(Admssible_dose_range, admissible_extra)
  return(Admssible_dose_range)
}




plot_UBLRM_fit <- function(data, mcmc_results) {
  # 실제 관찰된 비율 계산
  observed_tox <- sapply(1:nrow(data), function(i) {
    (data$n10[i] + data$n11[i]) / data$npt[i]
  })
  
  observed_eff <- sapply(1:nrow(data), function(i) {
    (data$n01[i] + data$n11[i]) / data$npt[i]
  })
  
  if (!is.null(mcmc_results$tox_sample)) {
    tox_ci <- t(apply(mcmc_results$tox_sample, 2, quantile, probs = c(0.05, 0.95)))
  }
  
  if (!is.null(mcmc_results$eff_sample)) {
    eff_ci <- t(apply(mcmc_results$eff_sample, 2, quantile, probs = c(0.05, 0.95)))
  }
  
  # 플롯 데이터 준비
  dose_levels <- 1:length(mcmc_results$tox)
  plot_data <- data.frame(
    dose_level = dose_levels,
    observed_tox = observed_tox,
    predicted_tox = mcmc_results$tox,
    observed_eff = observed_eff,
    predicted_eff = mcmc_results$eff
  )
  
  # Combined plot with both toxicity and efficacy
  p <- ggplot(plot_data) +
    # Toxicity components (Red)
    geom_line(aes(x = dose_level, y = predicted_tox), 
              color = "blue", size = 1) +
    geom_point(aes(x = dose_level, y = observed_tox), 
               color = "blue", size = 4) +
    
    # Efficacy components (Blue)
    geom_line(aes(x = dose_level, y = predicted_eff), 
              color = "red", size = 1) +
    geom_point(aes(x = dose_level, y = observed_eff), 
               color = "red", size = 4) +
    
    # Labels and theme
    labs(title = "Toxicity (Blue) and Efficacy (Red): Observed vs Predicted",
         x = "Dose Level",
         y = "Probability") +
    scale_x_continuous(breaks = dose_levels) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  
  # 리본 플롯 추가 (샘플이 있는 경우에만)
  if (!is.null(mcmc_results$tox_sample)) {
    p <- p + geom_ribbon(aes(x = dose_level, ymin = tox_ci[,1], ymax = tox_ci[,2]),
                         fill = "blue", alpha = 0.2)
  }
  
  if (!is.null(mcmc_results$eff_sample)) {
    p <- p + geom_ribbon(aes(x = dose_level, ymin = eff_ci[,1], ymax = eff_ci[,2]),
                         fill = "red", alpha = 0.2)
  }
  
  # 예측값과 관찰값의 차이 계산 (샘플이 있는 경우에만)
  if (!is.null(mcmc_results$tox_sample) & !is.null(mcmc_results$eff_sample)) {
    differences <- data.frame(
      dose_level = dose_levels,
      tox_diff = mcmc_results$tox - observed_tox,
      eff_diff = mcmc_results$eff - observed_eff,
      tox_within_ci = observed_tox >= tox_ci[,1] & observed_tox <= tox_ci[,2],
      eff_within_ci = observed_eff >= eff_ci[,1] & observed_eff <= eff_ci[,2]
    )
  } else {
    differences <- data.frame(
      dose_level = dose_levels,
      tox_diff = mcmc_results$tox - observed_tox,
      eff_diff = mcmc_results$eff - observed_eff
    )
  }
  p <- p + coord_cartesian(ylim = c(0, 1))
  return(list(
    plot = p,
    differences = differences,
    plot_data = plot_data
  ))
}

