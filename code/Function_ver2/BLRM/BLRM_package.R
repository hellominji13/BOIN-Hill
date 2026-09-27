## prior 모델이랑 모델 잘 맞추기
## NAN 값이 계속 생성됨 (error는 아니고 fitting 이 잘 안됨 문제 -> 사용시 고쳐두기)


blrm_mono_ss <- function(para,prior, data, current_dose, output_excel=FALSE, output_pdf=FALSE){
  ################### generic functions to be used ###################
  # probability of toxicity
  tox_prob <- function(d, para){  
    # d: a numeric scalar
    # para: the original form of alpha & beta
    # Once we have (posterior) estimates of parameters, we are able to calculate probability of toxicity for each provisional dose
    
    logit_pi1 <- para[1] + exp(para[2]) * d
    pi1 <- exp(logit_pi1)/(1 + exp(logit_pi1))
    return(list(pi1=pi1))
  }
  
  # prior for the log(alpha) and log(beta)
  prior_summary <- function(prior){ # `prior` is a list composed of three elements, and 
    # each element is a numeric vector containing information 
    # for log(alpha) and log(beta), respectively, and their correlation
    mu <- prior[[1]]
    se <- prior[[2]]
    corr <- prior[[3]]
    covariance_matrix <- matrix(c(se[1]**2, se[1]*se[2]*corr, se[1]*se[2]*corr, se[2]**2), ncol=2, byrow=T)
    return(list(mean=mu, var=covariance_matrix))
  }
  
  # go through cumulated observed cohorts
  single_cohort <- function(prior_para){
    
    # data for JAGS
    mydata <- list(nb_pat=sum(cohort_size), 
                   s=data_dlt,    
                   d=data_sdose,
                   p1=prior_para$mean,
                   p2=prior_para$var)
    
    niters <- (1+burn_in)*nsamples/2 # number of iterations, e.g., (1+0.2)*10000, since 20% of 10000 will be chopped off
    
    # justify the model for library(rjags)
    # modelstring <- "model
    # {
    # Omega1[1:2, 1:2] <- inverse(p2[, ])
    # log.alpha[1:2] ~ dmnorm(p1[], Omega1[, ])
    # alpha[1] <- exp(log.alpha[1])
    # alpha[2] <- exp(log.alpha[2])
    # for(i in 1:nb_pat){
    # logit(pi1[i]) <- log.alpha[1] + alpha[2] * d[i]
    # s[i] ~ dbern(pi1[i])
    # }
    # }"
    modelstring <- "model
    {
    Omega1[1:2, 1:2] <- inverse(p2[, ])
    alpha[1:2] ~ dmnorm(p1[], Omega1[, ])

    for(i in 1:nb_pat){
        logit(pi1[i]) <- alpha[1] + exp(alpha[2]) * d[i]
        s[i] ~ dbern(pi1[i])
    }
    }"
    # initialize the parameters in the log form, and set up the RNG seeds
    inits.list <- list(list(alpha=c(1,0), .RNG.seed=seeds[1], .RNG.name="base::Wichmann-Hill"),
                  list(alpha=c(1,0), .RNG.seed=seeds[2], .RNG.name="base::Wichmann-Hill"))
    jagsobj <- rjags::jags.model(textConnection(modelstring), data=mydata, n.chains=2, quiet=TRUE, inits=inits.list)
    update(jagsobj, n.iter=niters, progress.bar="none")
    res <- rjags::jags.samples(jagsobj, "alpha", n.iter=niters, progress.bar="none") # `alpha`: alpha & beta
    # chop off the burn-in iterations, e.g., 20% of the target number of samples (e.g., 10000)
    # for each parameter; there are 2 chains generated in parallel, and thus, 10% of each will be chopped off
    alpha1 <- res$alpha[1,,][-c(1:(burn_in*nsamples/2)),]
    beta1 <- res$alpha[2,,][-c(1:(burn_in*nsamples/2)),]
    # the posterior samples of interest
    posterior_param <- cbind(c(alpha1), c(beta1)) 
    
    # (posterior) parameter estimates: log(alpha.hat) & log(beta.hat)
    log_para <- posterior_param
    para_hat <- apply(log_para, 2, mean) # point estimates: alpha.hat and beta.hat
    para_sd <- apply(log_para, 2, sd)    # standard deviation of the point estimates
    para_corr <- boot::corr(log_para)          # correlation coefficient between alpha.hat and beta.hat
    posterior_para_summary <- list(para_hat=para_hat, para_sd=para_sd, para_corr=para_corr)
    
    # collect posterior estimates for pi1, i.e., Pr(DLT | data): for each pair of 
    # posterior parameter estimates, we have a posterior estimate of DLT rate (probability of toxicity)
    samples_sdose <- matrix(NA, nrow=nprov_dose, ncol=nsamples) # go through each standardized provisional dose
    for(i in 1:nprov_dose){
      for(j in 1:nsamples){
        samples_sdose[i,j] <- tox_prob(sprov_dose[i], posterior_param[j,])$pi1  
      }
    }
    
    # interval probabilities by dose: obtain the probability that Pr(DLT | data) lies in each category interval
    posterior_prob_summary <- matrix(NA, nrow=length(category_bound)+1, ncol=nprov_dose)
    for(i in 1:nprov_dose){
      posterior_prob_summary[,i] <- as.numeric(table(cut(samples_sdose[i,], breaks=c(0, category_bound[1], category_bound[2], 1), right=TRUE))/nsamples)
    }
    
    # (posterior) pi1 estimates: mean & sd and 0.025, 0.5 & 0.975 quantile
    posterior_pi_summary <- matrix(0, 5, nprov_dose)
    for(i in 1:nprov_dose){
      posterior_pi_summary[,i] <- c(mean(samples_sdose[i,]), sd(samples_sdose[i,]), quantile(samples_sdose[i,], c(0.025, 0.5, 0.975)))
    }
    
    # the next dose: select the dose with the highest probability of targeted-toxicity among safe doses
    ## justify safe doses first
    safe_dose_range <- (posterior_prob_summary[3,] <= ewoc) # posterior_prob_summary[3,]: (0.33, 1]\
    
    
    #조건 추가
    if (data$n_pat[current_dose] > para$cohort) {
      # Allow safe_dose_range to include up to current_dose + 1
      safe_dose_range <- safe_dose_range & (seq_along(safe_dose_range) <= current_dose + 1)
    } else {
      # Restrict safe_dose_range to current_dose only
      safe_dose_range <- safe_dose_range & (seq_along(safe_dose_range) <= current_dose)
    }
    # safe_dose_range <- safe_dose_range & (seq_along(safe_dose_range) <= current_dose + 1)
    
    
    
    if(sum(safe_dose_range) != 0){
      
      ## subset safe doses
      safe_dose <- prov_dose[safe_dose_range]
      safe_dose_prob <- posterior_prob_summary[,safe_dose_range]
      safe_dose_pi <- posterior_pi_summary[,safe_dose_range]
      # justify the dose with the highest probability of targeted toxicity
      if (is.vector(safe_dose_prob)) {
        next_dose_index <- which(safe_dose_range)
        next_dose_posterior_prob_summary <- safe_dose_prob
        next_dose_posterior_pi_summary <- safe_dose_pi
      }else {
        next_dose_index <- which.max(safe_dose_prob[2,]) 
        next_dose_posterior_prob_summary <- safe_dose_prob[,next_dose_index] # interval probabilities by dose
        next_dose_posterior_pi_summary <- safe_dose_pi[,next_dose_index]     # posterior summary of DLT rate
      }
      next_dose_level <- safe_dose[next_dose_index]     
      
      # the recommended dose level
      # combine them to a list
      next_dose <- list(dose=next_dose_level, 
                        posterior_prob_summary=list(next_dose_posterior_prob_summary), 
                        posterior_pi_summary=list(next_dose_posterior_pi_summary))     
    }else{
      next_dose <- NULL
    }
    # return the posterior summary & the recommended next dose
    return(list(posterior_pi_summary=posterior_pi_summary,
                posterior_prob_summary=posterior_prob_summary,
                posterior_para_summary=posterior_para_summary, 
                next_dose=next_dose))
  }
  ## end of single cohort
  
  # extract prior
  prior_para <- prior_summary(prior)
  
  # extract parameters
  seeds <- data$seeds
  nsamples <- data$nsamples 
  burn_in <- data$burn_in
  category_bound <- data$category_bound
  category_name <- data$category_name
  ewoc <- data$ewoc
  
  # extract drug information
  drug_name <- data$drug_name
  prov_dose <- data$prov_dose
  # order provisional doses from low to high (for purpose of index tracking)
  prov_dose <- sort(prov_dose, decreasing=FALSE)
  ref_dose <- data$ref_dose
  dose_unit <- data$dose_unit
  
  # extract observed cohorts
  dose <- data$dose
  n_pat <- data$n_pat
  dlt_test <- data$dlt
  
  # further process drug-related information
  sdose <- log(dose/ref_dose)           # log(scaled tested doses)
  ndose <- length(sdose)                # number of tested doses
  
  sprov_dose <- log(prov_dose/ref_dose) # log(scaled provisional doses)
  nprov_dose <- length(sprov_dose)      # number of provisional doses
  
  # go through all the tested doses cumulatively
  current_dose_index <- 1
  dose_index <- cohort_size <- dlt <- NULL
  data_sdose <- data_dlt <- NULL
  
  # save the results for each tested dose as a list
  cohort_all <- list()
  
  
  # Combine all data first
  for(current_dose_index in 1:ndose) {
    current_cohort_size <- n_pat[current_dose_index]
    current_dlt <- dlt_test[current_dose_index]
    
    # Accumulate data
    dose_index <- c(dose_index, current_dose_index)
    cohort_size <- c(cohort_size, current_cohort_size)
    dlt <- c(dlt, current_dlt)
    
    current_data_sdose <- rep(sdose[current_dose_index], current_cohort_size)
    data_sdose <- c(data_sdose, current_data_sdose)
    
    current_data_dlt <- ifelse(rep(current_dlt==0, current_cohort_size),
                               rep(0, current_cohort_size),
                               c(rep(1, current_dlt), rep(0, current_cohort_size-current_dlt)))
    data_dlt <- c(data_dlt, current_data_dlt)
  }
  
  # Run analysis with complete data
  cohort <- single_cohort(prior_para)
  
  
  if(is.null(cohort$next_dose)) {
    stop("There is no RND even not all the tested doses are gone through!")
  }
  cohort_final <- list(
    posterior_prob = cohort$posterior_prob_summary,
    posterior_para = cohort$posterior_para_summary,
    posterior_pi = cohort$posterior_pi_summary,
    next_dose = cohort$next_dose
  )
  
  # while(current_dose_index <= ndose){
  # 
  #   current_cohort_size <- n_pat[current_dose_index]
  #   current_dlt <- dlt_test[current_dose_index]
  # 
  #   # cumulatively
  #   dose_index <- c(dose_index, current_dose_index)
  #   cohort_size <- c(cohort_size, current_cohort_size)
  #   dlt <- c(dlt, current_dlt)
  # 
  #   # for each cohort, statistical analysis is based on all the cumulated observed cohorts information
  #   current_data_sdose <- rep(sdose[current_dose_index], current_cohort_size)
  #   data_sdose <- c(data_sdose, current_data_sdose)
  #   current_data_dlt <- ifelse(rep(current_dlt==0, current_cohort_size), rep(0, current_cohort_size), c(rep(1, current_dlt), rep(0, current_cohort_size-current_dlt)))
  #   data_dlt <- c(data_dlt, current_data_dlt)
  # 
  #   # go through the observed cohorts cumulatively
  #   cohort <- single_cohort(prior_para)
  #   cohort_all[[current_dose_index]] <- list(posterior_prob=cohort$posterior_prob_summary,
  #                                            posterior_para=cohort$posterior_para_summary,
  #                                            posterior_pi=cohort$posterior_pi_summary,
  #                                            next_dose=cohort$next_dose)
  # 
  #   if(is.null(cohort$next_dose)){ # no recommended next dose
  #     if(current_dose_index < ndose){
  #       stop("There is no RND even not all the tested doses are gone through!")
  #     }else{  # current_dose_index == ndose
  #       break
  #     }
  #   }else{ # continue to the next tested dose
  #     current_dose_index <- current_dose_index + 1
  #   }
  # }
  # # done with going through all the tested doses
  # print(cohort_all)
  # # we only take interest to the final-round information
  # cohort_final <- cohort_all[[length(cohort_all)]]
  
  
  # extract posterior summary and next dose
  prob_posterior <- cohort_final$posterior_prob
  pi_posterior <- cohort_final$posterior_pi
  para_posterior <- cohort_final$posterior_para
  next_dose <- cohort_final$next_dose
  
  # construct the framework for the simulation output
  return(list(prob_posterior=prob_posterior,    # interval probabilities by dose
              para_posterior=para_posterior,    # posterior parameter estimates
              pi_posterior=pi_posterior,        # posterior DLT rate estimates
              next_dose=next_dose              # recommended next dose
             ))           
}