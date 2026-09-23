troisPtrois.kim <- function (para,data = data, lastdose) 
{
  idx <- which(data$dose == lastdose)
  ndlt <- data$ndlt[idx]
  npt <- data$npt[idx]
  npt.previous<-NA; if(idx>1) {npt.previous <- data$npt[idx-1]}
  npt.next <- data$npt[min(para$J,idx+1)]
  
  mtd <- NA
  if (npt == 3) {
    nextdose <- ifelse(ndlt == 0, lastdose + 1, ifelse(ndlt > 1, NA, lastdose))
  }
  if (npt == 6) {
    nextdose <- ifelse(ndlt <= 1, lastdose + 1, NA)
  }

  if (is.na(nextdose) & npt.previous==6  & !is.na(npt.previous)) {nextdose<-NA; mtd<-lastdose-1}
  if (is.na(nextdose) & npt.previous < 6 & !is.na(npt.previous)) {nextdose<-lastdose-1; mtd=NA}
  if (is.na(nextdose) & is.na(npt.previous)) {nextdose<-NA; mtd=NA} 
  
  if (!is.na(nextdose) & nextdose == lastdose+1 & npt.next>0) {nextdose<-NA; mtd=lastdose} 
  if (!is.na(nextdose) & nextdose > max(data$dose) & npt==6) {nextdose<-NA; mtd<-lastdose} 
  if (!is.na(nextdose) & nextdose > max(data$dose) & npt<6)  {nextdose<-lastdose; mtd=NA} 
  list(nextdose = nextdose, mtd = mtd)
  
}

CreData <- function (ndose = 3, dosenames = paste("dose", 1:ndose, sep = " ")) 
{
  data <- data.frame(dose = 1:ndose, npt = rep(0, ndose), ndlt = rep(0, 
                                                                     ndose), row.names = dosenames)
  return(data)
}
dose3P3 <- function(para,outcome_table,time_table,T_prior=NULL){
  dt.tox <- outcome_table[grep("Tox", rownames(outcome_table)),]
  
  data <- CreData(nrow(dt.tox))
  nextdose <- 1
  index <- 1:3
  trace<-NULL
  while (nextdose %in% data$dose) {
    lastdose <- nextdose
    ndlt <- sum(dt.tox[lastdose,index]); 
    temp <- paste(paste("D",lastdose,"S", index,sep=""),collapse=" "); trace<-paste(trace,temp,collapse=" ")
    index<-index+3
    data <- updata(data, lastdose, 3, ndlt)
    nextdose <- troisPtrois.kim(para,data, lastdose)$nextdose
    mtd      <- troisPtrois.kim(para,data, lastdose)$mtd
    
    if (max(data$npt) > para$S1) {
      break
    }
  }
  list(obd = mtd, n = min(index) - 1, trace = trace, data = data)
}



UBLRM.stage1 <- function(para,dt.tox,T_prior){
  dt.tox <- outcome_table[grep("Tox", rownames(outcome_table)),]
  data <- CreData(nrow(dt.tox))
  nextdose <- 1
  index <- 1:para$cohort
  trace<-NULL
  
  while (nextdose %in% data$dose) {
    lastdose <- nextdose
    ndlt <- sum(dt.tox[1:para$J == lastdose, index])
    temp <- paste(paste("D", lastdose, "S", index, sep = ""), collapse = " ")
    trace <- paste(trace, temp, collapse = " ")
    index <- index + para$cohort
    data <- updata(data, lastdose, para$cohort, ndlt)
    
    #BLRM Data 입력, Dose assign 할때 category 기준으로 posterior 값을 나누어서 
    #ewoc는 over dose의 확률이 0.25 미만이 되도록 stopping rule 
    BLRM_data <- list(seeds=sample(1:1000, 2), nsamples=1000, burn_in=0.2, drug_name="test",dose_unit="mg", prov_dose=para$doses, ref_dose=para$ref_dose,
                      dose=para$doses, n_pat=data$npt, dlt=data$ndlt, category_bound=c(para$T_admissible_low, para$T_admissible),
                      category_name=c("under-dosing", "targeted-toxicity", "over-dosing"), ewoc = 0.9)
    
    # MCMC 실행 
    trial <- tryCatch({
      blrm_mono_ss(para=para, prior = T_prior, data = BLRM_data, current_dose =nextdose, output_excel = FALSE, output_pdf = FALSE)
    }, error = function(e) {
      if (grepl("There is no RND even not all the tested doses are gone through!", e$message)) {
        return(NULL)  # Return NULL if the specific error is caught
      } else {
        stop(e)  # Rethrow other errors
      }
    })    
    
    # 중단 규칙 추가로 
    if (is.null(trial)) { #만약 선택된 next dose가 없을경우, dose 그대로 유지 
      nextdose <- nextdose
      #return(list(obd = NA, n = min(index) - 1, trace = trace, data = data)) 또는 중단 (이경우는 논문에서 사용안한듯 )
    }
    else {
      # prior update 하기 
      #names(trial$para_posterior) <- c("mean", "se", "corr")
      #T_prior <- trial$para_posterior 
      
      # Next dose assign
      # next_dose_index <- which max(targeted-toxicity_prob)
      nextdose <- which(para$doses == trial$next_dose$dose)
    }
    
    
    if (max(data$npt) >= para$S1 ) {
      break
    }
    
  }
  obd  <- nextdose
  list(obd = obd, n = min(index) - 1, trace = trace, data = data)
}

  

UBOIN.stage1 <- function(para, outcome_table,time_table,T_prior=NULL) {
  dt.tox <- outcome_table[grep("Tox", rownames(outcome_table)),]
  data <- CreData(nrow(dt.tox))
  nextdose <- 1
  index <- 1:para$cohort
  trace <- NULL
  dose_range <- data$dose
  while (nextdose %in% dose_range) {
    boundary <- get.boundary(target = para$DLT, ncohort = (para$N / para$cohort), cohortsize = para$cohort)
    lastdose <- nextdose
    ndlt <- sum(dt.tox[which(1:para$J == lastdose), index])
    temp <- paste(paste("D", lastdose, "S", index, sep = ""), collapse = " ")
    trace <- paste(trace, temp, collapse = " ")
    index <- index + para$cohort
    data <- updata(data, lastdose, para$cohort, ndlt)
    
    m_j <- data$ndlt[lastdose]
    n_j <- data$npt[lastdose]
    # Next dose assign 
    if (m_j <= boundary$boundary_tab["Escalate if # of DLT <=",n_j/para$cohort]) {
      nextdose <- min(lastdose + 1, para$J)  # 증가
      }
    else if (m_j >= boundary$boundary_tab["Deescalate if # of DLT >=",n_j/para$cohort]){
      nextdose <- max(lastdose - 1, 1)  # 감소
      if (is.na(boundary$boundary_tab["Eliminate if # of DLT >=", n_j/para$cohort])) {
        boundary$boundary_tab["Eliminate if # of DLT >=", n_j/para$cohort] <- boundary$boundary_tab["Number of patients treated", 1]
      }
      if (m_j >= boundary$boundary_tab["Eliminate if # of DLT >=", n_j/para$cohort] && n_j>=para$cohort ) {
        dose_range <- dose_range[dose_range < lastdose]
        if (length(dose_range) == 0){
          mtd <- 99
          break
        }
      }
    }

    else if (m_j < boundary$boundary_tab["Deescalate if # of DLT >=", n_j/para$cohort] &&
             m_j > boundary$boundary_tab["Escalate if # of DLT <=", n_j/para$cohort])
      {nextdose <- lastdose}  # 유지
    if (max(data$npt) >= para$S1) {break}
  }
  # reach S1, stop
  
  result <- select.mtd(target = para$DLT, npts = data$npt, ntox = data$ndlt)
  mtd <- result$MTD
  #early stop 경우 정리
  if (mtd==99){
    mtd <- NA 
  }
  
  return(list(obd = mtd, n = min(index) - 1, trace = trace, data = data))
}


