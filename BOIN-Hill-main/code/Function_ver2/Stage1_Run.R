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


