# Libraries
library(quantreg);library(MASS);library(grf);library(FKSUM);library(SuperLearner);library(MonteCarlo);library(parTreat)

# Data-generating mechanism experiment 1, setting 1
generate_data <- function(n){
  sigm <- matrix(c(1,.5,.2,.3,.5,1,.7,0,.2,.7,1,0,.3,0,0,1),nrow = 4,ncol = 4)
  M <- matrix(0.2*c(-1,2,2,-1), nrow = 4)    
  L <- mvrnorm(n,rep(0,4),sigm)
  A <- rbinom(n,1,1/(1+exp(0.5+L%*%M))) 
  mu <- 1 + A + sin(L[,1]) + L[,2]^2 + L[,3] + L[,4] + L[,3]*L[,4]
  Y <- mu + rgamma(n,shape=1,scale=2)
  data <- data.frame(cbind(Y,A,L))
  colnames(data) <- c("Y","A","L1","L2","L3","L4")
  return(data)
}

# Data-generating mechanism experiment 1, setting 2
generate_heterosc_data <- function(n){
  sigm <- matrix(c(1,.5,.2,.3,.5,1,.7,0,.2,.7,1,0,.3,0,0,1),nrow = 4,ncol = 4)
  M <- matrix(0.2*c(-1,2,2,-1), nrow = 4)    
  L <- mvrnorm(n,rep(0,4),sigm)
  A <- rbinom(n,1,1/(1+exp(0.5+L%*%M))) 
  mu <- 1 + A + sin(L[,1]) + L[,2]^2 + L[,3] + L[,4] + L[,3]*L[,4]
  Y <- mu + rgamma(n,shape=1,scale=2+A)
  data <- data.frame(cbind(Y,A,L))
  colnames(data) <- c("Y","A","L1","L2","L3","L4")
  return(data)
}

# Data-generating mechanism experiment 2
generate_complex_data <- function(n){
  sigm <- matrix(c(1,.5,.2,.3,.5,1,.7,0,.2,.7,1,0,.3,0,0,1),nrow = 4,ncol = 4)
  M <- matrix(0.2*c(-1,2,2,-1), nrow = 4)    
  L <- mvrnorm(n,rep(0,4),sigm)
  A <- rbinom(n,1,1/(1+exp(0.5+L%*%M+0.5*L[,1]^2-0.5*L[,2]^2+0.5*L[,3]*L[,4]))) 
  mu <- 1 + A + sin(L[,1]) + L[,2]^2 + L[,3] + L[,4] + L[,3]*L[,4]
  Y <- mu + rgamma(n,shape=1,scale=3)
  data <- data.frame(cbind(Y,A,L))
  colnames(data) <- c("Y","A","L1","L2","L3","L4")
  return(data)
}

# Data-generating mechanism experiment 3
generate_randtreat_data <- function(n){
  sigm <- matrix(c(1,.5,.2,.3,.5,1,.7,0,.2,.7,1,0,.3,0,0,1),nrow = 4,ncol = 4)
  M <- matrix(0.2*c(-1,2,2,-1), nrow = 4)    
  
  L <- mvrnorm(n,rep(0,4),sigm)
  A <- rbinom(n,1,0.5)     
  mu <- 1 + A + sin(L[,1]) + L[,2]^2 + L[,3] + L[,4] + L[,3]*L[,4]
  Y <- mu + rgamma(n,shape=1,scale=2)
  data <- data.frame(cbind(Y,A,L))
  colnames(data) <- c("Y","A","L1","L2","L3","L4")
  #ggplot(data, aes(x=Y, color=as.factor(A)))  + geom_density()
  return(data)
}

# Main function
estimate <- function(data,tau){
  output<-list()
  output[["oracle"]]=list("estimate"=NA,"se"=NA)
  output[["QR"]]=list("estimate"=NA,"se"=NA)
  
  # Method: standard quantile regression
  tryCatch(
    #try to do this
    {
      oracle <- rq(Y~A + I(sin(L1)) + I(L2^2) + L3 + L4 + L3*L4, data = data, tau = tau)
      fitrq <- rq(Y~A + L1 + L2 + L3 + L4, data = data, tau = tau)
      output[["oracle"]]=list("estimate"=as.numeric(oracle$coefficients[2]),"se"=summary(oracle, se="boot")$coef[2,2])
      coef <- as.numeric(fitrq$coefficients[2])
      coef_se <- summary(fitrq, se="boot")$coef[2,2]
      output[["QR"]]=list("estimate"=coef,"se"=coef_se)
    },
    #if an error occurs, tell me the error
    error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
    }
  )
  
  
  # Method: Plug-In estimator
  
  n <- dim(data)[1]
  AL <- data[,-1]
  L <- AL[,-1]
  A <- data[,2]
  Y <- data[,1]
  SL.library <- c("SL.glm", "SL.glmnet" , "SL.randomForest", "SL.gam")
  sl.a <- SuperLearner(Y=A,X=L,SL.library = SL.library,family = binomial())
  m.A <- predict(sl.a, newdata=L)$pred
  
  # model for Q_tau(Y|A,L)
  qrf <- quantile_forest(X=AL, Y=Y, quantiles = c(tau))
  qAL <- predict(qrf, newdata=AL)$pred
  
  # if A is binary, E[Q_tau(Y|A,L)|L] = E[Q_tau(Y|A=1,L)|L]E(A|L) + E[Q_tau(Y|A=0,L)|L](1-E(A|L))
  A1L <- AL
  A1L$A <- 1
  A0L <- AL
  A0L$A <- 0
  
  # Predict Q_tau(Y|A=1,L) and Q_tau(Y|A=0,L)
  q1L <- predict(qrf, newdata=A1L)$pred
  q0L <- predict(qrf, newdata=A0L)$pred
  
  # Compute psi for every individual
  psi1 <- (A-m.A)*(qAL-q1L*m.A-q0L*(1-m.A))/mean((A-m.A)^2)
  
  # Averaging over the whole data set yields estimator for psi
  PlugIn <- mean(psi1)
  
  # EIF for every individual
  EIF1 <- (A-m.A)/mean((A-m.A)^2)*(qAL-q1L*m.A-q0L*(1-m.A)-PlugIn*(A-m.A))
  
  # SE(PlugIn)
  PlugIn_SE <- sd(EIF1)/sqrt(n)
  output[["plugin"]]=list("estimate"=PlugIn,"se"=PlugIn_SE)
  
  
  
  # method: DML without targeting
  
  # Estimate conditional density Y|A,L
  f <- fk_density(Y-qAL, x_eval = rep(0,n))$y
  
  # Compute psi for every individual
  psi3 <- (A-m.A)*(qAL-q1L*m.A-q0L*(1-m.A)+(tau-ifelse(Y<=qAL,1,0))/f)/mean((A-m.A)^2)
  
  # DML estimator for Psi
  DML <- mean(psi3)
  
  # EIF for every individual
  EIF3 <- (A-m.A)/mean((A-m.A)^2)*(qAL-q1L*m.A-q0L*(1-m.A)+(tau-ifelse(Y<=qAL,1,0))/f-DML*(A-m.A))
  
  # SE(DML estimator)
  DML_SE <- sd(EIF3)/sqrt(n)
  output[["DML"]]=list("estimate"=DML,"se"=DML_SE,"max1f"=max(1/f))
  
  
  
  # method: TMLE
  # Updating/targeting the model Q_tau(Y|A,L) (iteratively, using the initial model)
  eps.updates <- rep(0,20)
  w <- matrix(0,nrow = 20,ncol=n)
  w0 <- matrix(0,nrow = 20,ncol=n)
  w1 <- matrix(0,nrow = 20,ncol=n)
  k<-1
  crit <- c(100000,99999)
  while (k<=20 & crit[k+1]<crit[k]) {
    eps.w <- as.vector(eps.updates%*%w)
    eps.w0 <- as.vector(eps.updates%*%w0)
    eps.w1 <- as.vector(eps.updates%*%w1)
    qAL.update <- qAL + eps.w
    
    w[k,] <- (A-m.A)/fk_density(Y-qAL.update, x_eval = rep(0,n))$y
    w0[k,] <- -m.A/fk_density(Y-q0L-as.vector(eps.updates%*%w0), x_eval = rep(0,n))$y
    w1[k,] <- (1-m.A)/fk_density(Y-q1L-as.vector(eps.updates%*%w1), x_eval = rep(0,n))$y
    
    r <- Y-qAL.update                             # residuals
    
    # What we want to make zero:
    fYAL <- function(eps){
      return(sum(w[k,]*(tau-ifelse(r <= eps*w[k,],1,0))))
    }
    
    eps.updates[k] <- uniroot(fYAL, c(-1,1))$root    # Approximate root
    
    # Keeping track of the absolute value of the term we want to make zero.
    # Once this term is not decreasing anymore, we exit the loop
    crit[k+2] <- abs(fYAL(eps.updates[k]))
    k = k+1
  }
  k <- k-2
  
  eps.w <- as.vector(eps.updates[1:k]%*%w[1:k,])
  eps.w0 <- as.vector(eps.updates[1:k]%*%w0[1:k,])
  eps.w1 <- as.vector(eps.updates[1:k]%*%w1[1:k,])
  qAL.update <- qAL + eps.w
  q1L.update <- q1L + eps.w1
  q0L.update <- q0L + eps.w0
  
  # Compute psi for every individual
  f5 <- fk_density(Y-qAL.update, x_eval = rep(0,n))$y # evaluate density in updated quantile
  psi5 <- (A-m.A)*(qAL.update-q1L.update*m.A-q0L.update*(1-m.A)+(tau-ifelse(Y<=qAL.update,1,0))/f5)/mean((A-m.A)^2)
  
  # TMLE estimator for psi
  TMLE <- mean(psi5)
  
  # EIF for every individual
  EIF5 <- (A-m.A)/mean((A-m.A)^2)*(qAL.update-q1L.update*m.A-q0L.update*(1-m.A)+(tau-ifelse(Y<=qAL.update,1,0))/f5-TMLE*(A-m.A))
  
  # SE(TMLE)
  TMLE_SE <- sd(EIF5)/sqrt(n)
  output[["TMLE"]]=list("estimate"=TMLE,"se"=TMLE_SE,"term"=mean((A-m.A)*(tau-ifelse(Y<=qAL.update,1,0))/f5),"max1f"=max(1/f5))
  
  
  
  #Method: update model but then set additional term exactly equal to zero
  
  # psi for every individual
  psi6 <- (A-m.A)*(qAL.update-q1L.update*m.A-q0L.update*(1-m.A))/mean((A-m.A)^2)
  
  # Compute psi averaging over the whole data set
  TMLEbis <- mean(psi6)
  
  # EIF for every individual
  EIF6 <- (A-m.A)/mean((A-m.A)^2)*(qAL.update-q1L.update*m.A-q0L.update*(1-m.A)+(tau-ifelse(Y<=qAL.update,1,0))/f5-TMLEbis*(A-m.A))
  
  # SE(Psi)
  TMLEbis_SE <- sd(EIF6)/sqrt(n)
  output[["TMLE_bis"]]=list("estimate"=TMLEbis,"se"=TMLEbis_SE)
  
  # # 1 step TMLE
  # qAL.1step <- qAL + as.vector(eps.updates[1]%*%w[1,])
  # q1L.1step <- q1L + as.vector(eps.updates[1]%*%w0[1,])
  # q0L.1step <- q0L + as.vector(eps.updates[1]%*%w1[1,])
  # f5.1step <- fk_density(Y-qAL.1step, x_eval = rep(0,n))$y
  # psi5.1step <- (A-m.A)*(qAL.1step-q1L.1step*m.A-q0L.1step*(1-m.A)+(tau-ifelse(Y<=qAL.1step,1,0))/f5.1step)/mean((A-m.A)^2)
  # TMLE.1step <- mean(psi5.1step)
  # EIF5.1step <- (A-m.A)/mean((A-m.A)^2)*(qAL.1step-q1L.1step*m.A-q0L.1step*(1-m.A)+(tau-ifelse(Y<=qAL.1step,1,0))/f5.1step-TMLE.1step*(A-m.A))
  # TMLE_SE.1step <- sd(EIF5.1step)/sqrt(n)
  # output[["TMLE.1step"]]=list("estimate"=TMLE.1step,"se"=TMLE_SE.1step,"term"=mean((A-m.A)*(tau-ifelse(Y<=qAL.1step,1,0))/f5.1step),"max1f"=max(1/f5.1step))
  # psi6.1step <- (A-m.A)*(qAL.1step-q1L.1step*m.A-q0L.1step*(1-m.A))/mean((A-m.A)^2)
  # TMLEbis.1step <- mean(psi6.1step)
  # EIF6.1step <- (A-m.A)/mean((A-m.A)^2)*(qAL.1step-q1L.1step*m.A-q0L.1step*(1-m.A)+(tau-ifelse(Y<=qAL.1step,1,0))/f5.1step-TMLEbis.1step*(A-m.A))
  # TMLEbis_SE.1step <- sd(EIF6.1step)/sqrt(n)
  # output[["TMLEbis.1step"]]=list("estimate"=TMLEbis.1step,"se"=TMLEbis_SE.1step)
  
  
  # Cross fitting: DML and TMLE
  
  I1 <- sort(sample(1:n,n/5))
  I2 <- sort(sample(setdiff(1:n,I1),n/5))
  I3 <- sort(sample(setdiff(1:n,c(I1,I2)),n/5))
  I4 <- sort(sample(setdiff(1:n,c(I1,I2,I3)),n/5))
  I5 <- setdiff(1:n,c(I1,I2,I3,I4))
  II <- cbind(I1,I2,I3,I4,I5)
  
  psi3CF <- numeric(n)
  EIF3CF <- numeric(n)
  psi5CF <- numeric(n)
  EIF5CF <- numeric(n)
  
  m.A.I <- numeric(n)
  qAL.I <- numeric(n)
  q1L.I <- numeric(n)
  q0L.I <- numeric(n)
  qAL.update.I <- numeric(n)
  q1L.update.I <- numeric(n)
  q0L.update.I <- numeric(n)
  f.I <- numeric(n)
  f5CF <- numeric(n)
  
  for(l in 1:5){
    
    # DMl and TMLE with cross-fitting
    
    I <- II[,l]
    IC <- setdiff(1:n,I)
    
    # fit models for nuisance parameters
    sl.aIC <- SuperLearner(Y=A[IC],X=data.frame(L[IC,]),newX=data.frame(L[I,]),SL.library = SL.library,family = binomial())
    qrf.IC <- quantile_forest(X=AL[IC,], Y=Y[IC], quantiles = c(tau))
    
    
    # Make predictions using the fitted model (using model fitted on the other part of the data)
    m.A.I[I] <- predict(sl.aIC, newdata=data.frame(L[I,]))$pred
    qAL.I[I] <- predict(qrf.IC, newdata=AL[I,])$pred
    q1L.I[I] <- predict(qrf.IC, newdata=A1L[I,])$pred
    q0L.I[I] <- predict(qrf.IC, newdata=A0L[I,])$pred
  }
  for(l in 1:5){
    I <- II[,l]
    IC <- setdiff(1:n,I)
    f.I[I] <- fk_density(Y[IC]-qAL.I[IC], x_eval = rep(0,n/5))$y
  }
  # psi for every individual
  psi3CF <- (A-m.A.I)*(qAL.I-q1L.I*m.A.I-q0L.I*(1-m.A.I)+(tau-ifelse(Y<=qAL.I,1,0))/f.I)/mean((A-m.A.I)^2)
  
  # Compute psi1 averaging over the whole data set
  DMLCF <- mean(psi3CF)
  
  # EIF for every individual
  EIF3CF <- (A-m.A.I)/mean((A-m.A.I)^2)*(qAL.I-q1L.I*m.A.I-q0L.I*(1-m.A.I)+(tau-ifelse(Y<=qAL.I,1,0))/f.I-DMLCF*(A-m.A.I))
  
  # SE(DMLCF)
  DMLCF_SE <- sd(EIF3CF)/sqrt(n)
  output[["DML_CF"]]=list("estimate"=DMLCF,"se"=DMLCF_SE,"max1f"=max(1/f.I))
  
  
  # TMLE
  
  eps.updatesCF <- rep(0,20)
  wCF <- matrix(0,nrow = 20,ncol=n)
  w0CF <- matrix(0,nrow = 20,ncol=n)
  w1CF <- matrix(0,nrow = 20,ncol=n)
  k<-1
  crit <- c(100000,99999)
  while (k<=20 & crit[k+1]<crit[k]){
    eps.wCF <- as.vector(eps.updatesCF%*%wCF)
    eps.w0CF <- as.vector(eps.updatesCF%*%w0CF)
    eps.w1CF <- as.vector(eps.updatesCF%*%w1CF)
    qAL.updateCF <- qAL.I + eps.wCF
    
    wCF[k,] <- (A-m.A.I)/fk_density(Y-qAL.updateCF, x_eval = rep(0,n))$y    
    w0CF[k,] <- -m.A.I/fk_density(Y-q0L.I-as.vector(eps.updatesCF%*%w0CF), x_eval = rep(0,n))$y
    w1CF[k,] <- (1-m.A.I)/fk_density(Y-q1L.I-as.vector(eps.updatesCF%*%w1CF), x_eval = rep(0,n))$y
    
    rCF <- Y-qAL.updateCF                             # residuals
    
    fYALCF <- function(eps){
      return(sum(wCF[k,]*(tau-ifelse(rCF <= eps*wCF[k,],1,0))))
    }
    
    eps.updatesCF[k] <- uniroot(fYALCF, c(-1,1))$root       
    
    # Keeping track of the absolute value of the term we want to make zero.
    # Once this term is not decreasing anymore, we exit the loop
    crit[k+2] <- abs(fYALCF(eps.updates[k]))
    k = k+1
  }
  k <- k-2
  
  eps.wCF <- as.vector(eps.updatesCF[1:k]%*%wCF[1:k,])
  eps.w0CF <- as.vector(eps.updatesCF[1:k]%*%w0CF[1:k,])
  eps.w1CF <- as.vector(eps.updatesCF[1:k]%*%w1CF[1:k,])
  qAL.updateCF <- qAL.I + eps.wCF
  q1L.updateCF <- q1L.I + eps.w1CF
  q0L.updateCF <- q0L.I + eps.w0CF
  
  # psi for every individual
  for(l in 1:5){
    I <- II[,l]
    IC <- setdiff(1:n,I)
    f5CF[I] <- fk_density(Y[IC]-qAL.updateCF[IC], x_eval = rep(0,n/5))$y
  }
  psi5CF <- (A-m.A.I)*(qAL.updateCF-q1L.updateCF*m.A.I-q0L.updateCF*(1-m.A.I)+(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF)/mean((A-m.A.I)^2)
  
  # TMLE estimator
  TMLECF <- mean(psi5CF)
  
  # EIF for every individual
  EIF5CF <- (A-m.A.I)/mean((A-m.A.I)^2)*(qAL.updateCF-q1L.updateCF*m.A.I-q0L.updateCF*(1-m.A.I)+(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF-TMLECF*(A-m.A.I))
  
  # SE(TMLECF)
  TMLECF_SE <- sd(EIF5CF)/sqrt(n)
  output[["TMLE_CF"]]=list("estimate"=TMLECF,"se"=TMLECF_SE,"term"=mean((A-m.A.I)*(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF),"max1f"=max(1/f5CF))
  
  
  # method: update model but then set additional term exactly equal to zero
  
  # psi for every individual
  psi6CF <- (A-m.A.I)*(qAL.updateCF-q1L.updateCF*m.A.I-q0L.updateCF*(1-m.A.I))/mean((A-m.A.I)^2)
  
  # Compute psi averaging over the whole data set
  TMLEbisCF <- mean(psi6CF)
  
  # EIF for every individual
  EIF6CF <- (A-m.A.I)/mean((A-m.A.I)^2)*(qAL.updateCF-q1L.updateCF*m.A.I-q0L.updateCF*(1-m.A.I)+(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF-TMLEbisCF*(A-m.A.I))
  
  # SE(Psi)
  TMLEbisCF_SE <- sd(EIF6CF)/sqrt(n)
  output[["TMLE_CF_bis"]]=list("estimate"=TMLEbisCF,"se"=TMLEbisCF_SE)
  
  
  # # 1 step TMLE
  # qAL.1stepCF <- qAL.I + as.vector(eps.updatesCF[1]%*%wCF[1,])
  # q1L.1stepCF <- q1L.I + as.vector(eps.updatesCF[1]%*%w0CF[1,])
  # q0L.1stepCF <- q0L.I + as.vector(eps.updatesCF[1]%*%w1CF[1,])
  # f5.1stepCF <- fk_density(Y-qAL.1step, x_eval = rep(0,n))$y
  # psi5.1stepCF <- (A-m.A)*(qAL.1stepCF-q1L.1stepCF*m.A-q0L.1stepCF*(1-m.A)+(tau-ifelse(Y<=qAL.1stepCF,1,0))/f5.1stepCF)/mean((A-m.A)^2)
  # TMLE.1stepCF <- mean(psi5.1stepCF)
  # EIF5.1stepCF <- (A-m.A)/mean((A-m.A)^2)*(qAL.1stepCF-q1L.1stepCF*m.A-q0L.1stepCF*(1-m.A)+(tau-ifelse(Y<=qAL.1stepCF,1,0))/f5.1stepCF-TMLE.1stepCF*(A-m.A))
  # TMLE_SE.1stepCF <- sd(EIF5.1stepCF)/sqrt(n)
  # output[["TMLE.1stepCF"]]=list("estimate"=TMLE.1stepCF,"se"=TMLE_SE.1stepCF,"term"=mean((A-m.A)*(tau-ifelse(Y<=qAL.1stepCF,1,0))/f5.1stepCF),"max1f"=max(1/f5.1stepCF))
  # psi6.1stepCF <- (A-m.A)*(qAL.1stepCF-q1L.1stepCF*m.A-q0L.1stepCF*(1-m.A))/mean((A-m.A)^2)
  # TMLEbis.1stepCF <- mean(psi6.1stepCF)
  # EIF6.1stepCF <- (A-m.A)/mean((A-m.A)^2)*(qAL.1stepCF-q1L.1stepCF*m.A-q0L.1stepCF*(1-m.A)+(tau-ifelse(Y<=qAL.1stepCF,1,0))/f5.1stepCF-TMLEbis.1stepCF*(A-m.A))
  # TMLEbis_SE.1stepCF <- sd(EIF6.1stepCF)/sqrt(n)
  # output[["TMLEbis.1stepCF"]]=list("estimate"=TMLEbis.1stepCF,"se"=TMLEbis_SE.1stepCF)
  
  return(output)
}


# Choose one of the data generating functions defined above within this function
MonteCarlo_sim<-function(n,tau){
  data <- generate_data(n)
  output <- estimate(data,tau)
  
  # Add the following 3 lines if you also want to use the methods of Athey et al. (2021) and add them to the output list by adding: "eifad_est"=imbens$eifad$tau, "eifad_SE"=imbens$eifad$se, "waq_est"=imbens$waq$tau, "waq_SE"=imbens$waq$se
  # y1 <- data$Y[which(data$A==1)]
  # y0 <- data$Y[which(data$A==0)]
  # athey <- list("eifad"=eif_additive(y0,y1),"waq"=waq(y0,y1))
  
  ## What the function will return
  return(list("oracle_est"=output$oracle$estimate, "oracle_se"=output$oracle$se, "QR_est"=output$QR$estimate, "QR_SE"=output$QR$se, "PlugIn_est"=output$plugin$estimate, "PlugIn_SE"=output$plugin$se,
              "DML_est"=output$DML$estimate, "DML_SE"=output$DML$se, "DML_max1f" = output$DML$max1f, "DML_CF_est"=output$DML_CF$estimate,
              "DML_CF_SE"=output$DML_CF$se, "DML_CF_max1f"=output$DML_CF$max1f, "TMLE_est"=output$TMLE$estimate, "TMLE_SE"=output$TMLE$se,
              "TMLE_max1f"=output$TMLE$max1f, "TMLE_term"=output$TMLE$term, "TMLE_CF_est"=output$TMLE_CF$estimate, "TMLE_CF_SE"=output$TMLE_CF$se,
              "TMLE_CF_max1f"=output$TMLE_CF$max1f, "TMLE_CF_term"=output$TMLE_CF$term, "TMLE_bis_est"=output$TMLE_bis$estimate, "TMLE_bis_se"=output$TMLE_bis$se,
              "TMLE_CF_bis_est"=output$TMLE_CF_bis$estimate, "TMLE_CF_bis_se"=output$TMLE_CF_bis$se))
}


it <- 1000
n <- c(500,1000)
tau <- c(0.5,0.75,0.9)
set.seed(123)

# Using the MonteCarlo package/function, the simulation can be run in parallel with ncpus the number of CPU's
MC_result<-MonteCarlo(MonteCarlo_sim,nrep=it,param_list = list("n"=n,"tau"=tau), ncpus = 16)
summary(MC_result)



