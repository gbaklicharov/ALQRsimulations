# Libraries
library(quantreg);library(MASS);library(grf);library(FKSUM);library(SuperLearner);library(MonteCarlo)


# Data-generating mechanism experiment 1, setting 3
generate_cont_data <- function(n){
  sigm <- matrix(c(1,.5,.2,.3,.5,1,.7,0,.2,.7,1,0,.3,0,0,1),nrow = 4,ncol = 4)
  M <- matrix(c(-1,2,2,-1), nrow = 4)  
  L <- mvrnorm(n,rep(0,4),sigm)
  A <- rnorm(n,-0.5-L%*%M,5)     
  mu <- 1 + A + sin(L[,1]) + L[,2]^2 + L[,3] + L[,4] + L[,3]*L[,4]
  Y <- mu + rgamma(n,shape=1,scale=4)
  data <- data.frame(cbind(Y,A,L))
  colnames(data) <- c("Y","A","L1","L2","L3","L4")
  return(data)
}


# Main function
estimate <- function(data,tau,nfolds=5){
  output<-list()
  output[["oracle"]]=list("estimate"=NA,"se"=NA)

  # Method: standard parametric quantile regression
  tryCatch(
    #try to do this
    {
      oracle <- rq(Y~A + I(sin(L1)) + I(L2^2) + L3 + L4 + L3*L4, data = data, tau = tau)
      output[["oracle"]]=list("estimate"=as.numeric(oracle$coefficients[2]),"se"=summary(oracle, se="boot")$coef[2,2])
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
  
  # propensity score model
  sl.a <- SuperLearner(Y=A,X=L,SL.library = SL.library)
  m.A <- predict(sl.a, newdata=L)$pred
  
  # model for Q_tau(Y|A,L)
  qrf <- quantile_forest(X=AL, Y=Y, quantiles = c(tau))
  qAL <- predict(qrf, newdata=AL)$pred
  
  # If A continuous: model for E[q(Y|A,L)|L]
  sl.qAL <- SuperLearner(Y=qAL,X=L,SL.library = SL.library)
  EqAL <- predict(sl.qAL, newdata=L)$pred
  
  # Compute psi for every individual
  psi1 <- (A-m.A)*(qAL-EqAL)/mean((A-m.A)^2)
  
  # Averaging over the whole data set yields estimator for psi
  PlugIn <- mean(psi1)
  
  # EIF for every individual
  EIF1 <- (A-m.A)/mean((A-m.A)^2)*(qAL-EqAL-PlugIn*(A-m.A))
  
  # SE(PlugIn)
  PlugIn_SE <- sd(EIF1)/sqrt(n)
  output[["plugin"]]=list("estimate"=PlugIn,"se"=PlugIn_SE)
  
  
  
  # method: DML without targeting
  
  # Estimate conditional density Y|A,L
  f <- fk_density(Y-qAL, x_eval = rep(0,n))$y
  
  # Compute psi for every individual
  psi3 <- (A-m.A)*(qAL-EqAL+(tau-ifelse(Y<=qAL,1,0))/f)/mean((A-m.A)^2)
  
  # DML estimator for Psi
  DML <- mean(psi3)
  
  # EIF for every individual
  EIF3 <- (A-m.A)/mean((A-m.A)^2)*(qAL-EqAL+(tau-ifelse(Y<=qAL,1,0))/f-DML*(A-m.A))
  
  # SE(DML estimator)
  DML_SE <- sd(EIF3)/sqrt(n)
  output[["DML"]]=list("estimate"=DML,"se"=DML_SE,"max1f"=max(1/f))
  
  
  
  # method: TMLE
  # Updating/targeting the model Q_tau(Y|A,L) (only one step)
  
  w <- (A-m.A)/fk_density(Y-qAL, x_eval = rep(0,n))$y   # weights
  
  fYAL <- function(eps){
    return(sum(w*(tau-ifelse(Y-qAL <= eps*w,1,0))))
  }
  eps <- uniroot(fYAL, c(-1,1))$root
  
  qAL.update <- qAL + eps*w
  
  sl.qAL.update <- SuperLearner(Y=qAL.update, X=L,SL.library = SL.library)
  EqAL.update <- predict(sl.qAL.update, newdata=L)$pred
  
  # Compute psi for every individual
  f5 <- fk_density(Y-qAL.update, x_eval = rep(0,n))$y
  psi5 <- (A-m.A)*(qAL.update-EqAL.update+(tau-ifelse(Y<=qAL.update,1,0))/f5)/mean((A-m.A)^2)
  
  # TMLE estimator for psi
  TMLE <- mean(psi5)
  
  # EIF for every individual
  EIF5 <- (A-m.A)/mean((A-m.A)^2)*(qAL.update-EqAL.update+(tau-ifelse(Y<=qAL.update,1,0))/f5-TMLE*(A-m.A))
  
  # SE(TMLE)
  TMLE_SE <- sd(EIF5)/sqrt(n)
  output[["TMLE"]]=list("estimate"=TMLE,"se"=TMLE_SE,"term"=mean((A-m.A)*(tau-ifelse(Y<=qAL.update,1,0))/f5),"max1f"=max(1/f5))
  
  
  
  #Method: update model but then set additional term exactly equal to zero
  
  # psi for every individual
  psi6 <- (A-m.A)*(qAL.update-EqAL.update)/mean((A-m.A)^2)
  
  # Compute psi averaging over the whole data set
  TMLEbis <- mean(psi6)
  
  # EIF for every individual
  EIF6 <- (A-m.A)/mean((A-m.A)^2)*(qAL.update-EqAL.update+(tau-ifelse(Y<=qAL.update,1,0))/f5-TMLEbis*(A-m.A))
  
  # SE(Psi)
  TMLEbis_SE <- sd(EIF6)/sqrt(n)
  output[["TMLE_bis"]]=list("estimate"=TMLEbis,"se"=TMLEbis_SE)

  
  
  
  # Cross fitting: DML and TMLE
  
  set.seed(123)
  folds <- sample(rep(1:nfolds, length.out = n))
  
  psi3CF <- numeric(n)
  EIF3CF <- numeric(n)
  psi5CF <- numeric(n)
  EIF5CF <- numeric(n)
  
  m.A.I <- numeric(n)
  qAL.I <- numeric(n)
  EqAL.I <- numeric(n)
  Ew.I <- numeric(n)
  f.I <- numeric(n)
  f5CF <- numeric(n)
  
  for(l in 1:nfolds){
    
    # DMl and TMLE with cross-fitting
    
    I <- folds == l
    IC <- folds != l
    
    # fit models for nuisance parameters
    sl.aIC <- SuperLearner(Y=A[IC],X=data.frame(L[IC,]),SL.library = SL.library)
    qrf.IC <- quantile_forest(X=AL[IC,], Y=Y[IC], quantiles = c(tau))
    qAL.ICtrain <- predict(qrf.IC, newdata=AL[IC,])$pred
    
    
    # Make predictions using the fitted model (using model fitted on the other part of the data)
    m.A.I[I] <- predict(sl.aIC, newdata=data.frame(L[I,]))$pred
    qAL.I[I] <- predict(qrf.IC, newdata=AL[I,])$pred
    
    m.A.ICtrain <- predict(sl.aIC, newdata=data.frame(L[IC,]))$pred
    f.ICtrain <- fk_density(Y[IC]-qAL.ICtrain, x_eval = rep(0,sum(IC)))$y
    sl.wIC <- SuperLearner(Y=(A[IC]-m.A.ICtrain)/f.ICtrain, X=L[IC,],SL.library = SL.library)
    Ew.I[I] <- predict(sl.wIC, newdata=L[I,])$pred
    
    sl.qAL.IC <- SuperLearner(Y=qAL.ICtrain,X=L[IC,],SL.library = SL.library)
    EqAL.I[I] <- predict(sl.qAL.IC, newdata=L[I,])$pred
  }
  for(l in 1:nfolds){
    I <- folds == l
    IC <- folds != l
    f.I[I] <- fk_density(Y[IC]-qAL.I[IC], x_eval = rep(0,sum(I)))$y
  }
  # psi for every individual
  psi3CF <- (A-m.A.I)*(qAL.I-EqAL.I+(tau-ifelse(Y<=qAL.I,1,0))/f.I)/mean((A-m.A.I)^2)
  
  # Compute psi1 averaging over the whole data set
  DMLCF <- mean(psi3CF)
  
  # EIF for every individual
  EIF3CF <- (A-m.A.I)/mean((A-m.A.I)^2)*(qAL.I-EqAL.I+(tau-ifelse(Y<=qAL.I,1,0))/f.I-DMLCF*(A-m.A.I))
  
  # SE(DMLCF)
  DMLCF_SE <- sd(EIF3CF)/sqrt(n)
  output[["DML_CF"]]=list("estimate"=DMLCF,"se"=DMLCF_SE,"max1f"=max(1/f.I))
  
  
  # TMLE (one step)
  
  wCF <- (A-m.A.I)/fk_density(Y-qAL.I, x_eval = rep(0,sum(I)))$y   # weights
  
  fYALCF <- function(eps){
    return(sum(wCF*(tau-ifelse(Y-qAL.I <= eps*wCF,1,0))))
  }
  
  eps.CF <- uniroot(fYALCF, c(-1,1))$root
  
  qAL.updateCF <- qAL.I + eps.CF*wCF
  
  # psi for every individual
  for(l in 1:nfolds){
    I <- folds == l
    IC <- folds != l
    f5CF[I] <- fk_density(Y[IC]-qAL.updateCF[IC], x_eval = rep(0,sum(I)))$y
  }
  EqAL.updateCF <- EqAL.I + eps.CF*Ew.I
  psi5CF <- (A-m.A.I)*(qAL.updateCF-EqAL.updateCF+(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF)/mean((A-m.A.I)^2)
  
  # TMLE estimator
  TMLECF <- mean(psi5CF)
  
  # EIF for every individual
  EIF5CF <- (A-m.A.I)/mean((A-m.A.I)^2)*(qAL.updateCF-EqAL.updateCF+(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF-TMLECF*(A-m.A.I))
  
  # SE(TMLECF)
  TMLECF_SE <- sd(EIF5CF)/sqrt(n)
  output[["TMLE_CF"]]=list("estimate"=TMLECF,"se"=TMLECF_SE,"term"=mean((A-m.A.I)*(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF),"max1f"=max(1/f5CF))
  
  
  # method: update model but then set additional term exactly equal to zero
  
  # psi for every individual
  psi6CF <- (A-m.A.I)*(qAL.updateCF-EqAL.updateCF)/mean((A-m.A.I)^2)
  
  # Compute psi averaging over the whole data set
  TMLEbisCF <- mean(psi6CF)
  
  # EIF for every individual
  EIF6CF <- (A-m.A.I)/mean((A-m.A.I)^2)*(qAL.updateCF-EqAL.updateCF+(tau-ifelse(Y<=qAL.updateCF,1,0))/f5CF-TMLEbisCF*(A-m.A.I))
  
  # SE(Psi)
  TMLEbisCF_SE <- sd(EIF6CF)/sqrt(n)
  output[["TMLE_CF_bis"]]=list("estimate"=TMLEbisCF,"se"=TMLEbisCF_SE)
  
  return(output)
}


# Choose one of the data generating functions defined above within this function
MonteCarlo_sim<-function(n,tau){
  data <- generate_cont_data(n)
  output <- estimate(data,tau)
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



