# Libraries
library(quantreg);library(MASS);library(grf);library(FKSUM);library(SuperLearner);library(MonteCarlo)


# Data-generating mechanism experiment 4
generate_hd_data <- function(n){
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  L <- mvrnorm(n,rep(0,50),ar1_cor(50,0.5))
  A <- rnorm(n,L[,1:10]%*%matrix(1/(1:10)))
  Y <- rnorm(n,A+L[,1:5]%*%matrix(1/(1:5))+L[,11:15]%*%matrix(1/(1:5)),2)
  data <- data.frame(cbind(Y,A,L))
  colnames(data) <- c("Y","A",'L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16','L17','L18','L19','L20','L21','L22','L23','L24','L25','L26','L27','L28','L29','L30','L31','L32','L33','L34','L35','L36','L37','L38','L39','L40','L41','L42','L43','L44','L45','L46','L47','L48','L49','L50')
  return(data)
}


# Main function
estimate <- function(data,tau,nfolds=5){
  output<-list()
  output[["QRoracle"]]=list("estimate"=NA,"se"=NA)
  output[["QRvs"]]=list("estimate"=NA,"se"=NA)
  
  # Method: standard parametric quantile regression with stepwise variable selection
  tryCatch(
    #try to do this
    {
      fitoracle <- rq(Y~A+L1+L2+L3+L4+L5+L11+L12+L13+L14+L15, data = data, tau = tau)
      output[["QRoracle"]]=list("estimate"=as.numeric(fitoracle$coefficients["A"]),"se"=summary(fitoracle, se="boot")$coef["A",2])
      model.main <- rq(Y~., data = data, tau = tau)
      fitrq <- step(model.main, trace=FALSE, scope = list(upper=formula(model.main), lower="Y~A"))
      coef <- as.numeric(fitrq$coefficients["A"])
      coef_se <- summary(fitrq, se="boot")$coef["A",2]
      output[["QRvs"]]=list("estimate"=coef,"se"=coef_se)
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
  
  
  # Method: Plug-In estimator (only with quantile regression forests (QRF))
  
  n <- dim(data)[1]
  AL <- data[,-1]
  L <- AL[,-1]
  A <- data[,2]
  Y <- data[,1]
  SL.library <- c("SL.glm", "SL.glmnet" , "SL.randomForest", "SL.gam")
  
  # propensity score model
  sl.a <- SuperLearner(Y=A,X=L,SL.library = SL.library)
  m.A <- predict(sl.a, newdata=L)$pred
  
  # model for Q_tau(Y|A,L) (with QRF)
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
  
  # QRF
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
  
  # DML based on parametric quantile regression with stepwise variable selection
  qALvs <- fitrq$fitted.values
  # Estimate conditional density Y|A,L
  fvs <- fk_density(Y-qALvs, x_eval = rep(0,n))$y
  # Compute psi for every individual
  psi3vs <- coef + (A-m.A)*((tau-ifelse(Y<=qALvs,1,0))/fvs)/mean((A-m.A)^2)
  # DML estimator for Psi
  DMLvs <- mean(psi3vs)
  # EIF for every individual
  EIF3vs <- (A-m.A)/mean((A-m.A)^2)*(coef*(A-m.A)+(tau-ifelse(Y<=qALvs,1,0))/fvs-DMLvs*(A-m.A))
  # SE(DML estimator)
  DMLvs_SE <- sd(EIF3vs)/sqrt(n)
  output[["DML.vs"]]=list("estimate"=DMLvs,"se"=DMLvs_SE,"max1f"=max(1/fvs))
  
  
  
  # method: TMLE with QRF
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
  
  
  
  # method: TMLE with stepwise variable selection
  
  w.vs <- (A-m.A)/fk_density(Y-qALvs, x_eval = rep(0,n))$y   # weights
  
  fYALvs <- function(eps){
    return(sum(w.vs*(tau-ifelse(Y-qALvs <= eps*w.vs,1,0))))
  }
  eps.vs <- uniroot(fYALvs, c(-1,1))$root
  
  qALvs.update <- qALvs + eps.vs*w.vs
  
  sl.qALvs.update <- SuperLearner(Y=qALvs.update, X=L,SL.library = SL.library)
  EqALvs.update <- predict(sl.qALvs.update, newdata=L)$pred
  
  
  # psi for every individual
  f5vs <- fk_density(Y-qALvs.update, x_eval = rep(0,n))$y
  psi5vs <- (A-m.A)*(qALvs.update-EqALvs.update+(tau-ifelse(Y<=qALvs.update,1,0))/f5vs)/mean((A-m.A)^2)
  
  # Compute psi averaging over the whole data set
  TMLEvs <- mean(psi5vs)
  
  # EIF for every individual
  EIF5vs <- (A-m.A)/mean((A-m.A)^2)*(qALvs.update-EqALvs.update+(tau-ifelse(Y<=qALvs.update,1,0))/f5vs-TMLEvs*(A-m.A))
  
  # SE(Psi)
  TMLEvs_SE <- sd(EIF5vs)/sqrt(n)
  output[["TMLEvs"]]=list("estimate"=TMLEvs,"se"=TMLEvs_SE,"term"=mean((A-m.A)*(tau-ifelse(Y<=qALvs.update,1,0))/f5vs),"max1f"=max(1/f5vs))
  
  
  #Method: update model but then set additional term equal to zero
  
  # psi for every individual
  psi6vs <- (A-m.A)*(qALvs.update-EqALvs.update)/mean((A-m.A)^2)
  
  # Compute psi averaging over the whole data set
  TMLEvsbis <- mean(psi6vs)
  
  # EIF for every individual
  EIF6vs <- (A-m.A)/mean((A-m.A)^2)*(qALvs.update-EqALvs.update+(tau-ifelse(Y<=qALvs.update,1,0))/f5vs-TMLEvsbis*(A-m.A))
  
  # SE(Psi)
  TMLEvsbis_SE <- sd(EIF6vs)/sqrt(n)
  output[["TMLEvs_bis"]]=list("estimate"=TMLEvsbis,"se"=TMLEvsbis_SE)
  
  
  
  
  # Cross fitting: DML and TMLE
  
  set.seed(123)
  folds <- sample(rep(1:nfolds, length.out = n))
  
  psi3CF <- numeric(n)
  EIF3CF <- numeric(n)
  psi5CF <- numeric(n)
  EIF5CF <- numeric(n)
  
  psi3CFvs <- numeric(n)
  EIF3CFvs <- numeric(n)
  psi5CFvs <- numeric(n)
  EIF5CFvs <- numeric(n)
  
  m.A.I <- numeric(n)
  qAL.I <- numeric(n)
  EqAL.I <- numeric(n)
  qAL.Ivs <- numeric(n)
  coef.I <- numeric(n)
  f.I <- numeric(n)
  f.Ivs <- numeric(n)
  f5CF <- numeric(n)
  f5CFvs <- numeric(n)
  Ew.I <- numeric(n)
  Ew.Ivs <- numeric(n)
  
  for(l in 1:nfolds){
    
    # DMl and TMLE with cross-fitting
    
    I <- folds == l
    IC <- folds != l
    
    # fit models for nuisance parameters
    sl.aIC <- SuperLearner(Y=A[IC],X=data.frame(L[IC,]),SL.library = SL.library)
    qrf.IC <- quantile_forest(X=AL[IC,], Y=Y[IC], quantiles = c(tau))
    qAL.ICtrain <- predict(qrf.IC, newdata=AL[IC,])$pred
    
    model.mainIC <- rq(Y~., data = data[IC,], tau = tau)
    fitrqIC <- step(model.mainIC, trace=FALSE, scope = list(upper=formula(model.mainIC), lower="Y~A"))
    qAL.ICtrainvs <- fitrqIC$fitted.values 
    
    
    # Make predictions using the fitted model (using model fitted on the other part of the data)
    m.A.I[I] <- predict(sl.aIC, newdata=data.frame(L[I,]))$pred
    qAL.I[I] <- predict(qrf.IC, newdata=AL[I,])$pred
    sl.qAL.IC <- SuperLearner(Y=qAL.ICtrain,X=L[IC,],SL.library = SL.library)
    EqAL.I[I] <- predict(sl.qAL.IC, newdata=L[I,])$pred
    
    qAL.Ivs[I] <- predict(fitrqIC, newdata=AL[I,])
    coef.I[I] <- as.numeric(fitrqIC$coefficients["A"])
    
    m.A.ICtrain <- predict(sl.aIC, newdata=data.frame(L[IC,]))$pred
    f.ICtrain <- fk_density(Y[IC]-qAL[IC], x_eval = rep(0,sum(IC)))$y
    f.ICtrainvs <- fk_density(Y[IC]-qALvs[IC], x_eval = rep(0,sum(IC)))$y
    sl.wIC <- SuperLearner(Y=(A[IC]-m.A.ICtrain)/f.ICtrain, X=L[IC,],SL.library = SL.library)
    Ew.I[I] <- predict(sl.wIC, newdata=L[I,])$pred
    sl.wICvs <- SuperLearner(Y=(A[IC]-m.A.ICtrain)/f.ICtrainvs, X=L[IC,],SL.library = SL.library)
    Ew.Ivs[I] <- predict(sl.wICvs, newdata=L[I,])$pred
  }
  for(l in 1:nfolds){
    I <- folds == l
    IC <- folds != l
    f.I[I] <- fk_density(Y[IC]-qAL.I[IC], x_eval = rep(0,sum(I)))$y
    f.Ivs[I] <- fk_density(Y[IC]-qAL.Ivs[IC], x_eval = rep(0,sum(I)))$y
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
  
  # DML variable selection
  # psi for every individual
  psi3CFvs <- coef.I + (A-m.A.I)*((tau-ifelse(Y<=qAL.Ivs,1,0))/f.Ivs)/mean((A-m.A.I)^2)
  # Compute psi1 averaging over the whole data set
  DMLvs_CF <- mean(psi3CFvs)
  # EIF for every individual
  EIF3CFvs <- (A-m.A.I)/mean((A-m.A.I)^2)*(coef.I*(A-m.A.I)+(tau-ifelse(Y<=qAL.Ivs,1,0))/f.Ivs-DMLvs_CF*(A-m.A.I))
  # SE(Psi3CF)
  DMLvs_CF_SE <- sd(EIF3CFvs)/sqrt(n)
  output[["DML_CF.vs"]]=list("estimate"=DMLvs_CF,"se"=DMLvs_CF_SE,"max1f"=max(1/f.Ivs))
  
  
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
  
  
  # TMLE variable selection
  wCFvs <- (A-m.A.I)/fk_density(Y-qAL.Ivs, x_eval = rep(0,sum(I)))$y   # weights
  fYALCFvs <- function(eps){
    return(sum(wCFvs*(tau-ifelse(Y-qAL.Ivs <= eps*wCFvs,1,0))))
  }
  eps.CFvs <- uniroot(fYALCFvs, c(-1,1))$root
  qALvs.updateCF <- qAL.Ivs + eps.CFvs*wCFvs
  # psi for every individual
  for(l in 1:nfolds){
    I <- folds == l
    IC <- folds != l
    f5CFvs[I] <- fk_density(Y[IC]-qALvs.updateCF[IC], x_eval = rep(0,sum(I)))$y
  }
  psi5CFvs <- coef.I + eps.CFvs/f5CFvs + (A-m.A.I)*(-eps.CFvs*Ew.I+(tau-ifelse(Y<=qALvs.updateCF,1,0))/f5CFvs)/mean((A-m.A.I)^2)
  # Compute psi averaging over the whole data set
  TMLEvs_CF <- mean(psi5CFvs)
  # EIF for every individual
  EIF5CFvs <- (A-m.A.I)/mean((A-m.A.I)^2)*((coef.I-eps.CFvs*Ew.I)*(A-m.A.I)+(tau-ifelse(Y<=qALvs.updateCF,1,0))/f5CFvs-TMLEvs_CF*(A-m.A.I))
  # SE(Psi)
  TMLEvs_CF_SE <- sd(EIF5CFvs)/sqrt(n)
  output[["TMLEvs_CF"]]=list("estimate"=TMLEvs_CF,"se"=TMLEvs_CF_SE,"term"=mean((A-m.A.I)*(tau-ifelse(Y<=qALvs.updateCF,1,0))/f5CFvs),"max1f"=max(1/f5CFvs))
  
  
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
  
  # psi for every individual
  psi6CFvs<- coef.I + eps.CFvs/f5CFvs - (A-m.A.I)*(-eps.CFvs*Ew.I)/mean((A-m.A.I)^2)
  # Compute psi averaging over the whole data set
  TMLEvsbisCF <- mean(psi6CFvs)
  # EIF for every individual
  EIF6CFvs <- (A-m.A.I)/mean((A-m.A.I)^2)*((coef.I-eps.CFvs*Ew.I)*(A-m.A.I)+(tau-ifelse(Y<=qALvs.updateCF,1,0))/f5CFvs-TMLEvsbisCF*(A-m.A.I))
  # SE(Psi)
  TMLEvsbisCF_SE <- sd(EIF6CFvs)/sqrt(n)
  output[["TMLEvs_CF_bis"]]=list("estimate"=TMLEvsbisCF,"se"=TMLEvsbisCF_SE)
  
  
  return(output)
}


# Choose one of the data generating functions defined above within this function
MonteCarlo_sim<-function(n,tau){
  data <- generate_hd_data(n)
  output <- estimate(data,tau)
  ## What the function will return
  return(list("QRoracle_est"=output$QRoracle$estimate, "QRoracle_SE"=output$QRoracle$se, "QRvs_est"=output$QRvs$estimate, "QRvs_SE"=output$QRvs$se,
              "PlugIn_est"=output$plugin$estimate, "PlugIn_SE"=output$plugin$se,
              "DML_est"=output$DML$estimate, "DML_SE"=output$DML$se, "DML_max1f" = output$DML$max1f, "DML_CF_est"=output$DML_CF$estimate,
              "DML_CF_SE"=output$DML_CF$se, "DML_CF_max1f"=output$DML_CF$max1f,
              "DML.vs_est"=output$DML.vs$estimate, "DML.vs_SE"=output$DML.vs$se, "DML.vs_max1f" = output$DML.vs$max1f, "DML.vs_CF_est"=output$DML_CF.vs$estimate,
              "DML.vs_CF_SE"=output$DML_CF.vs$se, "DML.vs_CF_max1f"=output$DML_CF.vs$max1f, "TMLE_est"=output$TMLE$estimate, "TMLE_SE"=output$TMLE$se,
              "TMLE_max1f"=output$TMLE$max1f, "TMLE_term"=output$TMLE$term, "TMLE_CF_est"=output$TMLE_CF$estimate, "TMLE_CF_SE"=output$TMLE_CF$se,
              "TMLE_CF_max1f"=output$TMLE_CF$max1f, "TMLE_CF_term"=output$TMLE_CF$term, "TMLE_bis_est"=output$TMLE_bis$estimate, "TMLE_bis_se"=output$TMLE_bis$se,
              "TMLE_CF_bis_est"=output$TMLE_CF_bis$estimate, "TMLE_CF_bis_se"=output$TMLE_CF_bis$se,
              "TMLEvs_est"=output$TMLEvs$estimate, "TMLEvs_SE"=output$TMLEvs$se,
              "TMLEvs_max1f"=output$TMLEvs$max1f, "TMLEvs_term"=output$TMLEvs$term, "TMLEvs_CF_est"=output$TMLEvs_CF$estimate, "TMLEvs_CF_SE"=output$TMLEvs_CF$se,
              "TMLEvs_CF_max1f"=output$TMLEvs_CF$max1f, "TMLEvs_CF_term"=output$TMLEvs_CF$term, "TMLEvs_bis_est"=output$TMLEvs_bis$estimate, "TMLEvs_bis_se"=output$TMLEvs_bis$se,
              "TMLEvs_CF_bis_est"=output$TMLEvs_CF_bis$estimate, "TMLEvs_CF_bis_se"=output$TMLEvs_CF_bis$se))
}


it <- 1000
n <- c(250,500)
tau <- c(0.5,0.75,0.9)
set.seed(123)

# Using the MonteCarlo package/function, the simulation can be run in parallel with ncpus the number of CPU's
MC_result<-MonteCarlo(MonteCarlo_sim,nrep=it,param_list = list("n"=n,"tau"=tau), ncpus = 16)
summary(MC_result)



