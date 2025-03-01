library(caret)
library(dplyr)


#### Delta = 0.5 ####
R <- 100
valtabRMSEs <- array(NA, dim = c(R, 2, 3))
valtabGDEVs <- array(NA, dim = c(R, 2, 3))
parmtabs    <- array(NA, dim = c(R, 4, 3))

tablist_delta0.5 <- list(valtabRMSEs=valtabRMSEs, valtabGDEVs=valtabGDEVs,
                         parmtabs=parmtabs)


for (r in 1:100) {
  
  Tt = 5
  I = 5000
  psi = 1
  a10 = 3
  Delta = 0.5
  set.seed(r)
  ap <- a10
  bp <- a10
  theta1 <- rgamma(I, shape=ap+1, rate=bp)   
  mu1    <- runif( I, 2000, 4000)
  nu1    <- rpois( I, 0.4) + 1
  Y1     <- rgamma(I, shape=nu1/psi, scale=mu1/theta1*psi)
  a      <- (ap + nu1/psi)
  b      <- (bp + Y1/mu1/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta2 <- rgamma(I, shape=ap+1, rate=bp)   
  mu2    <- runif( I, 2000, 4000)
  nu2    <- rpois( I, 0.6) + rbinom(I, 1, 0.8)
  Y2     <- rgamma(I, shape=nu2/psi, scale=mu2/theta2*psi)   
  a      <- (ap + nu2/psi)
  b      <- (bp + Y2/mu2/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta3 <- rgamma(I, shape=ap+1, rate=bp)   
  mu3    <- runif( I, 2000, 4000)
  nu3    <- rpois( I, 0.8) + rbinom(I, 1, 0.6)
  Y3     <- rgamma(I, shape=nu3/psi, scale=mu3/theta3*psi)   
  a      <- (ap + nu3/psi)
  b      <- (bp + Y3/mu3/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta4 <- rgamma(I, shape=ap+1, rate=bp)   
  mu4    <- runif( I, 2000, 4000)
  nu4    <- rpois( I, 0.5) + rbinom(I, 1, 0.4)
  Y4     <- rgamma(I, shape=nu4/psi, scale=mu4/theta4*psi)   
  a      <- (ap + nu4/psi)
  b      <- (bp + Y4/mu4/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta5 <- rgamma(I, shape=ap+1, rate=bp)   
  mu5    <- runif( I, 2000, 4000)
  nu5    <- rpois( I, 1.2) + rbinom(I, 1, 0.2)
  Y5     <- rgamma(I, shape=nu5/psi, scale=mu5/theta5*psi)   
  a      <- (ap + nu5/psi)
  b      <- (bp + Y5/mu5/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta6 <- rgamma(I, shape=ap+1, rate=bp)   
  mu6    <- runif( I, 2000, 4000)
  nu6    <- rpois( I, 1.4)
  Y6     <- rgamma(I, shape=nu6/psi, scale=mu6/theta6*psi)   
  
  
  
  mdata  <-
    data.frame(paste0("P", rep(formatC(1:I, width = 4, format = "d", flag = "0"), each=Tt+1)),
               rep(1:6, I),
               as.vector(t(cbind(nu1, nu2, nu3, nu4, nu5, nu6))),
               as.vector(t(cbind(mu1, mu2, mu3, mu4, mu5, mu6))),
               as.vector(t(cbind(Y1 , Y2 , Y3 , Y4 , Y5 , Y6))))
  colnames(mdata) <- c("ID", "Time", "nu", "mu", "Y")
  mdata$Ypermu    <- mdata$Y / mdata$mu
  
  mdata_train <- mdata[mdata$Time != 6,]
  mdata_test  <- mdata[mdata$Time == 6,]
  
  hmmean <- sum(mdata_train$Y) / sum(mdata_train$nu)
  
  dGP <- function(y, k, a, bmp, log=FALSE) {
    temp   <- lgamma(k+a+1) - lgamma(k) -lgamma(a+1) +
      k*log(y / (y+bmp)) + (a+1)*log(bmp / (y+bmp)) - log(y)
    result <- (log)*temp + (!log)*exp(temp)
    return(result) }
  
  ggSM <- function(a10, psi, delta, heterogenous = TRUE, skip = FALSE) {
    
    M  <- length(unique(mdata_train$ID))
    Tt <- nrow(mdata_train)/length(unique(mdata_train$ID))
    Y  <- mdata_train$Y
    v  <- mdata_train$nu
    mu <- (  heterogenous   * mdata_train$mu   + 
               (-heterogenous+1)*(mdata_train$mu>0)*hmmean          )
    Ys <- mdata_train$Ypermu*(  heterogenous   * 1 + 
                                  (-heterogenous+1)*mdata_train$mu/hmmean)
    # standardized Y by mu
    a  <- rep(NA, nrow(mdata_train))
    b  <- rep(NA, nrow(mdata_train))
    q  <- rep(NA, nrow(mdata_train))
    ap <- rep(NA, nrow(mdata_train))
    bp <- rep(NA, nrow(mdata_train))
    
    for (i in 1:nrow(mdata_train)) {
      
      ap[i] <- if (i %% Tt == 1 || sum(v[((i-1) %/% Tt *Tt +1):(i-1)])==0) {a10}
      else if (v[i-1]==0 && skip==TRUE) {a[i-1]} 
      else {q[i-1]    /delta   *a[i-1]}
      bp[i] <- if (i %% Tt == 1 || sum(v[((i-1) %/% Tt *Tt +1):(i-1)])==0) {a10}  
      else if (v[i-1]==0 && skip==TRUE) {b[i-1]} 
      else {q[i-1]*((1/delta-1)*a[i-1]+b[i-1])}
      # for each policyholder, ap[i] and bp[i] are reset at the first period
      # (assuming that the dataframe is sorted by ID and HY)
      
      a[i] <- ap[i] + v[i]/psi 
      b[i] <- bp[i] + Ys[i]/psi
      
      q[i] <- delta*a10 / (a[i] - delta^2*a[i] + delta^2*a10) 
      
    }
    neglik <- -sum(dGP(Y, v/psi, ap, bp*mu*psi, log=TRUE)[v>0])
    return(list(a=a, b=b, q=q, ap=ap, bp=bp, neglik=neglik)) }
  
  homo_BM_neglik <- function(parm) { # only 2 parameters, delta = 1 by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- 1
    result <- ggSM(a10, psi, delta, heterogenous = FALSE)$neglik
    return(result) }
  
  homo_ScM_neglik <- function(parm) { # 3 parameters, delta in [0,1] by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- parm[3]
    result <- ggSM(a10, psi, delta, heterogenous = FALSE, skip = FALSE)$neglik
    return(result) }
  
  hetero_BM_neglik <- function(parm) { # only 2 parameters, delta = 1 by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- 1
    result <- ggSM(a10, psi, delta, heterogenous = TRUE)$neglik
    return(result) }
  
  hetero_ScM_neglik <- function(parm) { # 3 parameters, delta in [0,1] by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- parm[3]
    result <- ggSM(a10, psi, delta, heterogenous = TRUE, skip = FALSE)$neglik
    return(result) }
  
  
  homo_BM_fit <- optim(c(1, 1), homo_BM_neglik, method="L-BFGS-B",
                       lower = c(0,0), upper = c(Inf, Inf))
  
  homo_ScM_fit <- optim(c(homo_BM_fit$par/2, 1), homo_ScM_neglik, method="L-BFGS-B",
                        lower = c(0,0,0), upper = c(Inf, Inf,1))
  
  hetero_BM_fit <- optim(c(1,1), hetero_BM_neglik, method="L-BFGS-B",
                         lower = c(0,0), upper = c(Inf, Inf))
  
  hetero_ScM_fit <- optim(c(hetero_BM_fit$par/2, 1), hetero_ScM_neglik, method="L-BFGS-B",
                          lower = c(0,0,0), upper = c(Inf, Inf, 1))
  
  homo_BM_fit$value
  homo_ScM_fit$value
  hetero_BM_fit$value
  hetero_ScM_fit$value
  
  endindex <- (1:nrow(mdata_train)) %% Tt == 0
  
  # Crediblity factors for t+1:   b_{t+1|t}/a_{t+1|t}             (from (2.4))
  #                            = (b_{t}/a_{t}-1) * Delta_t + 1    (from (3.6))
  
  homo_BM_cred <-(ggSM(a10=homo_BM_fit$par[1], psi=homo_BM_fit$par[2], 
                       delta=1, heterogenous=FALSE)$b[endindex]/
                    ggSM(a10=homo_BM_fit$par[1], psi=homo_BM_fit$par[2], 
                         delta=1, heterogenous=FALSE)$a[endindex])
  
  hetero_BM_cred <-(ggSM(a10=hetero_BM_fit$par[1], psi=hetero_BM_fit$par[2], 
                         delta=1, heterogenous=TRUE)$b[endindex]/
                      ggSM(a10=hetero_BM_fit$par[1], psi=hetero_BM_fit$par[2], 
                           delta=1, heterogenous=TRUE)$a[endindex])
  
  homo_ScM_cred <-(ggSM(a10=homo_ScM_fit$par[1], psi=homo_ScM_fit$par[2], 
                        delta=homo_ScM_fit$par[3], heterogenous=FALSE, skip = FALSE)$b[endindex]/
                     ggSM(a10=homo_ScM_fit$par[1], psi=homo_ScM_fit$par[2], 
                          delta=homo_ScM_fit$par[3], heterogenous=FALSE, skip = FALSE)$a[endindex]-1)*
    homo_ScM_fit$par[3] + 1
  
  hetero_ScM_cred <-(ggSM(a10=hetero_ScM_fit$par[1], psi=hetero_ScM_fit$par[2], 
                          delta=hetero_ScM_fit$par[3], heterogenous=TRUE, skip = FALSE)$b[endindex]/
                       ggSM(a10=hetero_ScM_fit$par[1], psi=hetero_ScM_fit$par[2], 
                            delta=hetero_ScM_fit$par[3], heterogenous=TRUE, skip = FALSE)$a[endindex]-1)*
    hetero_ScM_fit$par[3] + 1
  
  dzd <- aggregate(mdata_train[,c(3,5,6)], by=list(mdata_train$ID), FUN = sum)
  
  homo_CM_cred <- (dzd$Y/hmmean + homo_BM_fit$par[1]*homo_BM_fit$par[2]) /
    (dzd$nu       + homo_BM_fit$par[1]*homo_BM_fit$par[2])
  
  hetero_CM_cred <- (dzd$Ypermu + homo_BM_fit$par[1]*homo_BM_fit$par[2]) /
    (dzd$nu     + homo_BM_fit$par[1]*homo_BM_fit$par[2])
  
  
  zzzz <- data.frame(mdata_test$ID, homo_BM_cred, hetero_BM_cred, 
                     homo_ScM_cred, hetero_ScM_cred)
  colnames(zzzz)[1] <- "ID"
  
  zdz <- merge(x = zzzz, y = mdata_test,
               by = "ID", all.x = TRUE)
  
  
  valtable_RMSE <- 
    rbind(
      c(RMSE(zdz$Y, hmmean                       * zdz$nu),
        RMSE(zdz$Y, hmmean * zdz$homo_BM_cred    * zdz$nu),
        RMSE(zdz$Y, hmmean * zdz$homo_ScM_cred   * zdz$nu)),
      c(RMSE(zdz$Y, zdz$mu                       * zdz$nu),
        RMSE(zdz$Y, zdz$mu * zdz$hetero_BM_cred  * zdz$nu),
        RMSE(zdz$Y, zdz$mu * zdz$hetero_ScM_cred * zdz$nu)))
  
  colnames(valtable_RMSE) <- c("Indep", "Buhlmann","ScM")
  rownames(valtable_RMSE) <- c("Homogenous", "Heterogeneous")
  
  
  
  GDEV <- function(y, mu, nu) {
    vec    <- -nu*log(y/mu/nu) + (y-mu*nu)/mu
    result <- 2*sum( vec[nu>0] ) }
  
  valtable_GDEV <- 
    rbind(
      c(GDEV(zdz$Y, hmmean                      , zdz$nu),
        GDEV(zdz$Y, hmmean * zdz$homo_BM_cred   , zdz$nu),
        GDEV(zdz$Y, hmmean * zdz$homo_ScM_cred  , zdz$nu)),
      c(GDEV(zdz$Y, zdz$mu                      , zdz$nu),
        GDEV(zdz$Y, zdz$mu * zdz$hetero_BM_cred , zdz$nu),
        GDEV(zdz$Y, zdz$mu * zdz$hetero_ScM_cred, zdz$nu)))
  
  colnames(valtable_GDEV) <- c("Indep", "Buhlmann", "ScM")
  rownames(valtable_GDEV) <- c("Homogenous", "Heterogeneous")
  
  valtable_RMSE
  valtable_GDEV
  
  parmtable_dep <- round(rbind(c(
    homo_BM_fit$par , 1),
    homo_ScM_fit$par,
    c( hetero_BM_fit$par , 1),
    hetero_ScM_fit$par), 4 )
  
  colnames(parmtable_dep) <- c("$a_{10}$", "$\\psi$", "$\\Delta$")
  rownames(parmtable_dep) <- c("Homo_BM"  , "Homo_Scm"  ,
                               "Hetero_BM", "Hetero_Scm")
  
  print(r)
  tablist_delta0.5$valtabRMSEs[r,,] <- valtable_RMSE
  tablist_delta0.5$valtabGDEVs[r,,] <- valtable_GDEV
  tablist_delta0.5$parmtabs[   r,,] <- parmtable_dep
}



#### Delta = 1.0 ####
R <- 100
valtabRMSEs <- array(NA, dim = c(R, 2, 3))
valtabGDEVs <- array(NA, dim = c(R, 2, 3))
parmtabs    <- array(NA, dim = c(R, 4, 3))

tablist_delta1.0 <- list(valtabRMSEs=valtabRMSEs, valtabGDEVs=valtabGDEVs,
                          parmtabs=parmtabs)


for (r in 1:100) {
  
  Tt = 5
  I = 5000
  psi = 1
  a10 = 3
  Delta = 1.0
  set.seed(r)
  ap <- a10
  bp <- a10
  theta1 <- rgamma(I, shape=ap+1, rate=bp)   
  mu1    <- runif( I, 2000, 4000)
  nu1    <- rpois( I, 0.4) + 1
  Y1     <- rgamma(I, shape=nu1/psi, scale=mu1/theta1*psi)
  a      <- (ap + nu1/psi)
  b      <- (bp + Y1/mu1/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta2 <- rgamma(I, shape=ap+1, rate=bp)   
  mu2    <- runif( I, 2000, 4000)
  nu2    <- rpois( I, 0.6) + rbinom(I, 1, 0.8)
  Y2     <- rgamma(I, shape=nu2/psi, scale=mu2/theta2*psi)   
  a      <- (ap + nu2/psi)
  b      <- (bp + Y2/mu2/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta3 <- rgamma(I, shape=ap+1, rate=bp)   
  mu3    <- runif( I, 2000, 4000)
  nu3    <- rpois( I, 0.8) + rbinom(I, 1, 0.6)
  Y3     <- rgamma(I, shape=nu3/psi, scale=mu3/theta3*psi)   
  a      <- (ap + nu3/psi)
  b      <- (bp + Y3/mu3/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta4 <- rgamma(I, shape=ap+1, rate=bp)   
  mu4    <- runif( I, 2000, 4000)
  nu4    <- rpois( I, 1.0) + rbinom(I, 1, 0.4)
  Y4     <- rgamma(I, shape=nu4/psi, scale=mu4/theta4*psi)   
  a      <- (ap + nu4/psi)
  b      <- (bp + Y4/mu4/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta5 <- rgamma(I, shape=ap+1, rate=bp)   
  mu5    <- runif( I, 2000, 4000)
  nu5    <- rpois( I, 1.2) + rbinom(I, 1, 0.2)
  Y5     <- rgamma(I, shape=nu5/psi, scale=mu5/theta5*psi)   
  a      <- (ap + nu5/psi)
  b      <- (bp + Y5/mu5/psi)
  q      <- Delta*a10 / (a-Delta^2*a+Delta^2*a10)
  
  ap     <- q*( 1/Delta  )*a
  bp     <- q*((1/Delta-1)*a+b)
  theta6 <- rgamma(I, shape=ap+1, rate=bp)   
  mu6    <- runif( I, 2000, 4000)
  nu6    <- rpois( I, 1.4)
  Y6     <- rgamma(I, shape=nu6/psi, scale=mu6/theta6*psi)   
  
                 
                   
  mdata  <-
    data.frame(paste0("P", rep(formatC(1:I, width = 4, format = "d", flag = "0"), each=Tt+1)),
               rep(1:6, I),
               as.vector(t(cbind(nu1, nu2, nu3, nu4, nu5, nu6))),
               as.vector(t(cbind(mu1, mu2, mu3, mu4, mu5, mu6))),
               as.vector(t(cbind(Y1 , Y2 , Y3 , Y4 , Y5 , Y6))))
  colnames(mdata) <- c("ID", "Time", "nu", "mu", "Y")
  mdata$Ypermu    <- mdata$Y / mdata$mu
  
  mdata_train <- mdata[mdata$Time != 6,]
  mdata_test  <- mdata[mdata$Time == 6,]
  
  hmmean <- sum(mdata_train$Y) / sum(mdata_train$nu)
  
  dGP <- function(y, k, a, bmp, log=FALSE) {
    temp   <- lgamma(k+a+1) - lgamma(k) -lgamma(a+1) +
      k*log(y / (y+bmp)) + (a+1)*log(bmp / (y+bmp)) - log(y)
    result <- (log)*temp + (!log)*exp(temp)
    return(result) }
  
  ggSM <- function(a10, psi, delta, heterogenous = TRUE, skip = FALSE) {
    
    M  <- length(unique(mdata_train$ID))
    Tt <- nrow(mdata_train)/length(unique(mdata_train$ID))
    Y  <- mdata_train$Y
    v  <- mdata_train$nu
    mu <- (  heterogenous   * mdata_train$mu   + 
               (-heterogenous+1)*(mdata_train$mu>0)*hmmean          )
    Ys <- mdata_train$Ypermu*(  heterogenous   * 1 + 
                                  (-heterogenous+1)*mdata_train$mu/hmmean)
    # standardized Y by mu
    a  <- rep(NA, nrow(mdata_train))
    b  <- rep(NA, nrow(mdata_train))
    q  <- rep(NA, nrow(mdata_train))
    ap <- rep(NA, nrow(mdata_train))
    bp <- rep(NA, nrow(mdata_train))
    
    for (i in 1:nrow(mdata_train)) {
      
      ap[i] <- if (i %% Tt == 1 || sum(v[((i-1) %/% Tt *Tt +1):(i-1)])==0) {a10}
      else if (v[i-1]==0 && skip==TRUE) {a[i-1]} 
      else {q[i-1]    /delta   *a[i-1]}
      bp[i] <- if (i %% Tt == 1 || sum(v[((i-1) %/% Tt *Tt +1):(i-1)])==0) {a10}  
      else if (v[i-1]==0 && skip==TRUE) {b[i-1]} 
      else {q[i-1]*((1/delta-1)*a[i-1]+b[i-1])}
      # for each policyholder, ap[i] and bp[i] are reset at the first period
      # (assuming that the dataframe is sorted by ID and HY)
      
      a[i] <- ap[i] + v[i]/psi 
      b[i] <- bp[i] + Ys[i]/psi
      
      q[i] <- delta*a10 / (a[i] - delta^2*a[i] + delta^2*a10) 
      
    }
    neglik <- -sum(dGP(Y, v/psi, ap, bp*mu*psi, log=TRUE)[v>0])
    return(list(a=a, b=b, q=q, ap=ap, bp=bp, neglik=neglik)) }
  
  homo_BM_neglik <- function(parm) { # only 2 parameters, delta = 1 by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- 1
    result <- ggSM(a10, psi, delta, heterogenous = FALSE)$neglik
    return(result) }
  
  homo_ScM_neglik <- function(parm) { # 3 parameters, delta in [0,1] by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- parm[3]
    result <- ggSM(a10, psi, delta, heterogenous = FALSE, skip = FALSE)$neglik
    return(result) }
  
  hetero_BM_neglik <- function(parm) { # only 2 parameters, delta = 1 by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- 1
    result <- ggSM(a10, psi, delta, heterogenous = TRUE)$neglik
    return(result) }
  
  hetero_ScM_neglik <- function(parm) { # 3 parameters, delta in [0,1] by definition
    a10   <- parm[1]
    psi   <- parm[2]
    delta <- parm[3]
    result <- ggSM(a10, psi, delta, heterogenous = TRUE, skip = FALSE)$neglik
    return(result) }
  
  
  homo_BM_fit <- optim(c(1, 1), homo_BM_neglik, method="L-BFGS-B",
                       lower = c(0,0), upper = c(Inf, Inf))
  
  homo_ScM_fit <- optim(c(homo_BM_fit$par/2, 1), homo_ScM_neglik, method="L-BFGS-B",
                        lower = c(0,0,0), upper = c(Inf, Inf,1))
  
  hetero_BM_fit <- optim(c(1,1), hetero_BM_neglik, method="L-BFGS-B",
                         lower = c(0,0), upper = c(Inf, Inf))
  
  hetero_ScM_fit <- optim(c(hetero_BM_fit$par/2, 1), hetero_ScM_neglik, method="L-BFGS-B",
                          lower = c(0,0,0), upper = c(Inf, Inf, 1))
  
  homo_BM_fit$value
  homo_ScM_fit$value
  hetero_BM_fit$value
  hetero_ScM_fit$value
  
  endindex <- (1:nrow(mdata_train)) %% Tt == 0
  
  # Crediblity factors for t+1:   b_{t+1|t}/a_{t+1|t}             (from (2.4))
  #                            = (b_{t}/a_{t}-1) * Delta_t + 1    (from (3.6))
  
  homo_BM_cred <-(ggSM(a10=homo_BM_fit$par[1], psi=homo_BM_fit$par[2], 
                       delta=1, heterogenous=FALSE)$b[endindex]/
                    ggSM(a10=homo_BM_fit$par[1], psi=homo_BM_fit$par[2], 
                         delta=1, heterogenous=FALSE)$a[endindex])
  
  hetero_BM_cred <-(ggSM(a10=hetero_BM_fit$par[1], psi=hetero_BM_fit$par[2], 
                         delta=1, heterogenous=TRUE)$b[endindex]/
                      ggSM(a10=hetero_BM_fit$par[1], psi=hetero_BM_fit$par[2], 
                           delta=1, heterogenous=TRUE)$a[endindex])
  
  homo_ScM_cred <-(ggSM(a10=homo_ScM_fit$par[1], psi=homo_ScM_fit$par[2], 
                        delta=homo_ScM_fit$par[3], heterogenous=FALSE, skip = FALSE)$b[endindex]/
                     ggSM(a10=homo_ScM_fit$par[1], psi=homo_ScM_fit$par[2], 
                          delta=homo_ScM_fit$par[3], heterogenous=FALSE, skip = FALSE)$a[endindex]-1)*
    homo_ScM_fit$par[3] + 1
  
  hetero_ScM_cred <-(ggSM(a10=hetero_ScM_fit$par[1], psi=hetero_ScM_fit$par[2], 
                          delta=hetero_ScM_fit$par[3], heterogenous=TRUE, skip = FALSE)$b[endindex]/
                       ggSM(a10=hetero_ScM_fit$par[1], psi=hetero_ScM_fit$par[2], 
                            delta=hetero_ScM_fit$par[3], heterogenous=TRUE, skip = FALSE)$a[endindex]-1)*
    hetero_ScM_fit$par[3] + 1
  
  dzd <- aggregate(mdata_train[,c(3,5,6)], by=list(mdata_train$ID), FUN = sum)
  
  homo_CM_cred <- (dzd$Y/hmmean + homo_BM_fit$par[1]*homo_BM_fit$par[2]) /
                  (dzd$nu       + homo_BM_fit$par[1]*homo_BM_fit$par[2])
    
  hetero_CM_cred <- (dzd$Ypermu + homo_BM_fit$par[1]*homo_BM_fit$par[2]) /
                    (dzd$nu     + homo_BM_fit$par[1]*homo_BM_fit$par[2])
  
  
  zzzz <- data.frame(mdata_test$ID, homo_BM_cred, hetero_BM_cred, 
                     homo_ScM_cred, hetero_ScM_cred)
  colnames(zzzz)[1] <- "ID"
  
  zdz <- merge(x = zzzz, y = mdata_test,
               by = "ID", all.x = TRUE)
  
  
  valtable_RMSE <- 
    rbind(
      c(RMSE(zdz$Y, hmmean                       * zdz$nu),
        RMSE(zdz$Y, hmmean * zdz$homo_BM_cred    * zdz$nu),
        RMSE(zdz$Y, hmmean * zdz$homo_ScM_cred   * zdz$nu)),
      c(RMSE(zdz$Y, zdz$mu                       * zdz$nu),
        RMSE(zdz$Y, zdz$mu * zdz$hetero_BM_cred  * zdz$nu),
        RMSE(zdz$Y, zdz$mu * zdz$hetero_ScM_cred * zdz$nu)))
  
  colnames(valtable_RMSE) <- c("Indep", "Buhlmann","ScM")
  rownames(valtable_RMSE) <- c("Homogenous", "Heterogeneous")
  
  
  
  GDEV <- function(y, mu, nu) {
    vec    <- -nu*log(y/mu/nu) + (y-mu*nu)/mu
    result <- 2*sum( vec[nu>0] ) }
  
  valtable_GDEV <- 
    rbind(
      c(GDEV(zdz$Y, hmmean                      , zdz$nu),
        GDEV(zdz$Y, hmmean * zdz$homo_BM_cred   , zdz$nu),
        GDEV(zdz$Y, hmmean * zdz$homo_ScM_cred  , zdz$nu)),
      c(GDEV(zdz$Y, zdz$mu                      , zdz$nu),
        GDEV(zdz$Y, zdz$mu * zdz$hetero_BM_cred , zdz$nu),
        GDEV(zdz$Y, zdz$mu * zdz$hetero_ScM_cred, zdz$nu)))
  
  colnames(valtable_GDEV) <- c("Indep", "Buhlmann", "ScM")
  rownames(valtable_GDEV) <- c("Homogenous", "Heterogeneous")
  
  valtable_RMSE
  valtable_GDEV
  
  parmtable_dep <- round(rbind(c(
    homo_BM_fit$par , 1),
    homo_ScM_fit$par,
    c( hetero_BM_fit$par , 1),
    hetero_ScM_fit$par), 4 )
  
  colnames(parmtable_dep) <- c("$a_{10}$", "$\\psi$", "$\\Delta$")
  rownames(parmtable_dep) <- c("Homo_BM"  , "Homo_Scm"  ,
                               "Hetero_BM", "Hetero_Scm")
  
  print(r)
  tablist_delta1.0$valtabRMSEs[r,,] <- valtable_RMSE
  tablist_delta1.0$valtabGDEVs[r,,] <- valtable_GDEV
  tablist_delta1.0$parmtabs[   r,,] <- parmtable_dep
}



#### summary of the tables ####

delta0.5_RMSE_mean <- rbind(
  c(mean(tablist_delta0.5$valtabRMSEs[,1,1]),
    mean(tablist_delta0.5$valtabRMSEs[,1,2]),
    mean(tablist_delta0.5$valtabRMSEs[,1,3])),
  c(mean(tablist_delta0.5$valtabRMSEs[,2,1]),
    mean(tablist_delta0.5$valtabRMSEs[,2,2]),
    mean(tablist_delta0.5$valtabRMSEs[,2,3])))

delta0.5_RMSE_sd <- rbind(
  c(sd(tablist_delta0.5$valtabRMSEs[,1,1]),
    sd(tablist_delta0.5$valtabRMSEs[,1,2]),
    sd(tablist_delta0.5$valtabRMSEs[,1,3])),
  c(sd(tablist_delta0.5$valtabRMSEs[,2,1]),
    sd(tablist_delta0.5$valtabRMSEs[,2,2]),
    sd(tablist_delta0.5$valtabRMSEs[,2,3])))

delta0.5_GDEV_mean <- rbind(
  c(mean(tablist_delta0.5$valtabGDEVs[,1,1]),
    mean(tablist_delta0.5$valtabGDEVs[,1,2]),
    mean(tablist_delta0.5$valtabGDEVs[,1,3])),
  c(mean(tablist_delta0.5$valtabGDEVs[,2,1]),
    mean(tablist_delta0.5$valtabGDEVs[,2,2]),
    mean(tablist_delta0.5$valtabGDEVs[,2,3])))

delta0.5_GDEV_sd <- rbind(
  c(sd(tablist_delta0.5$valtabGDEVs[,1,1]),
    sd(tablist_delta0.5$valtabGDEVs[,1,2]),
    sd(tablist_delta0.5$valtabGDEVs[,1,3])),
  c(sd(tablist_delta0.5$valtabGDEVs[,2,1]),
    sd(tablist_delta0.5$valtabGDEVs[,2,2]),
    sd(tablist_delta0.5$valtabGDEVs[,2,3])))


delta0.5simvaltab <- 
cbind(
rbind(      format(round(delta0.5_RMSE_mean[c(1,3,5)], 2),nsmall=2),
paste0("(", format(round(delta0.5_RMSE_sd[  c(1,3,5)] ,2),nsmall=2), ")"),
            format(round(delta0.5_RMSE_mean[c(2,4,6)], 2),nsmall=2),
paste0("(", format(round(delta0.5_RMSE_sd[  c(2,4,6)] ,2),nsmall=2), ")")),
rbind(      format(round(delta0.5_GDEV_mean[c(1,3,5)], 2),nsmall=2),
paste0("(", format(round(delta0.5_GDEV_sd[  c(1,3,5)] ,2),nsmall=2), ")"),
            format(round(delta0.5_GDEV_mean[c(2,4,6)], 2),nsmall=2),
paste0("(", format(round(delta0.5_GDEV_sd[  c(2,4,6)] ,2),nsmall=2), ")")))


rownames(delta0.5simvaltab) <- c("Homogenous","","Heterogenous","")
colnames(delta0.5simvaltab) <- c("Independent", "Buhlmann", "SSM",
                                 "Independent", "Buhlmann", "SSM")


delta1.0_RMSE_mean <- rbind(
  c(mean(tablist_delta1.0$valtabRMSEs[,1,1]),
    mean(tablist_delta1.0$valtabRMSEs[,1,2]),
    mean(tablist_delta1.0$valtabRMSEs[,1,3])),
  c(mean(tablist_delta1.0$valtabRMSEs[,2,1]),
    mean(tablist_delta1.0$valtabRMSEs[,2,2]),
    mean(tablist_delta1.0$valtabRMSEs[,2,3])))

delta1.0_RMSE_sd <- rbind(
  c(sd(tablist_delta1.0$valtabRMSEs[,1,1]),
    sd(tablist_delta1.0$valtabRMSEs[,1,2]),
    sd(tablist_delta1.0$valtabRMSEs[,1,3])),
  c(sd(tablist_delta1.0$valtabRMSEs[,2,1]),
    sd(tablist_delta1.0$valtabRMSEs[,2,2]),
    sd(tablist_delta1.0$valtabRMSEs[,2,3])))

delta1.0_GDEV_mean <- rbind(
  c(mean(tablist_delta1.0$valtabGDEVs[,1,1]),
    mean(tablist_delta1.0$valtabGDEVs[,1,2]),
    mean(tablist_delta1.0$valtabGDEVs[,1,3])),
  c(mean(tablist_delta1.0$valtabGDEVs[,2,1]),
    mean(tablist_delta1.0$valtabGDEVs[,2,2]),
    mean(tablist_delta1.0$valtabGDEVs[,2,3])))

delta1.0_GDEV_sd <- rbind(
  c(sd(tablist_delta1.0$valtabGDEVs[,1,1]),
    sd(tablist_delta1.0$valtabGDEVs[,1,2]),
    sd(tablist_delta1.0$valtabGDEVs[,1,3])),
  c(sd(tablist_delta1.0$valtabGDEVs[,2,1]),
    sd(tablist_delta1.0$valtabGDEVs[,2,2]),
    sd(tablist_delta1.0$valtabGDEVs[,2,3])))

delta1.0simvaltab <- 
cbind(
rbind(      format(round(delta1.0_RMSE_mean[c(1,3,5)], 2),nsmall=2),
paste0("(", format(round(delta1.0_RMSE_sd[  c(1,3,5)] ,2),nsmall=2), ")"),
            format(round(delta1.0_RMSE_mean[c(2,4,6)], 2),nsmall=2),
paste0("(", format(round(delta1.0_RMSE_sd[  c(2,4,6)] ,2),nsmall=2), ")")),
rbind(      format(round(delta1.0_GDEV_mean[c(1,3,5)], 2),nsmall=2),
paste0("(", format(round(delta1.0_GDEV_sd[  c(1,3,5)] ,2),nsmall=2), ")"),
            format(round(delta1.0_GDEV_mean[c(2,4,6)], 2),nsmall=2),
paste0("(", format(round(delta1.0_GDEV_sd[  c(2,4,6)] ,2),nsmall=2), ")")))

rownames(delta1.0simvaltab) <- c("Homogeneous","","Heterogeneous","")
colnames(delta1.0simvaltab) <- c("Independent", "Buhlmann", "SSM",
                                 "Independent", "Buhlmann", "SSM")

sim_esttab <- cbind(
  rbind(format(round(rbind(
    c(mean(tablist_delta0.5$parmtabs[,1,1]),
      mean(tablist_delta0.5$parmtabs[,1,2]),
      mean(tablist_delta0.5$parmtabs[,1,3])),
    c(mean(tablist_delta0.5$parmtabs[,2,1]),
      mean(tablist_delta0.5$parmtabs[,2,2]),
      mean(tablist_delta0.5$parmtabs[,2,3])),
    c(mean(tablist_delta0.5$parmtabs[,3,1]),
      mean(tablist_delta0.5$parmtabs[,3,2]),
      mean(tablist_delta0.5$parmtabs[,3,3])),
    c(mean(tablist_delta0.5$parmtabs[,4,1]),
      mean(tablist_delta0.5$parmtabs[,4,2]),
      mean(tablist_delta0.5$parmtabs[,4,3]))),4),nsmall=4),
    matrix(ncol=3, 
           paste0("(", format(round(rbind(
             c(sd(tablist_delta0.5$parmtabs[,1,1]),
               sd(tablist_delta0.5$parmtabs[,1,2]),
               sd(tablist_delta0.5$parmtabs[,1,3])),
             c(sd(tablist_delta0.5$parmtabs[,2,1]),
               sd(tablist_delta0.5$parmtabs[,2,2]),
               sd(tablist_delta0.5$parmtabs[,2,3])),
             c(sd(tablist_delta0.5$parmtabs[,3,1]),
               sd(tablist_delta0.5$parmtabs[,3,2]),
               sd(tablist_delta0.5$parmtabs[,3,3])),
             c(sd(tablist_delta0.5$parmtabs[,4,1]),
               sd(tablist_delta0.5$parmtabs[,4,2]),
               sd(tablist_delta0.5$parmtabs[,4,3]))), 4),nsmall=4), ")")))[c(1,5,2,6, 3,7,4,8), ],
  rbind(format(round(rbind(
    c(mean(tablist_delta1.0$parmtabs[,1,1]),
      mean(tablist_delta1.0$parmtabs[,1,2]),
      mean(tablist_delta1.0$parmtabs[,1,3])),
    c(mean(tablist_delta1.0$parmtabs[,2,1]),
      mean(tablist_delta1.0$parmtabs[,2,2]),
      mean(tablist_delta1.0$parmtabs[,2,3])),
    c(mean(tablist_delta1.0$parmtabs[,3,1]),
      mean(tablist_delta1.0$parmtabs[,3,2]),
      mean(tablist_delta1.0$parmtabs[,3,3])),
    c(mean(tablist_delta1.0$parmtabs[,4,1]),
      mean(tablist_delta1.0$parmtabs[,4,2]),
      mean(tablist_delta1.0$parmtabs[,4,3]))),4),nsmall=4),
    matrix(ncol=3, 
           paste0("(", format(round(rbind(
             c(sd(tablist_delta1.0$parmtabs[,1,1]),
               sd(tablist_delta1.0$parmtabs[,1,2]),
               sd(tablist_delta1.0$parmtabs[,1,3])),
             c(sd(tablist_delta1.0$parmtabs[,2,1]),
               sd(tablist_delta1.0$parmtabs[,2,2]),
               sd(tablist_delta1.0$parmtabs[,2,3])),
             c(sd(tablist_delta1.0$parmtabs[,3,1]),
               sd(tablist_delta1.0$parmtabs[,3,2]),
               sd(tablist_delta1.0$parmtabs[,3,3])),
             c(sd(tablist_delta1.0$parmtabs[,4,1]),
               sd(tablist_delta1.0$parmtabs[,4,2]),
               sd(tablist_delta1.0$parmtabs[,4,3]))), 4),nsmall=4), ")")))[c(1,5,2,6, 3,7,4,8), ])
sim_esttab[2,3] <- "-"
sim_esttab[2,6] <- "-"
sim_esttab[6,3] <- "-"
sim_esttab[6,6] <- "-"
rownames(sim_esttab) <- c("Homogeneous", "Buhlmann", "Homogeneous","SSM",
                          "Heterogeneous","Buhlmann","Heterogeneous","SSM")
colnames(sim_esttab) <- c("$\\hat{a}_{1|0}$", "$\\hat{\\psi}$", "$\\hat{\\Delta}$",
                          "$\\hat{a}_{1|0}$", "$\\hat{\\psi}$", "$\\hat{\\Delta}$")

library(knitr)
library(kableExtra)
options(knitr.table.format = "latex")
options(knitr.kable.NA = '')

kable(delta0.5simvaltab, digits=2, booktabs = T, align='c',
      caption = "Summary of out-of-sample validation with simulations ($\\Delta=0.5$)",
      linesep = c("", "", "", "",  "\\hline"),
      escape=FALSE)  %>%
  add_header_above(c(" "=1, "RMSE" = 3,"GDEV" = 3))



kable(delta1.0simvaltab, digits=2, booktabs = T, align='c',
      caption = "Summary of out-of-sample validation with simulations ($\\Delta=1.0$)",
      linesep = c("", "", "", "",  "\\hline"),
      escape=FALSE)  %>%
  add_header_above(c(" "=1, "RMSE" = 3,"GDEV" = 3))


kable(sim_esttab, digits=2, booktabs = T, align='c',
      caption = "Summary of estimation with simulations where $a_{1|0}=3, \\psi=1$",
      linesep = c("", "", "", "", "", "","", "", "\\hline"),
      escape=FALSE)  %>%
  add_header_above(c(" "=1, "$\\Delta=0.5$" = 3,"$\\Delta=1.0$" = 3), escape=FALSE)

