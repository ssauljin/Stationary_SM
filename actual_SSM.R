library(caret)
library(dplyr)

library(knitr)
library(kableExtra)
options(knitr.table.format = "latex")
options(knitr.kable.NA = '')

rdata <- read.csv("MEPS_marginal.csv", header=TRUE,
                          stringsAsFactors=TRUE)
rdata$DUPERSID <- paste0(rdata$DUPERSID,rdata$OPDATEYR,"H",ceiling(rdata$OPDATEMM/6))

head(rdata)
str(rdata)

round(summary(rdata$Doctor_Type)/nrow(rdata)*100,2)
round(summary(rdata$Care_Category)/nrow(rdata)*100,2)
round(summary(rdata$Special_Cond)/nrow(rdata)*100,2)
round(summary(rdata$Surgery)/nrow(rdata)*100,2)
round(summary(rdata$Prescription)/nrow(rdata)*100,2)
round(summary(rdata$Telehealth)/nrow(rdata)*100,2)

gamma_fit <- glm(Total_Charge ~ Doctor_Type+ Care_Category + Special_Cond +
                   Surgery + Prescription + Telehealth,
                 data=rdata, family=Gamma(link="log"))


marginal_esttable <- summary(gamma_fit)$coef[,c(1,4)]
colnames(marginal_esttable) <- c("Estimate", "p-value")

kable(marginal_esttable, digits=4, booktabs = T, align='c',
      caption = "Summary of the estimated regression coefficients for the working severity model \\label{tab:marginal_est}",
      linesep = c("", "", "", "", "", "", "", "", "", "", "", "",  "\\hline"),
      escape=FALSE) 



rdata <- data.frame(gamma_fit$data, exp(predict(gamma_fit)))
colnames(rdata)[ncol(rdata)] <- "mu"
rdata$nu    <- 1

zz <- aggregate(rdata[,12:14], by = list(rdata$DUPERSID), FUN = sum)
colnames(zz)[1] <- "IDnHY"

uIDs <- unique(substr(zz$IDnHY, 1, 10))

zzz <- data.frame(paste0(rep(uIDs, each=8),
rep(paste0(rep(2019:2022, each=2), c("H1", "H2")), times = length(uIDs))))

colnames(zzz)[1] <- "IDnHY"

mdata <- merge(x=zzz, y=zz, 
               by=c("IDnHY"), all.x=TRUE)
colnames(mdata)[2] <- "Y"
mdata$mu <- mdata$mu / mdata$nu
mdata$Ypermu <- mdata$Y / mdata$mu
mdata <- mdata %>% replace(is.na(.), 0)

rm(zz, zzz)

mdata$ID <- substr(mdata$IDnHY,1,10)
mdata$HY <- substr(mdata$IDnHY,11,16)

mdata_train <- mdata[mdata$HY != "2022H2",]
mdata_test  <- mdata[mdata$HY == "2022H2",]

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

homo_SkM_neglik <- function(parm) { # 3 parameters, delta in [0,1] by definition
  a10   <- parm[1]
  psi   <- parm[2]
  delta <- parm[3]
  result <- ggSM(a10, psi, delta, heterogenous = FALSE, skip = TRUE)$neglik
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

hetero_SkM_neglik <- function(parm) { # 3 parameters, delta in [0,1] by definition
  a10   <- parm[1]
  psi   <- parm[2]
  delta <- parm[3]
  result <- ggSM(a10, psi, delta, heterogenous = TRUE, skip = TRUE)$neglik
  return(result) }

homo_BM_fit <- optim(c(3, 1/3), homo_BM_neglik, method="L-BFGS-B",
                     lower = c(0,0), upper = c(Inf, Inf))

homo_ScM_fit <- optim(c(homo_BM_fit$par, 0.99), homo_ScM_neglik, method="L-BFGS-B",
                     lower = c(0,0,0), upper = c(Inf, Inf,1))

homo_SkM_fit <- optim(c(homo_BM_fit$par, 0.99), homo_SkM_neglik, method="L-BFGS-B",
                      lower = c(0,0,0), upper = c(Inf, Inf,1))

hetero_BM_fit <- optim(c(3,1/3), hetero_BM_neglik, method="L-BFGS-B",
                     lower = c(0,0), upper = c(Inf, Inf))

hetero_ScM_fit <- optim(c(hetero_BM_fit$par, 0.99), hetero_ScM_neglik, method="L-BFGS-B",
                       lower = c(0,0,0), upper = c(Inf, Inf, 1))

hetero_SkM_fit <- optim(c(hetero_BM_fit$par, 0.99), hetero_SkM_neglik, method="L-BFGS-B",
                        lower = c(0,0,0), upper = c(Inf, Inf, 1))


homo_BM_fit$value
homo_ScM_fit$value
homo_SkM_fit$value
hetero_BM_fit$value
hetero_ScM_fit$value
hetero_SkM_fit$value

endindex <- (1:nrow(mdata_train)) %% 7 == 0

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

homo_SkM_cred <-(ggSM(a10=homo_SkM_fit$par[1], psi=homo_SkM_fit$par[2], 
                      delta=homo_SkM_fit$par[3], heterogenous=FALSE, skip = TRUE)$b[endindex]/
                 ggSM(a10=homo_SkM_fit$par[1], psi=homo_SkM_fit$par[2], 
                      delta=homo_SkM_fit$par[3], heterogenous=FALSE, skip = TRUE)$a[endindex]-1)*
                            homo_SkM_fit$par[3] + 1  

hetero_SkM_cred <-(ggSM(a10=hetero_SkM_fit$par[1], psi=hetero_SkM_fit$par[2], 
                      delta=hetero_SkM_fit$par[3], heterogenous=TRUE, skip = TRUE)$b[endindex]/
                   ggSM(a10=hetero_SkM_fit$par[1], psi=hetero_SkM_fit$par[2], 
                      delta=hetero_SkM_fit$par[3], heterogenous=TRUE, skip = TRUE)$a[endindex]-1)*
                            hetero_SkM_fit$par[3] + 1  

zzzz <- data.frame(mdata_test$ID, homo_BM_cred, hetero_BM_cred,
                   homo_ScM_cred, hetero_ScM_cred, homo_SkM_cred, hetero_SkM_cred)
colnames(zzzz)[1] <- "ID"

zdz <- merge(x = zzzz, y = mdata_test,
             by = "ID", all.x = TRUE)


valtable_RMSE <- 
  rbind(
    c(RMSE(zdz$Y, hmmean                       * zdz$nu),
      RMSE(zdz$Y, hmmean * zdz$homo_BM_cred    * zdz$nu),
      RMSE(zdz$Y, hmmean * zdz$homo_ScM_cred   * zdz$nu),
      RMSE(zdz$Y, hmmean * zdz$homo_SkM_cred   * zdz$nu)),
    c(RMSE(zdz$Y, zdz$mu                       * zdz$nu),
      RMSE(zdz$Y, zdz$mu * zdz$hetero_BM_cred  * zdz$nu),
      RMSE(zdz$Y, zdz$mu * zdz$hetero_ScM_cred * zdz$nu),
      RMSE(zdz$Y, zdz$mu * zdz$hetero_SkM_cred * zdz$nu)))

colnames(valtable_RMSE) <- c("Indep", "Buhlmann", "SM (w/o skip)", "SM (with skip)")
rownames(valtable_RMSE) <- c("Homogenous", "Heterogeneous")



GDEV <- function(y, mu, nu) {
  vec    <- -nu*log(y/mu/nu) + (y-mu*nu)/mu
  result <- 2*sum( vec[nu>0] ) }

valtable_GDEV <- 
  rbind(
    c(GDEV(zdz$Y, hmmean                      , zdz$nu),
      GDEV(zdz$Y, hmmean * zdz$homo_BM_cred   , zdz$nu),
      GDEV(zdz$Y, hmmean * zdz$homo_ScM_cred  , zdz$nu),
      GDEV(zdz$Y, hmmean * zdz$homo_SkM_cred  , zdz$nu)),
    c(GDEV(zdz$Y, zdz$mu                      , zdz$nu),
      GDEV(zdz$Y, zdz$mu * zdz$hetero_BM_cred , zdz$nu),
      GDEV(zdz$Y, zdz$mu * zdz$hetero_ScM_cred, zdz$nu),
      GDEV(zdz$Y, zdz$mu * zdz$hetero_SkM_cred, zdz$nu)))

colnames(valtable_GDEV) <- c("Indep", "Buhlmann", "SM (w/o skip)", "SM (with skip)")
rownames(valtable_GDEV) <- c("Homogenous", "Heterogeneous")

valtable_RMSE
valtable_GDEV

act_esttab <- round(rbind(c(
     homo_BM_fit$par , 1),
    homo_ScM_fit$par ,
    homo_SkM_fit$par ,
c( hetero_BM_fit$par , 1),
  hetero_ScM_fit$par ,
  hetero_SkM_fit$par), 4 )

colnames(act_esttab) <- c("$\\hat{a}_{1|0}$", "$\\hat{\\psi}$", "$\\hat{\\Delta}$")
rownames(act_esttab) <- c("Homogeneous Buhlmann"  , "Homogeneous SSM"  , "Homo_Skm" ,
                             "Heterogeneous Buhlmann", "Heterogeneous SSM", "Hetero_Skm")

actvaltab <- cbind(valtable_RMSE[,1:3], valtable_GDEV[,1:3])


kable(actvaltab, digits=2, booktabs = T, align='c',
      caption = "Summary of out-of-sample validation with the actual dataset",
      linesep = c("", "", "", "",  "\\hline"),
      escape=FALSE)  %>%
  add_header_above(c(" "=1, "RMSE" = 3,"GDEV" = 3))

kable(act_esttab[c(1:2,4:5),], digits=4, booktabs = T, align='c',
      caption = "Summary of estimation with the actual dataset \ref{tab:actest}",
      linesep = c("", "", "", "", "", "","", "", "\\hline"),
      escape=FALSE)  

