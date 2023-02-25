# library(gss)

################################################################################
################################### gss ########################################
################################################################################
### fit gss model with individual effect and treatment effect 

# source("./Biofluc_gss_dataLoader.R")

miRNAs = colnames(CountsWithGroupInfo.normalized)[12:dim(CountsWithGroupInfo.normalized)[2]]

pval.func = function(miRNA){
  warn <- err <- NULL
  res <- withCallingHandlers(
    tryCatch({
      count = get(miRNA, CountsWithGroupInfo.normalized)
      fit1 = gssanova(count ~ Time + Subject_Group + InfectionStatus + InfectionStatus:Subject_Group + Run,
                      data = CountsWithGroupInfo.normalized, seed = 3,
                      family = "poisson", random = ~1|Subject_ID) 
      fit0.1 = gssanova(count ~ Time + Subject_Group + Run,
                        data = CountsWithGroupInfo.normalized, seed = 3,
                        family = "poisson", random = ~1|Subject_ID)
      fit0.2 = gssanova(count ~ Time + InfectionStatus + Run,
                        data = CountsWithGroupInfo.normalized, seed = 3,
                        family = "poisson", random = ~1|Subject_ID)
      fit0.3 = gssanova(count ~ Time + InfectionStatus + InfectionStatus:Subject_Group + Run,
                        data = CountsWithGroupInfo.normalized, seed = 3,
                        family = "poisson", random = ~1|Subject_ID)
      fit0.4 = gssanova(count ~ Time + Subject_Group + InfectionStatus:Subject_Group + Run,
                        data = CountsWithGroupInfo.normalized, seed = 3,
                        family = "poisson", random = ~1|Subject_ID) 
      sum.exp.eta1 = 0
      sum.exp.eta0.1 = sum.exp.eta0.2 = sum.exp.eta0.3 = sum.exp.eta0.4 = 0
      for (j in 1:length(count)){
        sum.exp.eta1 = sum.exp.eta1 + exp(fit1$eta[j])
        sum.exp.eta0.1 = sum.exp.eta0.1 + exp(fit0.1$eta[j])
        sum.exp.eta0.2 = sum.exp.eta0.2 + exp(fit0.2$eta[j])
        sum.exp.eta0.3 = sum.exp.eta0.3 + exp(fit0.3$eta[j])
        sum.exp.eta0.4 = sum.exp.eta0.4 + exp(fit0.4$eta[j])
      }
      MLE.1 = 2*(-sum(fit1$eta * count) + sum.exp.eta1)
      MLE.0.1 = 2*(-sum(fit0.1$eta * count) + sum.exp.eta0.1)
      MLE.0.2 = 2*(-sum(fit0.2$eta * count) + sum.exp.eta0.2)
      MLE.0.3 = 2*(-sum(fit0.3$eta * count) + sum.exp.eta0.3)
      MLE.0.4 = 2*(-sum(fit0.4$eta * count) + sum.exp.eta0.4)
      
      delta2LL.1 = MLE.0.1 - MLE.1  
      delta2LL.2 = MLE.0.2 - MLE.1  
      delta2LL.3 = MLE.0.3 - MLE.1  
      delta2LL.4 = MLE.0.4 - MLE.1  
      
      pval.1 = pchisq(delta2LL.1, df = 2, lower.tail = FALSE)
      pval.2 = pchisq(delta2LL.2, df = 2, lower.tail = FALSE)
      pval.3 = pchisq(delta2LL.3, df = 1, lower.tail = FALSE)
      pval.4 = pchisq(delta2LL.4, df = 1, lower.tail = FALSE)
      
    },
    error=function(e){err <<- conditionMessage(e)
    NULL
    }), warning=function(w) {
      warn <<- append(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  return(list(pval.1, pval.2, pval.3, pval.4, err, warn))
}


