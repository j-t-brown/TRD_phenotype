### TRD Rotation Project ###
# optimism adjustment

# load bolasso coefficients
coeffs <- coeffs<-read.csv('',header=T)
coeffs <- coeffs[2:length(coeffs)] # remove X

# load dataset
dataset<-read.csv('',header=T)

# data preprocessing
dataset$DEMO__ADI<-as.numeric(dataset$DEMO__ADI)
colnames(dataset)[colnames(dataset)=="OUTCOME"]<-"outcome"
dataset$DEMO__GENDER__M<-ifelse(dataset$DEMO__GENDER=="M",1,0)
dataset$DEMO__GENDER__F<-ifelse(dataset$DEMO__GENDER=="F",1,0)
dataset$DEMO__GENDER__U<-ifelse(dataset$DEMO__GENDER=="U",1,0)
dataset$DEMO__RACE__A<-ifelse(dataset$DEMO__RACE=="A",1,0)
dataset$DEMO__RACE__B<-ifelse(dataset$DEMO__RACE=="B",1,0)
dataset$DEMO__RACE__H<-ifelse(dataset$DEMO__RACE=="H",1,0)
dataset$DEMO__RACE__I<-ifelse(dataset$DEMO__RACE=="I",1,0)
dataset$DEMO__RACE__O<-ifelse(dataset$DEMO__RACE=="O",1,0)
dataset$DEMO__RACE__U<-ifelse(dataset$DEMO__RACE=="U",1,0)
dataset$DEMO__RACE__W<-ifelse(dataset$DEMO__RACE=="W",1,0)

# transform ICD and RXNORM counts
log1p_feat=names(dataset)[grep("RXNORM|CCS",names(dataset))]
dataset[log1p_feat]<-log1p(dataset[log1p_feat])

# helper functions
extract_features <- function(dataframe,threshold){
  # returns a vector of features from dataframe that have equal to or
  # greater than threshold percentage of non-zero coefficients
  df <- data.frame(dataframe)
  df[is.na(df)] <- 0
  df[df!=0] <- 1
  
  features <- names(df)
  all_features <- c()
  for (feat in features){
    if (mean(df[,feat]) >= threshold){
      all_features <- c(all_features,feat)
    }
  }
  return(all_features) 
}

top_decile_pct <- function(preds,outcomes){
  # calculates the percentage of all cases
  # in the top decile of all predictions
  q<-quantile(preds, 0.9)
  total_cases <- sum(outcomes)
  n_cases <-sum(outcomes[which(preds>=q)])
  return(n_cases/total_cases)
}

roc_diff <- function(model,matrix=F){
  # returns the difference between the auc of the roc
  # between the model being applied to the bootstrapped dataset
  # and the original dataset
  if (matrix==T){
    p1 <- predict(model,data.matrix(data.boot[glm_features_100]),
                  type="response",standardize=TRUE)
    p2 <- predict(model,data.matrix(dataset[glm_features_100]),
                  type="response",standardize=TRUE)
  } else {
    p1 <- predict(model,data.boot[glm_features_100],
                  type="response",standardize=TRUE)
    p2 <- predict(model,dataset[glm_features_100],
                  type="response",standardize=TRUE)
  }
  ans <- as.numeric(unlist(calculate_auc(p1,data.boot$outcome))@y.values) -
    as.numeric(unlist(calculate_auc(p2,dataset$outcome))@y.values)
  return(ans)
}

pr_diff <- function(model,matrix=F){
  # returns the difference between the auc of the pr curve
  # between the model being applied to the bootstrapped dataset
  # and the original dataset
  if (matrix==T){
    p1 <- predict(model,data.matrix(data.boot[glm_features_100]),
                  type="response",standardize=TRUE)
    p2 <- predict(model,data.matrix(dataset[glm_features_100]),
                  type="response",standardize=TRUE)
  } else {
  p1 <- predict(model,data.boot[glm_features_100],
                type="response",standardize=TRUE)
  p2 <- predict(model,dataset[glm_features_100],
                type="response",standardize=TRUE)
  }
  ans <- PRAUC(p1,data.boot$outcome) -
    PRAUC(p2,dataset$outcome)
  return(ans)
}

dec_diff <- function(model,matrix=F){
  # returns the difference between the top decile pct
  # between the model being applied to the bootstrapped dataset
  # and the original dataset
  if (matrix==T){
    p1 <- predict(model,data.matrix(data.boot[glm_features_100]),
                  type="response",standardize=TRUE)
    p2 <- predict(model,data.matrix(dataset[glm_features_100]),
                  type="response",standardize=TRUE)
  } else {
    p1 <- predict(model,data.boot[glm_features_100],
                  type="response",standardize=TRUE)
    p2 <- predict(model,dataset[glm_features_100],
                  type="response",standardize=TRUE)
  }
  ans <- top_decile_pct(p1,data.boot$outcome) -
    top_decile_pct(p2,dataset$outcome)
  return(ans)
}

cal_diff <- function(model,matrix=F){
  # returns the difference between the Dxy of the calibration curves
  # between the model being applied to the bootstrapped dataset
  # and the original dataset
  if (matrix==T){
    p1 <- predict(model,data.matrix(data.boot[glm_features_100]),
                  type="response",standardize=TRUE)
    p2 <- predict(model,data.matrix(dataset[glm_features_100]),
                  type="response",standardize=TRUE)
  } else {
    p1 <- predict(model,data.boot[glm_features_100],
                  type="response",standardize=TRUE)
    p2 <- predict(model,dataset[glm_features_100],
                  type="response",standardize=TRUE)
  }
  ans <- as.numeric(val.prob(p1,data.boot$outcome,pl=F)[1])-
    as.numeric(val.prob(p2,dataset$outcome,pl=F)[1])
  return(ans)
}

# boLASSO selected features for models
glm_features_100 <- extract_features(coeffs,1.0) #significant features in 100% of boLASSO iterations

# generate base models on full dataset
# glm model w/o interactions
glm_100.base<-glm(outcome~.,dataset[c("outcome",glm_features_100)],family=binomial,model=FALSE)

# include interaction terms
glm_int.base<-glm(outcome~.^2,dataset[c("outcome",glm_features_100)],family=binomial,model=FALSE)

# ridge of selected predictors
ridge.base<-cv.glmnet(x=data.matrix(dataset[,names(dataset) %in% glm_features_100]),
                 y=data.matrix(dataset$outcome),
                 alpha=0,family=c("binomial"),
                 standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)

# lasso of selected predictors
lasso.base<-cv.glmnet(x=data.matrix(dataset[,names(dataset) %in% glm_features_100]),
                 y=data.matrix(dataset$outcome),
                 alpha=1,family=c("binomial"),
                 standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)

# save models
saveRDS(glm_100.base,file="trd_full_glm100.RDS")
saveRDS(glm_int.base,file="trd_full_glmint.RDS")
saveRDS(ridge.base,file="trd_full_ridge.RDS")
saveRDS(lasso.base,file="trd_full_lasso.RDS")

# calculate predictions
preds_glm100<-predict(glm_100.base,dataset[glm_features_100],
                               type="response",standardize=TRUE)
preds_glm_int<-predict(glm_int.base,dataset[glm_features_100],
                                type="response",standardize=TRUE)
preds_ridge<-predict(ridge.base,data.matrix(dataset[glm_features_100]),
                              type="response",standardize=TRUE)
preds_lasso<-predict(lasso.base,data.matrix(dataset[glm_features_100]),
                                  type="response",standardize=TRUE)

# calculate base values
# roc
glm_100.base_roc <- as.numeric(unlist(calculate_auc(preds_glm100,dataset$outcome)@y.values))
glm_int.base_roc <- as.numeric(unlist(calculate_auc(preds_glm_int,dataset$outcome)@y.values))
ridge.base_roc <- as.numeric(unlist(calculate_auc(preds_ridge,dataset$outcome)@y.values))
lasso.base_roc <- as.numeric(unlist(calculate_auc(preds_lasso,dataset$outcome)@y.values))

# pr
glm_100.base_pr <- PRAUC(preds_glm100,dataset$outcome)
glm_int.base_pr <- PRAUC(preds_glm_int,dataset$outcome)
ridge.base_pr <- PRAUC(preds_ridge,dataset$outcome)
lasso.base_pr <- PRAUC(preds_lasso,dataset$outcome)

# decile
glm_100.base_dec <- top_decile_pct(preds_glm100,dataset$outcome)
glm_int.base_dec <- top_decile_pct(preds_glm_int,dataset$outcome)
ridge.base_dec <- top_decile_pct(preds_ridge,dataset$outcome)
lasso.base_dec <- top_decile_pct(preds_lasso,dataset$outcome)

# Dxy
glm_100.base_dxy <- as.numeric(val.prob(preds_glm100,dataset$outcome,pl=F)[1])
glm_int.base_dxy <- as.numeric(val.prob(preds_glm_int,dataset$outcome,pl=F)[1])
ridge.base_dxy <- as.numeric(val.prob(preds_ridge,dataset$outcome,pl=F)[1])
lasso.base_dxy <- as.numeric(val.prob(preds_lasso,dataset$outcome,pl=F)[1])

# set parameters and create results vectors
bootstraps <- 200

# ROC
glm_100.roc <- rep(NA,bootstraps)
glm_int.roc <- rep(NA,bootstraps)
ridge.roc <- rep(NA,bootstraps)
lasso.roc <- rep(NA,bootstraps)

# PR
glm_100.pr <- rep(NA,bootstraps)
glm_int.pr <- rep(NA,bootstraps)
ridge.pr <- rep(NA,bootstraps)
lasso.pr <- rep(NA,bootstraps)

# Decile
glm_100.dec <- rep(NA,bootstraps)
glm_int.dec <- rep(NA,bootstraps)
ridge.dec <- rep(NA,bootstraps)
lasso.dec <- rep(NA,bootstraps)

# Dxy
glm_100.dxy <- rep(NA,bootstraps)
glm_int.dxy <- rep(NA,bootstraps)
ridge.dxy <- rep(NA,bootstraps)
lasso.dxy <- rep(NA,bootstraps)

# intiate experiment
for (idx in 1:bootstraps){
  
  # bootstrap
  data.boot=dataset[sample(nrow(dataset),replace = TRUE),]

  # glm model w/o interactions
  glm_100<-glm(outcome~.,data.boot[c("outcome",glm_features_100)],family=binomial,model=FALSE)

  # include interaction terms
  glm_int<-glm(outcome~.^2,data.boot[c("outcome",glm_features_100)],family=binomial,model=FALSE)

  # ridge of selected predictors
  ridge<-cv.glmnet(x=data.matrix(data.boot[,names(data.boot) %in% glm_features_100]),
                      y=data.matrix(data.boot$outcome),
                      alpha=0,family=c("binomial"),
                      standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)

  # lasso of selected predictors
  lasso<-cv.glmnet(x=data.matrix(data.boot[,names(data.boot) %in% glm_features_100]),
                      y=data.matrix(data.boot$outcome),
                      alpha=1,family=c("binomial"),
                      standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)
  
  # calculate roc aucs
  glm_100.roc[idx] <- roc_diff(glm_100)
  glm_int.roc[idx] <- roc_diff(glm_int)
  ridge.roc[idx] <- roc_diff(ridge,matrix=T)
  lasso.roc[idx] <- roc_diff(lasso,matrix=T)
  
  # calculate pr aucs
  glm_100.pr[idx] <- pr_diff(glm_100)
  glm_int.pr[idx] <- pr_diff(glm_int)
  ridge.pr[idx] <- pr_diff(ridge,matrix=T)
  lasso.pr[idx] <- pr_diff(lasso,matrix=T)
  
  # calculate dec pcts
  glm_100.dec[idx] <- dec_diff(glm_100)
  glm_int.dec[idx] <- dec_diff(glm_int)
  ridge.dec[idx] <- dec_diff(ridge,matrix=T)
  lasso.dec[idx] <- dec_diff(lasso,matrix=T)
  
  # calculate Dxys
  glm_100.dxy[idx] <- cal_diff(glm_100)
  glm_int.dxy[idx] <- cal_diff(glm_int)
  ridge.dxy[idx] <- cal_diff(ridge,matrix=T)
  lasso.dxy[idx] <- cal_diff(lasso,matrix=T)
  
  print(idx)
}

# save results
glm_100.results <- data.frame(glm_100.roc,glm_100.pr,glm_100.dec,glm_100.dxy)
glm_int.results <- data.frame(glm_int.roc,glm_int.pr,glm_int.dec,glm_int.dxy)
ridge.results <- data.frame(ridge.roc,ridge.pr,ridge.dec,ridge.dxy)
lasso.results <- data.frame(lasso.roc,lasso.pr,lasso.dec,lasso.dxy)

write.csv(glm_100.results,file='glm100_bootval_results.csv')
write.csv(glm_int.results,file='glmint_bootval_results.csv')
write.csv(ridge.results,file='ridge_bootval_results.csv')
write.csv(lasso.results,file='lasso_bootval_results.csv')
