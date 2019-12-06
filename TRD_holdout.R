### TRD Rotation Project ###
# holdout experiment

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

# set parameters for LASSO
featureVector<-'DEMO__AGE|DEMO__GENDER__|DEMO__ADI|DEMO__RACE__|ELIX|ICD_CCS_'
feat=names(dataset)[grep(featureVector,names(dataset))]
all_rx<-names(dataset)[grep('RXNORM_',names(dataset))]
log1p_feat=names(dataset)[grep("RXNORM|CCS",names(dataset))]
alph<-1 # lasso

# 80-20 hold out
dataset.train=dataset[sample(nrow(dataset),round(nrow(dataset)*80/100),replace=FALSE),]
dataset.test=dataset[!(dataset$PERSON_ID %in% dataset.train$PERSON_ID),]

# find features
thresh_codes<-function(input_data,threshold,binarize=FALSE){
  # returns the RxNorm variable names that include threshold % of all
  # of all the count values in the input_data
  # binarize converts all non-zero counts to 1
  if (binarize==TRUE){
    data=data.frame(input_data)
    data[data>0]=1
    data[data<1]=0 # in case of any negative values
    sum_col<-colSums(data)
    sum_col<-sort(sum_col,decreasing=TRUE)
    cdf<-cumsum(sum_col)/sum(sum_col)
    return(names(which(cdf<=threshold)))
  }
  sum_col<-colSums(input_data)
  sum_col<-sort(sum_col,decreasing=TRUE)
  cdf<-cumsum(sum_col)/sum(sum_col)
  
  return(names(which(cdf<=threshold)))
}

rx_feat<-thresh_codes(dataset.train[,all_rx],0.95,binarize=FALSE)
features<-c(feat,rx_feat)

# fit LASSO model before the log1p transform
lasso_pre<-cv.glmnet(x=data.matrix(dataset.train[,names(dataset.train) %in% features]),
                         y=data.matrix(dataset.train$outcome),
                         alpha=alph,family=c("binomial"),
                         standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)

# calculate LASSO predictions
maindata <- data.frame(dataset.test$outcome)
maindata$preds_lpre<-predict(lasso_pre,data.matrix(dataset.test[features]),
                             type="response",standardize=TRUE)

# transform ICD and RXNORM counts
dataset.train[log1p_feat]<-log1p(dataset.train[log1p_feat])
dataset.test[log1p_feat]<-log1p(dataset.test[log1p_feat])

# generate LASSO model
lasso_post<-cv.glmnet(x=data.matrix(dataset.train[,names(dataset.train) %in% features]),
                         y=data.matrix(dataset.train$outcome),
                         alpha=alph,family=c("binomial"),
                         standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)

# boLASSO selected features for glm model development
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

glm_features_100 <- extract_features(coeffs,1.0) #100%

# glm model w/o interactions
glm_model_100<-glm(outcome~.,dataset.train[c("outcome",glm_features_100)],family=binomial,model=FALSE)

# include interaction terms
glm_int_100<-glm(outcome~.^2,dataset.train[c("outcome",glm_features_100)],family=binomial,model=FALSE)

# ridge of selected predictors
ridge_100<-cv.glmnet(x=data.matrix(dataset.train[,names(dataset.train) %in% glm_features_100]),
                     y=data.matrix(dataset.train$outcome),
                     alpha=0,family=c("binomial"),
                     standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)

# lasso of selected predictors
lasso_100<-cv.glmnet(x=data.matrix(dataset.train[,names(dataset.train) %in% glm_features_100]),
                    y=data.matrix(dataset.train$outcome),
                    alpha=alph,family=c("binomial"),
                    standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)

# calculate remaining predictions
maindata$preds_lpost<-predict(lasso_post,data.matrix(dataset.test[features]),
                              type="response",standardize=TRUE)
maindata$preds_glm100<-predict(glm_model_100,dataset.test[glm_features_100],
                               type="response",standardize=TRUE)
maindata$preds_glm_int<-predict(glm_int_100,dataset.test[glm_features_100],
                              type="response",standardize=TRUE)
maindata$preds_ridge<-predict(ridge_100,data.matrix(dataset.test[glm_features_100]),
                              type="response",standardize=TRUE)
maindata$preds_lasso_100<-predict(lasso_100,data.matrix(dataset.test[glm_features_100]),
                              type="response",standardize=TRUE)

# save data
write.csv(maindata,file="trd_model_results.csv")

# save models
#saveRDS(lasso_pre,file="trd_lasso_pre.rds")
#saveRDS(lasso_post,file="trd_lasso_post.rds")
#saveRDS(glm_model_100,file="trd_glm_100.rds")
#saveRDS(glm_int_100,file="trd_glm_int.rds")
#saveRDS(ridge_100,file="trd_ridge_100.rds")
#saveRDS(lasso_100,file="trd_lasso_100.rds")
