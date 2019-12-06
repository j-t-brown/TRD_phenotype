### TRD Rotation Project ###
# boLASSO feature selection

# import TRD data
dataset<-read.csv('',header=T)

# data processing
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

# helper function to reduce RxNorm dimensionality
thresh_codes<-function(input_data,threshold,binarize=FALSE){
  # returns the index of column that will include a % of all
  # of all the values in the dataframe
  # the threshold is the %
  # binarize converts all non-zero counts to 1
  if (binarize==TRUE){
    data=data.frame(input_data)
    data[data!=0]=1
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

# create dataframe to store model coefficients
add_zeros<-function(names,size){
  # assigns all NA to empty dataframe
  X<-rep(NA,size)
  dframe<-data.frame(X)
  for (name in names){
    dframe[,name]<-NA
  }
  dframe[,"(Intercept)"]<-NA
  return(dframe)
}

# set parameters
loops<-100
featureVector<-'DEMO__AGE|DEMO__GENDER__|DEMO__ADI|DEMO__RACE__|ELIX|ICD_CCS_'
feat=names(dataset)[grep(featureVector,names(dataset))]
all_rx<-names(dataset)[grep('RXNORM_',names(dataset))]
log1p_feat=names(dataset)[grep("RXNORM|CCS",names(dataset))]
alph<-1 # lasso

# run experiment
coeffs<-add_zeros(names(dataset),loops)
for (i in 1:loops){
  
  # bootstrap training dataset
  d_train=dataset[sample(nrow(dataset),replace = TRUE),] 
  #d_test=dataset[!(dataset$PERSON_ID %in% d_train$PERSON_ID),]
  
  # find features
  rx_feat<-thresh_codes(d_train[,all_rx],0.95,binarize=FALSE)
  features<-c(feat,rx_feat)
  
  # transform ICD and RXNORM counts
  d_train[log1p_feat]<-log1p(d_train[log1p_feat])
  
  # generate model
  lasall_glmnet<-cv.glmnet(x=data.matrix(d_train[,names(d_train) %in% features]),
                           y=data.matrix(d_train$outcome),
                           alpha=alph,family=c("binomial"),
                           standardize=TRUE,parallel=TRUE,grouped=FALSE,nfolds=10)
  
  # save non-zero coefficients
  non_zero<-names(coef(lasall_glmnet)[which(coef(lasall_glmnet)!=0),1])
  for (c in 1:length(non_zero)){
    coeffs[i,non_zero[c]]<-coef(lasall_glmnet)[non_zero[c],]
  }
  #print(i)
  # save progress
  if (i %% 10 == 0){
    write.csv(coeffs,file="bolasso_coeffs_100.csv")}
}

# helper fucntion to extract relevant features 
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

# extract features
glm_features_100 <- extract_features(coeffs,1.0) #100%

# generate correlation plot
library(corrplot)
corr_matrix<- rcorr(as.matrix(dataset.train[glm_features_100]))
corrplot(corr_matrix$r, type="upper",order = "hclust",
         p.mat=corr_matrix$P, insig="blank", sig.level=0.01,
         tl.cex=0.75, tl.srt=45, tl.col='black')



