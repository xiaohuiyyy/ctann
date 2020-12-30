setwd("D:/nonpthesis")
cancer <- read.table("breast-cancer-wisconsin.data",sep=",")
colnames(cancer) <- c("CodeNum", "ClumpThick", 
                      "CellSize", "CellShape", 
                      "MA", "ECS", 
                      "BN", "BC",
                      "NN", 'Mitoses',"Class")
cancer <- cancer[,-1]
## ========== data preprocess ========== ##
## There are missing values in the 7th columns.
## First, we use nearest neighbors to complete the missing.
nonmissidx <- c(1:5,7:9)
missing_data <- cancer[cancer[,6]=="?",]
library(FNN)
nnformiss <- get.knnx(cancer[, nonmissidx], 
                      cancer[cancer[,6]=="?", nonmissidx],k=3)
comp <- cancer[nnformiss$nn.index[,2],6]
nn3 <- cancer[nnformiss$nn.index[,3],6]
for(i in 1:length(comp)){
  if(comp[i]=="?") comp[i] <- nn3[i]
}
cancer[cancer[,6]=="?", 6] <- comp
cancer[,6] <- as.integer(cancer[,6])

## ========== Binary Classifier ========== ##
library(caTools)
library(tcltk)

## ========== Kernal Regression ========== ##
library(mvtnorm)
K <- function(x, h) dmvnorm(x, mean=rep(0,9), diag(h,9)) # d-dimensional kernel
m.hat <- function(x, X, y, h){
  n <- nrow(X); 
  f <- sum(sapply(1:n, function(i){
    K(x-X[i,],h)
  }))
  r <- sum(sapply(1:n, function(i){
    K(x-X[i,],h)*y[i]
  }))
  return(r/f)
}
simreg <- function(o){
  split = sample.split(cancer$Class, SplitRatio = 0.75)
  training_set = subset(cancer, split == TRUE)
  testing_set = subset(cancer, split == FALSE)
  y_predkreg <- numeric(nrow(testing_set))
  for(i in 1:nrow(testing_set)){
    y_predkreg[i] <- m.hat(testing_set[i, -10], 
                           training_set[,-10],training_set[,10],
                           h=1)
  }
  y_predkreg <- ifelse(y_predkreg>3,4,2)
  t <- table(testing_set[, 10], y_predkreg)
  return(t)
}
errreg <- list(0)
pb <- tkProgressBar("kernel regression progress","completed %", 0, 100) 
star_time <- Sys.time()
for(j in 1:50){
  errreg[[j]] <- simreg(1)
  info <- sprintf("completed %d%%", round(j*100/50))
  setTkProgressBar(pb, j*100/50, sprintf("kernel regression process (%s)", info),info)
}
end_time <- Sys.time()
close(pb)  

## ========== k nearest neighbor regression ========== ##
library(FNN)
simknn <- function(o){
  split = sample.split(cancer$Class, SplitRatio = 0.75)
  training = subset(cancer, split == TRUE)
  testing_set = subset(cancer, split == FALSE)
  k <- 5
  nntesting <- get.knnx(training[, -10], 
                        testing_set[,-10],k=k)
  y_predknn <- numeric(nrow(testing_set))
  for(i in 1:nrow(testing_set)){
    y_predknn[i] <- mean(training[nntesting$nn.index[i,], 10])
  }
  y_predknn <- ifelse(y_predknn>3,4,2)
  t <- table(testing_set[, 10], y_predknn)
  return(t)
}

errknn <- list(0);
pb <- tkProgressBar("knn progress","completed %", 0, 100) 
star_time <- Sys.time()
for(j in 1:10){
  errknn[[j]] <- simknn(1)
  info <- sprintf("completed %d%%", round(j*100/10))
  setTkProgressBar(pb, j*100/10, sprintf("knn process (%s)", info),info)
}
end_time <- Sys.time()
close(pb)
end_time-star_time
errstat(errknn)

## ========== Neural Network ========== ##
library(neuralnet)

simnn <- function(o){
  split = sample.split(cancer$Class, SplitRatio = 0.75)
  training_set = subset(cancer, split == TRUE)
  testing_set = subset(cancer, split == FALSE)
  training_set[,10] <- training_set[,10]/2-1
  nn1=neuralnet(Class~.,data=training_set, 
                hidden=23 ,act.fct = "logistic", 
                linear.output = FALSE)
  prednn1=compute(nn1,testing_set[,-10])
  prob1 <- prednn1$net.result
  y_prednn1 <- ifelse(prob1>0.5, 1, 0)
  t1 <- table(testing_set[,10], y_prednn1)
  nn2=neuralnet(Class~.,data=training_set, 
                hidden=c(9,6) ,act.fct = "logistic", 
                linear.output = FALSE)
  prednn2=compute(nn2,testing_set[,-10])
  prob2 <- prednn2$net.result
  y_prednn2 <- ifelse(prob2>0.5, 1, 0)
  t2 <- table(testing_set[,10], y_prednn2)
  return(list(t1,t2))
}

errnn1 <- errnn2 <- list(0)
pb <- tkProgressBar("neural network progress","completed %", 0, 100) 
star_time <- Sys.time()
for(j in 1:10){
  rst <- simnn(1)
  errnn1[[j]] <- rst[[1]]
  errnn2[[j]] <- rst[[2]]
  info <- sprintf("completed %d%%", round(j*100/10))
  setTkProgressBar(pb, j*100/10, sprintf("neural network process (%s)", info),info)
}
end_time <- Sys.time()
close(pb)
end_time-star_time
errstat(errnn1)
errstat(errnn2)

## ========== Decision Tree ====== ##
library(rpart)

simct <- function(o){
  split = sample.split(cancer$Class, SplitRatio = 0.75)
  training_set = subset(cancer, split == TRUE)
  testing_set = subset(cancer, split == FALSE)
  training_set$Class <- as.factor(ifelse(training_set$Class==4, "Positive", "Negative"))
  testing_set$Class <- as.factor(ifelse(testing_set$Class==4, "Positive", "Negative"))
  output.tree <- rpart(Class~.,data=training_set,method="class")
  prob <- predict(output.tree, testing_set[,-10])
  y_predcttree <- as.integer(prob[,1]<=prob[,2])
  t <- table(testing_set$Class, y_predcttree)
  return(t)
}
errct <- list(0)
pb <- tkProgressBar("classifier tree progress","completed %", 0, 100) 
star_time <- Sys.time()
for(j in 1:10){
  rst <- simct(1)
  errct[[j]] <- rst
  info <- sprintf("completed %d%%", round(j*100/10))
  setTkProgressBar(pb, j*100/10, sprintf("classifier tree process (%s)", info),info)
}
end_time <- Sys.time()
close(pb)
end_time-star_time
errstat(errct)

## ========== Kernel SVM ====== ##
library(e1071)

simsvm <- function(o){
  split = sample.split(cancer$Class, SplitRatio = 0.75)
  training_set = subset(cancer, split == TRUE)
  testing_set = subset(cancer, split == FALSE)
  classifierR = svm(formula = Class ~ .,
                    data = training_set,
                    type = 'C-classification', 
                    # this is because we want to make a regression classification
                    kernel = 'radial')
  y_predR = predict(classifierR, newdata = testing_set[,-10])
  t <- table(testing_set[, 10], y_predR)
  return(t)
}
errsvm <- list(0)
pb <- tkProgressBar("svm progress","completed %", 0, 100) 
star_time <- Sys.time()
for(j in 1:10){
  errsvm[[j]] <- simsvm(1)
  info <- sprintf("completed %d%%", round(j*100/10))
  setTkProgressBar(pb, j*100/10, sprintf("svm process (%s)", info),info)
}
end_time <- Sys.time()
close(pb)
end_time-star_time
errstat(errsvm)

## ========== CT-ANN ========== ##
sigmoid = function(x) {
  1 / (1 + exp(-x))
}
simctann <- function(o){
  split <-  sample.split(cancer$Class, SplitRatio = 0.75)
  training = subset(cancer, split == TRUE)
  testing_set = subset(cancer, split == FALSE)
  training$Class <- as.factor(ifelse(training$Class==4, "Positive", "Negative"))
  testing_set$Class <- as.factor(ifelse(testing_set$Class==4, "Positive", "Negative"))
  dtree <- CVVS(training)
  y_predtrainct <- predict(dtree,newdata=training[,-10],type="class")
  y_predtrainct <- ifelse(y_predtrainct=="Positive",1,0) 
  training$ctclass <- y_predtrainct
  y_predtestct <- predict(dtree,newdata=testing_set[,-10],type="class")
  y_predtestct <- ifelse(y_predtestct=="Positive",1,0)
  testing_set$ctclass <- y_predtestct 
  impvar <- names(dtree$variable.importance)
  impvar <- impvar[1:min(7,length(impvar))]
  nn <- neuralnet(Class~., data=training_set[,c(impvar,"ctclass","Class")], 
                  hidden=numhn, act.fct = sigmoid,
                  linear.output = FALSE)
  prob_ctann <- predict(nn, testing_set[,c(impvar,"ctclass")])
  y_pred_ctann <- as.integer(prob_ctann[,1] <= prob_ctann[,2])
  table(testing_set$Class, y_pred_ctann)
}
CVVS <- function(training){
  CV <- list(0)
  for(i in 1:10){
    split2 <-  sample.split(training$Class, SplitRatio = 0.67)
    training_set = subset(training, split2 == TRUE)
    validation_set = subset(training, split2 == FALSE)
    dtree<-rpart(Class~.,data=training_set,method="class", 
                 parms=list(split="gini"), minsplit=10)
    y_predtrainct <- predict(dtree,newdata=training_set[,-10],type="class")
    y_predtrainct <- ifelse(y_predtrainct=="Positive",1,0) 
    training_set$ctclass <- y_predtrainct
    ttrain <- table(training_set$Class,training_set$ctclass)
    y_predvalidct <- predict(dtree,newdata=validation_set[,-10],type="class")
    y_predvalidct <- ifelse(y_predvalidct=="Positive",1,0) 
    validation_set$ctclass <- y_predvalidct
    tvalid <- table(validation_set$Class,validation_set$ctclass)
    L <- list(dtree, sum(tvalid)-sum(diag(tvalid)))
    CV[[i]] <- L
  }
  errnum <- sapply(1:length(CV), function(j){CV[[j]][[2]]})
  idx <- which.min(errnum)
  return(CV[[idx]][[1]])
}

errctann <- list(0)
pb <- tkProgressBar("CT-ANN progress","completed %", 0, 100) 
star_time <- Sys.time()
for(j in 1:10){
  errctann[[j]] <- simctann(1)
  info <- sprintf("completed %d%%", round(j*100/10))
  setTkProgressBar(pb, j*100/10, sprintf("CT-ANN process (%s)", info),info)
}
end_time <- Sys.time()
close(pb)
end_time-star_time
errstat(errctann)

## ========== Compare ========== ##
errcomp <- function(err){
  n <- length(err)
  TN <- FP <- TP <- FN <- numeric(n)
  for(j in 1:n){
    TN[j] <- err[[j]][1,1]
    TP[j] <- err[[j]][2,2]
    FN[j] <- err[[j]][2,1]
    FP[j] <- err[[j]][1,2]
  }
  return(data.frame(TN,TP,FN,FP))
}
errstat <- function(err){
  d <- errcomp(err)
  r <- (d[,3]+d[,4])/174
  return(100*c(1-mean(r),sd(r)))
}
errF <- function(err){
  d <- errcomp(err)
  mean(sapply(1:nrow(d), function(j){
    tab <- d[j,]
    precision <- tab$TP/(tab$TP+tab$FP)
    recall <- tab$TP/(tab$TP+tab$FN)
    2*precision*recall/(precision+recall)
  }))
}

acc <- t(matrix(c(errstat(errknn),
         errstat(errnn1),
         errstat(errnn2),
         errstat(errct),
         errstat(errsvm),
         errstat(errctann)),nrow=2))

Fmea <- c(errF(errknn),
       errF(errnn1),
       errF(errnn2),
       errF(errct),
       errF(errsvm),
       errF(errctann))

tab <- data.frame(acc[,1],acc[,2], Fmea)


