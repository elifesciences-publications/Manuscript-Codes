#pipeline: 
#1. run wards heirarchical clustering first, relabel cell class based on these results
#2. run decision tree, create model using 80% training data, test on 20% test data, log accuracy
#3. optionally cross validate, OR
#4. run random forest, log accuracy from confusion matrix, plot and log gini/accuracy indices, plot MDS plot
#based off of proximity matrix




install.packages("mclust", dependencies = TRUE)
install.packages("flexmix", dependencies = TRUE)
install.packages("prabclus", dependencies = TRUE)
install.packages("scales", dependencies = TRUE)
install.packages("lazyeval", dependencies = TRUE)
install.packages("robustbase", dependencies = TRUE)
install.packages("kernlab", dependencies = TRUE)
install.packages("dendextend", dependencies = TRUE)
#load libraries
library(RGtk2)
library(rattle)
citation(package = "rattle", lib.loc = NULL)
library(rpart.plot)
citation(package = "rpart.plot", lib.loc = NULL)
library(RColorBrewer)
library(party)
library(partykit)
library(caret)
library(rpart)
citation(package = "rpart", lib.loc = NULL)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(randomForest)
library(treeinterpreter)
citation(package = "randomForest", lib.loc = NULL)

library(scales)
library(dendextend)
library(dplyr)
library(flexmix)
#excluding RS
#setwd("/Users/chuph/Desktop")
setwd("C:/Users/pchu1207/Google Drive/Data/Mae/")
getwd()
alllhx6young <- read.csv("cluster results.csv")
younghet <- alllhx6young

#select parameters het
younghet2 <- subset(younghet, select= c(cell, cell.index, wards.class, mem.resistance, max.fr,delta.ahp,sag,	halfwidth,	ahp.magnitude,ahp.latency,	ap.threshold,	fr.accom,	ap.amp.accom,	latency.to.first.ap,	ap.broadening,	ap.amp, number.rebound.ap,	membrane.tau))
sum(is.na(younghet2))
younghet3 <- na.omit(younghet2)
head(younghet3)
length(younghet3)
# n = 
sum(younghet3$wards.class=="1")
sum(younghet3$wards.class=="2")



#decision tree
#Randomness is used for generating the cross validation errors. We might notice that on separate runs of rpart with exactly the 
#same settings the relative errors are consistent because there is no random sampling in this calculation. 
#However the xerror and xtsd will vary unless we set the random number seed to the same value each time. 
#So we can only accurately compare when we are sure we are using the same random number sequence. Rattle specifically sets the seed each 
#time so that a user will not be puzzled about slightly different results every time they build the tree. 
require(devtools)
install_version("arules", version = "1.5-4")
install.packages("RGtk2", dependencies=c("Depends", "Suggests"))
install.packages("rattle", repos="http://rattle.togaware.com", type="source")
install.packages("https://cran.r-project.org/bin/windows/contrib/3.3/RGtk2_2.20.31.zip", repos=NULL)
install.packages("treeinterpreter", dependencies = TRUE)


setwd("/Users/pchu1207/Google Drive/Data/Mae")
alllhx6young <- read.csv("cluster results.csv")

younghet2 <- subset(alllhx6young, select= c(cell, cell.index, wards.class, mem.resistance, max.fr,delta.ahp,sag,	halfwidth,	ahp.magnitude,ahp.latency,	ap.threshold,	fr.accom,	ap.amp.accom,	latency.to.first.ap,	ap.broadening,	ap.amp, number.rebound.ap,	membrane.tau))
#are there any NAs (blanks) in data set?
sum(is.na(younghet2))
#exclude rows with NAs
younghet3 <- na.omit(younghet2)
#take a look at df
head(younghet3)
length(younghet3$cell)
#set seed to ensure reproducible results
#set.seed(41), rattle does this
#split into training and test sets, 80% goes to training set, do not replace after sampling a row
train <- sample_frac(younghet3, 0.8, replace=FALSE)
#check number sampled
length(train$cell)
#create test set
test <- younghet3[!(younghet3$cell.index %in% train$cell.index),]
#test n of test set
length(test$cell.index)
set.seed(1234)
#choose minsplit=0 and cp=0 to determine best tree fit. want to limit overfitting. run prp and plot the complexity parameter.
form <- as.formula(wards.class~.)
#run recursive partitioning
#set.seed(1234) , rattle does this
#training data
younghet4 <- rpart(form, data=train[,-1:-2], control=rpart.control(minsplit=0,minbucket=1,cp=0, method="class")) 
#full data set
younghet4 <- rpart(form, data=younghet3[,-1:-2], control=rpart.control(minsplit=0,minbucket=1,cp=0, method="class")) 
prp(younghet4)
#set.seed(1234) , rattle does this
#table with classified (where) vs assumed class (target)
train$where <- younghet4$where
table(train$wards.class, train$where)

younghet4$cptable
train <- train[,-grep("where", colnames(train))]


#check training model against test, confusion matrix
predict <- predict(younghet4, test, type="class")  
test$wards.class
predict
table(pred=predict,true=test$wards.class)
# % averagecorrect prediction, unpruned
mean(predict==test$wards.class)


#plot complexity table for rpart fit, shows reduction in error after each split. start with cp=0,minsplit=0, choose relative error reduction above dotted line (represents standard error above minimum cross validated error)
#set.seed(1234) 
plotcp(younghet4)
printcp(younghet4)
#prune tree based on complexity parameter
pruned <- prune(younghet4, cp=0.044444)
par(cex=1.1)
prp(pruned)



#prune tree based on complexity parameter for training
pruned <- prune(younghet4, cp=0.027778)
par(cex=1.1)
prp(pruned)

predict <- predict(pruned, test, type="class")  
test$wards.class
predict
table(pred=predict,true=test$wards.class)
# % averagecorrect prediction, unpruned
mean(predict==test$wards.class)


#check training model of pruned against test, confusion matrix
table(train$class, pruned$where)
#check training model of pruned against test
predict <- predict(pruned, test, type="class")  
predict
# % averagecorrect prediction
mean(predict==test$wards.class)
#check training model of pruned against test, confusion matrix
table(pred=predict,true=test$wards.class)



#function to do multiple runs
multiple_runs_classification <- function(train_fraction,n,dataset,prune_tree=FALSE){
  fraction_correct <- rep(NA,n)
    set.seed(42)
    for (i in 1:n){
       train <- sample_frac(dataset, 0.8, replace=FALSE)
       test <- dataset[!(dataset$cell.index %in% train$cell.index),]
       form <- as.formula(wards.class~.)
      rpart_model <- rpart(form,data=train[,-1], control=rpart.control(minsplit=1,minbucket=1,cp=0, method="class")) 
    if(prune_tree==FALSE) {
      rpart_test_predict <- predict(rpart_model, test, type="class")  
      fraction_correct[i] <- mean(rpart_test_predict==test$wards.class)
    }else{
      opt <- which.min(rpart_model$cptable[,"xerror"])
      cp <- rpart_model$cptable[opt, "CP"]
      #pruned_model <- prune(rpart_model,cp)
      pruned_model <- prune(rpart_model,0.044444) # CP determined from previous
      rpart_pruned_predict <-predict(pruned_model, test, type="class")  
      fraction_correct[i] <- mean(rpart_pruned_predict==test$wards.class)
    }
  }
  return(fraction_correct)
}

unpruned <- multiple_runs_classification(0.8,500,younghet3,prune_tree = FALSE)
unpruned
mean(unpruned)
sd(unpruned)
pruned <- multiple_runs_classification(0.8,500,younghet3, prune_tree = TRUE)
mean(pruned)
sd(pruned)






#######################################################################################################

#random forest
setwd("/Users/pchu1207/Google Drive/Data/Mae")
alllhx6young <- read.csv("Copy of sstClustering.csv")
alllhx6young <- read.csv("cluster results.csv")

str(alllhx6young)
head(alllhx6young)
younghet <- alllhx6young

#select parameters het, using response to min stim as a factor
#responding.in.min.stim.experiment
younghet2 <- subset(alllhx6young, select= c(cell, cell.index, wards.class, mem.resistance, max.fr,delta.ahp,sag,	halfwidth,	ahp.magnitude,ahp.latency,	ap.threshold,	fr.accom,	ap.amp.accom,	latency.to.first.ap,	ap.broadening,	ap.amp, number.rebound.ap,	membrane.tau))
#are there any NAs (blanks) in data set?
sum(is.na(younghet2))
#exclude rows with NAs
younghet3 <- na.omit(younghet2)
#take a look at df
head(younghet3)
#training data
train <- sample_frac(younghet3, 0.8, replace=FALSE)
#check number sampled
length(train$cell)
#create test set
test <- younghet3[!(younghet3$cell.index %in% train$cell.index),]
#test n of test set
length(test$cell.index)
#set starting point for random sampling
set.seed(1234)

str(younghet3)
#tune the tree 
#AThe code below allows us to choose the number of variables that minimizes our out-of-bag error: 
#the rate of incorrect classifications for those left out of each tree. s you can see, using two or 
#three variables for each tree-chosen at random by randomForest-is the best choice 
# Therefore, we include this as our mtry in the fitting of the forest.
tuning <- tuneRF(younghet3[,c( "mem.resistance", "max.fr","delta.ahp","sag",	"halfwidth",	"ahp.magnitude","ahp.latency",	"ap.threshold",	"fr.accom",	"ap.amp.accom",	"latency.to.first.ap",	"ap.broadening",	"ap.amp", "number.rebound.ap",	"membrane.tau")], factor(younghet3[,c("wards.class")]), 
ntreeTry=10000, data=younghet3, plot=T)
tuning

# two variables is lowest with current data 4/22/19

#random forest
#random forest
#random forest
younghet4rf <- randomForest(wards.class~.,data=younghet3[,-1:-2], ntree=10000, importance=TRUE, proximity=TRUE, mtry=6)
younghet4rf$err.rate
plot(younghet4rf)
younghet4rf

#####################this is only for double validation#####################################
#confusion matrix
#check full model against test, confusion matrix
predict <- predict(younghet4rf, test, type="class")  
# % averagecorrect prediction, of course this will be 100% because model was built off of it + the rest
mean(predict==test$wards.class)
#
table(predict,test$wards.class)




#resample a test set to check accuracy
#can rerun this many times to get an idea of prediction accuracy
train <- sample_frac(younghet3, 0.8, replace=FALSE)
test <- younghet3[!(younghet3$cell.index %in% train$cell.index),]
predict <- predict(younghet4rf, test, type="class")  
predict
# % averagecorrect prediction
mean(predict==test$wards.class)
table(predict,test$wards.class)

?sample_frac


#function to do multiple runs, e.g. create multiple models from training data test them each time
#gives you an idea of how consistent the model creation across training sets are, dont really need this
#for randomforest because it is inherently testing against OOB 
multiple_runs_classification <- function(train_fraction,test_fraction,n,dataset,prune_tree=FALSE){
  fraction_correct <- rep(NA,n)
  set.seed(42)
  for (i in seq(1,n,1)){
    train <- sample_frac(dataset, train_fraction, replace=FALSE)
    rpart_model <- randomForest(wards.class~.,data=train[,-1:-2], ntree=10000, importance=TRUE, proximity=TRUE, mtry=2)
    test <- younghet3[!(younghet3$cell.index %in% train$cell.index),]
    #form <- as.formula(wards.class~.)
      rpart_test_predict <- predict(rpart_model, test, type="class")  
      fraction_correct[i] <- mean(rpart_test_predict==test$wards.class)

    }
  print(fraction_correct)
}
#not sure why this doesnt work, resulting in 100% correct predictions for all the newly created 50 test sets
unpruned <- multiple_runs_classification(train_fraction=0.8,test_fraction=0.2, n=20,dataset=younghet3)
unpruned
mean(unpruned)
#####################end of double validation#####################################



#simple explanation of margin: In ensemble classification, you mostly do a majority vote from all models in the 
#ensemble, with the class voted most becoming the finally predicted class. The margin thereby indicates the ratio
#of votes having been done for the correct class for all samples. 1 indicates that for one sample, all votes of 
#the ensemble were correct, while e.g. 0 indicates a draw between the correct and the next best classes. 
#Therefore, values >0 mean that the majority was right, hence this sample ended up as being predicted correct - 
#whilst all values <0 mean that the majority was wrong, hence the sample ended up as being classified wrong. 
#Again, colors indicate classes, so in the example above you see that nearly all setosa samples got classified 
#correctly, but for some of the virginca and versicolor samples the ensemble was not so sure anymore (but they 
#still got the final result correct), while for 4 or 5 of them the final result was plain wrong (which is to 
#expected this way with this dataset).
plot(margin(younghet4rf))
legend("topleft",
       legend=levels(younghet4rf$predicted),
       fill=brewer.pal(length(levels(younghet4rf$predicted)),
                       "Set1"))
#if proximity=TRUE when randomForest is called, a matrix of proximity measures among the input (based on the 
#frequency that pairs of data points are in the same terminal nodes)."
#The term "proximity" means the "closeness" or "nearness" between pairs of cases.
#Proximities are calculated for each pair of cases/observations/sample points. If two cases occupy the same 
#terminal node through one tree, their proximity is increased by one. At the end of the run of all trees, the 
#proximities are normalized by dividing by the number of trees. Proximities are used in replacing missing data,
#locating outliers, and producing illuminating low-dimensional views of the data. 
#https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
#MDS can only base its analysis on the output of your randomForest. If you're expecting a better separation, 
#then you might want to check the classification performance of your randomForest. Another thing to keep in mind 
#is that your PCA is mapping from 9-dimensional data to 2 dimensions, but the MDS is mapping from an NxN-dimensional 
#proximity matrix to 2 dimensions, where N is the number of datapoints.
younghet4rf$proximity
#distance plot
par(cex=1.0)
MDSplot(younghet4rf, younghet3$wards.class,k=2)
legend("topleft",
       legend=levels(younghet4rf$predicted),
       fill=brewer.pal(length(levels(younghet4rf$predicted)),
                       "Set1"), cex=0.75)

##
#manual distance plot with individual points named and colored
multi<- cmdscale(1-younghet4rf$proximity,k=2)
multi
younghet4rf$x <- multi[,1]
younghet4rf$y <- multi[,2]

#plot multidimensional scaled data, first two eigenvectors with cell indices
younghet3$classB <- younghet3$wards.class
plot(younghet4rf$x, younghet4rf$y, xlab="dim1", ylab="dim2", type="n", main="ok")
text(younghet4rf$x, younghet4rf$y, younghet3$cell.index, cex=.75, col=younghet3$classB)
legend(-0.5,0.45,c("a","b"), fill=c("black","red"), bty="n", cex=1.2)
##



#look at "tree" table
getTree(younghet4rf, k=2)
#same tables from above
younghet4rf$confusion
younghet4rf$predicted
#percent correct prediction from random forest
mean(younghet4rf$predicted==younghet3$responding.in.min.stim.experiment)


 #table of gini, gini is metric for impurity/ variability
#variable importance, prediction accuracy
importance(younghet4rf)
#plots OOB (out of bag observations, i.e. test cases), impurity measure (gini) as a function of parameters. high score means the parameter being removed decreased the purity by alot (its important)
varImpPlot(younghet4rf)

#predict, Here we can see the proportion of trees that each neuron is in a class. 
predicted.vote <- predict(younghet4rf, test, type="vote", norm.votes="TRUE")
predicted.vote
#attach predicted class and vote percentages from model
test$class2 <- as.numeric(predict)
test$votes <- predicted.vote[,2]
#
test[order(-test$votes),]


length(younghet3$cell.index)

#
#
#
#
#

#wards cluster analysis
setwd("/Users/pchu1207/Google Drive/Data/Mae")
alllhx6young <- read.csv("Copy of sstClustering.csv")
alllhx6young <- read.csv("cluster results.csv")

str(alllhx6young)
head(alllhx6young)
younghet <- alllhx6young

#select parameters het
younghet2 <- subset(younghet, select= c(cell, cell.index, wards.class, mem.resistance, max.fr,delta.ahp,sag,	halfwidth,	ahp.magnitude,ahp.latency,	ap.threshold,	fr.accom,	ap.amp.accom,	latency.to.first.ap,	ap.broadening,	ap.amp, number.rebound.ap,	membrane.tau))
#are there any NAs (blanks) in data set?
sum(is.na(younghet2))
#exclude rows with NAs
younghet3 <- na.omit(younghet2)
younghet2 <- subset(younghet, select= c(cell.index,wards.class, mem.resistance, max.fr,delta.ahp,sag,	halfwidth,	ahp.magnitude,ahp.latency,	ap.threshold,	fr.accom,	ap.amp.accom,	latency.to.first.ap,	ap.broadening,	ap.amp, number.rebound.ap,	membrane.tau))
younghet3 <- na.omit(younghet2)
length(unique(younghet3[1,]))
head(younghet3)

names <- younghet3$wards.class # change what is displayed as labels for the branch tips
rownames(younghet3) <- make.names(names,unique=TRUE)
younghet4 <- scale(younghet3[,-1:-3])
younghet4

younghet6 <- as.dendrogram(hclust(dist(younghet4, method="euclidean"), method="ward.D"))
colors <- as.numeric(younghet3[,"wards.class"])
colors <- colors[order.dendrogram(younghet6)]
labels_colors(younghet6) <- colors
par(cex=0.5)
plot(younghet6, lty=2)
younghet6 %>% set("branches_k_color", k = 2) %>%  
  plot()    


#% correctly classified by wards
sum(younghet3$wards.class=="b")
100-(2/sum(younghet3$wards.class=="b"))
sum(younghet3$wards.class=="a")
100-(3/sum(younghet3$wards.class=="a"))

#mann whitney
options(digits=10)
mannwhitney <- apply(younghet3[,-1:-4],2,function(x) wilcox.test(paired=FALSE, x~responding.in.min.stim.experiment, data=younghet3)$p.value)
mannwhitney1 <- p.adjust(mannwhitney, method = "fdr")
print(as.data.frame(mannwhitney1))


