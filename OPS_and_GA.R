

# reduce demensions manually 19*18*10 fold
clean7=read.csv("cleandata1.csv")
clean7[,1]<-NULL
clean7[,20]<-NULL
dd = NULL
for (a in seq(1,19,1)) {
  yourData1=clean7
  yourData1[,a]<-NULL
  for (b in seq(1,18,1)) {
    yourData=yourData1
    yourData[,b]<-NULL
    yourData<-yourData[sample(nrow(yourData)),]
    folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)
    d = NULL 
    for (x in seq(1,10,1)) {
      testIndexes <- which(folds==x,arr.ind=TRUE)
      test7 <- yourData[testIndexes, ]
      train7 <- yourData[-testIndexes, ]
      testdata<-test7 
      test7[ ,18]<-NULL
      library("e1071")
      svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = 0.057)
      pred <- predict(svm_model, newdata = test7)
      prediction<-data.frame(pred)
      rsq <- function(x, y) summary(lm(y~x))$r.squared
      R2=rsq(testdata$solubility, prediction[,1])
      d = rbind(d, data.frame(x,R2))
    }
    average_d=sum(d[1:10,2])/10
    dd = rbind(dd, data.frame(a,b,average_d))
  }
}
write.csv(dd,"reduce2.csv")

#use confusion matrix plot
reduce=read.csv("reduce2.csv")
clean8=read.csv("composition_difference.csv")
clean7=read.csv("rank.csv")
colnames(clean7)<-clean8[1:20,1]
clean7[,20]<-NULL
colnames(clean7)[1]<-"A0"
colnames(clean7)[2]<-"R+"
colnames(clean7)[3]<-"N0"
colnames(clean7)[4]<-"D-"
colnames(clean7)[5]<-"C0"
colnames(clean7)[6]<-"E-"
colnames(clean7)[7]<-"Q0"
colnames(clean7)[8]<-"G0"
colnames(clean7)[9]<-"H+(10%),0(90%)"
colnames(clean7)[10]<-"I0"
colnames(clean7)[11]<-"L0"
colnames(clean7)[12]<-"K+"
colnames(clean7)[13]<-"M0"
colnames(clean7)[14]<-"F0"
colnames(clean7)[15]<-"P0"
colnames(clean7)[16]<-"S0"
colnames(clean7)[17]<-"T0"
colnames(clean7)[18]<-"W0"
colnames(clean7)[19]<-"YO"
dd = NULL
for (a in seq(1,19,1)) {
  yourData1=clean7
  yourData1[,a]<-NULL
  for (b in seq(1,18,1)) {
    dd = rbind(dd, data.frame(a,b,colnames(clean7)[a],colnames(yourData1)[b]))
  }
}
install.packages("ggplot2")
library("ggplot2")
reduce[,2:3]<-dd[,3:4]
write.csv(reduce,"reduce2.csv")
# just for plot to change clumn name
write.csv(reduce,"reduce2colume.csv")
tile <- ggplot() + geom_tile(aes(x=first.column.removed, y=second.column.removed, fill=R2),data=reduce, color="black",size=0.1) +labs(x="First amino acid removed",y="Second amino acid removed")
tile = tile + geom_text(aes(x=first.column.removed,y=second.column.removed, label=sprintf("%.4f", R2)),data=reduce, size=3, colour="black") + scale_fill_gradient(low="red",high="grey")
tile = tile + geom_tile(aes(x=first.column.removed,y=second.column.removed),data=subset(reduce, as.character(first.column.removed)==as.character(second.column.removed)), color="black",size=0.3, fill="black", alpha=0) 

tile


#try svm model
clean7=read.csv("cleandata_1.csv")
yourData=clean7
yourData<-yourData[sample(nrow(yourData)),]
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

d = NULL 
for (x in seq(1,10,1)) {
  testIndexes <- which(folds==x,arr.ind=TRUE)
  test7 <- yourData[testIndexes, ]
  train7 <- yourData[-testIndexes, ]
  testdata<-test7 
  test7[ ,22]<-NULL
  train7[,1]<-NULL
  test7[,1]<-NULL
  library("e1071")
  svm_model <- svm(solubility ~ ., train7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = 0.057)
  pred <- predict(svm_model, newdata = test7)
  prediction<-data.frame(pred)
  d = rbind(d, data.frame(testdata,prediction[,1]))
}

rsq <- function(x, y) summary(lm(y~x))$r.squared
R2=rsq(d[,22], d[,23])  

#-------------------------------------------------------
#OPS algorithm

library("Surrogate")

#initial value for x
rank=read.csv("rank_2.csv")
new <- rank[ which(rank[ ,22]== 0.10),]
sample <- new[ which(new[ ,24]== min(new[,24])),]
sample<-sample[1,2:21]
colnames(sample)<-clean8[1:20,1]
x<-as.vector(t(sample))

d1=NULL
for (i in seq(1,10000,1)) {
  y=as.vector(t(RandVec(a=-0.002, b=0.002, s=0, n=20, m=1)$RandVecOutput))
  x1=x+y
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  d1 = rbind(d1, data.frame(i,test7,pred))
  if (i==1) {
    if (pred >= 0.1) {
      x=x1
    } else {
      x=x
    } 
  } else {
    if (pred >= d1[i-1,22]) {
      x=x1
    } else {
      x=x
    }
  }
}
# max 0.8145091, just 10 seconds, very efficient
max(d1[,22])


# plot iteration and solution
plot(x = d1$i,y = d1$pred,
     xlab = "Iteration",
     ylab = "Objective function",
     xlim = c(1,10000),
     ylim = c(min(d1$pred),max(d1$pred)),		 
     main = "Objective function VS. iteration"
)
# OPS IIII-we need control the distance between varible with the original vector, set it as penalty
#initial value for x
rank=read.csv("rank_2.csv")
new <- rank[ which(rank[ ,22]== 0.10),]
sample <- new[ which(new[ ,24]== min(new[,24])),]
sample<-sample[1,2:21]
colnames(sample)<-clean8[1:20,1]
x0<-as.vector(t(sample))

d3=NULL
x=x0
for (i in seq(1,10000,1)) {
  Seed=sample(1:10000, size = 1)
  y=as.vector(t(RandVec(a=-0.01, b=0.01, s=0, n=20, m=1, Seed=Seed)$RandVecOutput))
  x1=x+y
  mutation=sum(abs(y))*350
  distance = dist(rbind(x0, x1), method = "euclidean")
  distance=as.numeric(distance)
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  result=pred-distance/0.1073986
  d3 = rbind(d3, data.frame(i,Seed,test7,pred,distance,result,mutation))
  if (i==1) {
    if (pred >= 0.1) {
      x=x1
    } else {
      x=x
    } 
  } else {
    if (result >= d3[i-1,25] & all(x1>=0) & all(x1<=1)) {
      x=x1
      solution=result
    } else {
      x=x
    }
  }
}
#solution 0.833 at this case, the max pred is 1, it is not resonable, scale penalty
max(d3$result)
write.csv(d3,"optimization1.csv")
#plot pred and diatance
d5=read.csv("optimization1.csv")
d5[,1]=NULL
new1<- d5[which(d5[ ,3]>=0 & d5[ ,3]<=1),]
for (i in seq(4,22,1)){
  new1<- new1[which(new1[ ,i]>=0 & new1[ ,i]<=1),]
}
new=new1[order(new1$pred,decreasing = FALSE),]
for (i in seq(1,nrow(new),1)) {
  x1=as.vector(new[i,3:22])
  y=x1-x0
  mutation=sum(abs(y))*350
  new[i,26]=mutation
}
colnames(new)[26]<-"mutation"
plot(x = new$pred,y = new$mutation,
     xlab = "Predicted solubility",
     ylab = "Number of mutated amino acids",
     xlim = c(min(new$pred),max(new$pred)),
     ylim = c(min(new$mutation),max(new$mutation)),		 
     main = "Predicted solubility VS. mutation"
)

#d2-result
#result=pred-distance/0.1073986
#after multiply step size, the distance may be large

d4=NULL 
x=x0
for (i in seq(1,100000,1)) {
  Seed=sample(1:100000, size = 1)
  y=as.vector(t(RandVec(a=-0.01, b=0.01, s=0, n=20, m=1, Seed=Seed)$RandVecOutput))
  x1=x+y
  distance = dist(rbind(x0, x1), method = "euclidean")
  distance = as.numeric(distance)
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  d4 = rbind(d4, data.frame(i,Seed,test7,pred,distance))
  if (i==1) {
    if (pred >= 0.1) {
      x=x1
    } else {
      x=x
    } 
  } else {
    if (pred >= d4[i-1,22] & distance <=0.05 & all(x1>=0) & all(x1<=1)) {
      x=x1
      solution=pred
    } else {
      x=x
    }
  }
}
write.csv(d4,"optimization2.csv")

# we need screen the table because pred unsatified conditions are also record,solution is right
# b=0.01, distance=0.05, iteration=1000, solution=0.269,178,0.322
# b=0.005, distance=0.05, iteration=1000,solution=0.192,0.168,0.207
# b=0.01, distance=0.05, iteration=10000, solution=0.16
max(d4$pred)
plot(x = d4$i,y = d4$pred,
     xlab = "Iteration",
     ylab = "Objective function",
     xlim = c(1,10000),
     ylim = c(min(d4$pred),max(d4$pred)),		 
     main = "Objective function VS. iteration"
)
new.data <- new[ which(new[ ,24]<=0.1),]
#OPS II

d5=NULL 
x=x0
for (i in seq(1,10000,1)) {
  Seed=sample(1:10000, size = 1)
  y=as.vector(t(RandVec(a=-0.01, b=0.01, s=0, n=20, m=1, Seed=Seed)$RandVecOutput))
  x1=x+y
  distance = dist(rbind(x0, x1), method = "euclidean")
  distance = as.numeric(distance)
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  d5 = rbind(d5, data.frame(i,Seed,test7,pred,distance))
  if (i==1) {
    if (pred >= 0.1) {
      x=x1
    } else {
      x=x
    } 
  } else {
    if (pred >= d5[i-1,22] & all(x1>=0) & all(x1<=1)) {
      x=x1
      solution=pred
    } else {
      x=x
    }
  }
}
write.csv(d5,"optimization3.csv")
#solution 0.76, screen solution 0.278 so d4 and d5 are similar
d5=d4
new1<- d5[which(d5[ ,3]>=0 & d5[ ,3]<=1),]
for (i in seq(4,22,1)){
  new1<- new1[which(new1[ ,i]>=0 & new1[ ,i]<=1),]
}
new=new1[order(new1$pred,decreasing = TRUE),]
new2<- new[which(new$distance<=0.05),]
x1=as.vector(new2[1,3:22])
y=x1-x0
mutation=sum(abs(y))*350
write.csv(d4,"optimization4.csv")
write.csv(new2,"optimization4filter.csv")
#mutation 63.48943, pred 0.6563880, distance 0.04845306, run 40000, 
sum(x1)
#record number of mutated amino acid
mutation=sum(abs(y))*350

abs(y)
sum(abs(y))
distance = dist(rbind(LB, UB), method = "euclidean")

#select 10 proteins we want to mutate
rank=read.csv("rank_2.csv")
#there are 58 proteins with real solubility 0.1
new <- rank[ which(rank[ ,22]== 0.10),]
#whether predicted solubility is close with real-rank them
new1=new[order(new$difference,decreasing = FALSE),]
#find sequence of record length
total=read.csv("solubility2.csv")
for (i in seq(1,58,1)){
  sample1 <- total[ which(total[ ,1]== new1[i,1]),]
  sequence=sample1[1,4]
  sequence<-as.character(sequence)
  new1[i,25]=nchar(sequence)
}
# !!! here it should be sample1 <- total[ which(as.numeric(rownames(total[2980,]))== new1[i,1]),]
total=read.csv("solubility2.csv")
for (i in seq(1,58,1)){
  sample1 <- total[ which(as.numeric(rownames(total))== new1[i,1]),]
  sequence=sample1[1,4]
  sequence<-as.character(sequence)
  new1[i,25]=nchar(sequence)
}
#fileter length less than 1kb, 27 proteins, choose 10 according to difference
sample <- new1[ which(new1[ ,25]<=333),]
sample1=sample[1:10,]
#write.csv(sample1,"optimization_sample.csv") 
write.csv(sample1,"new_optimization_sample.csv") 
new3=new2[order(new2$distance,decreasing = FALSE),]
#set a parameter as length of protein after additon, 20 varibles and 20 equations, give adding length x>=0
length*x0+A=(length+addition)*x1
A=(length+addition)*x1-length*x0=(350+addition)*x1[1,1]-350*x0[1]
for (i in seq(1,20,1)){
  x1[3,i]=(350+587.686)*x1[1,i]-350*x0[i]
}

#if addition =0, the addition is the less one, we want additon >0
#calculate the min amino acid A>=0
for (i in seq(1,20,1)){
  x1[2,i]=350*x0[i]/x1[1,i]-350 
}
min_addition=max(x1[2,1:20])
#587.686
#mutation 

# allow both addition and mutation to find the minimum change 
new2=read.csv("optimization4filter.csv")
x2=as.vector(new2[1,4:23])
x2=as.vector(t(x2))

library("GA")
fun <- function(y) {
  mutation=sum(abs(x2*(350+y)-x0*350))
  return(-mutation)
}
GA <- ga(type = "real-valued", fitness = fun, lower = 0, upper = 100, popSize = 100, maxiter = 1000)
summary(GA)
#solution 0.0103415  objective -63.49042
# definition of mutation is different (abs) with number of addition for sequence
# deletion is hard to differentiate with mutation, when length is larger than previous one, it is mutation not deletion
# try to change the distance threshold 0.05->0.02
new<- new2[which(new2$distance<=0.03),]
# no results for 0.02, we need to change the threshold in coding
# choose the first one for 0.03
x2=as.vector(new[1,4:23])
x2=as.vector(t(x2))

library("GA")
fun <- function(y) {
  mutation=sum(abs(x2*(350+y)-x0*350))
  return(-mutation)
}
GA <- ga(type = "real-valued", fitness = fun, lower = 0, upper = 100, popSize = 100, maxiter = 1000)
summary(GA)
# solution 3.749983 objective -30.22575 solubility 0.277199114  distance 0.02558098

# definition of mutation is different (abs) with number of addition for sequence
# deletion is hard to differentiate with mutation, when length is larger than previous one, it is mutation not deletion
# 63.49083

# try my algorithm for 10 proteins
# use all data to try model
clean8=read.csv("composition_difference.csv")
clean7=read.csv("cleandata_1.csv")
clean7[,1]<-NULL
colnames(clean7)<-clean8[1:20,1]
colnames(clean7)[21]<-"solubility"
library("e1071")
svm_model <- svm(solubility ~ ., clean7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = 0.057)
#initial value for x
#sample=read.csv("optimization_sample.csv")
sample=read.csv("new_optimization_sample.csv")
sample<-sample[1:10,3:22]
colnames(sample)<-clean8[1:20,1]
#write.csv(sample,"sample_10.csv")
write.csv(sample,"new_sample_10.csv")
#It seems pred is column 23 not 22, check for protein 622 too!   
sample=read.csv("new_sample_10.csv")
sample[,1]=NULL
x0=sample[6,1:20]
x0=as.vector(t(x0))
d1=NULL 
x=x0
for (i in seq(1,100000,1)) {
  Seed=sample(1:100000, size = 1)
  y=as.vector(t(RandVec(a=-0.01, b=0.01, s=0, n=20, m=1, Seed=Seed)$RandVecOutput))
  x1=x+y
  distance = dist(rbind(x0, x1), method = "euclidean")
  distance = as.numeric(distance)
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  d1 = rbind(d1, data.frame(i,Seed,test7,pred,distance))
  if (i==1) {
    if (pred >= 0.1) {
      x=x1
    } else {
      x=x
    } 
  } else {
    if (pred >= d1[i-1,23] & distance <=0.05 & all(x1>=0) & all(x1<=1)) {
      x=x1
      solution=pred
    } else {
      x=x
    }
  }
}
write.csv(d1,"protein1.csv") 
#1 50000 iteration solution 0.471(pred23) #12:40 start  22:00 stop  100000 iteration solution 0.477
#2 100000 iteration  0.42
#3 100000 iteration  0.552
#4 100000 0.353  iteration record time for analysis
#5 100000 0.485
#6 100000 0.539 iteration record time for analysis
# protein7  0.618
# protein8  0.615
# protein9  0.349
# protein10  0.778

#2-10 protein
x0=sample[6,1:20]
x0=as.vector(t(x0))
d6=NULL 
x=x0
for (i in seq(1,100000,1)) {
  start_time <- Sys.time()
  Seed=sample(1:100000, size = 1)
  y=as.vector(t(RandVec(a=-0.01, b=0.01, s=0, n=20, m=1, Seed=Seed)$RandVecOutput))
  x1=x+y
  distance = dist(rbind(x0, x1), method = "euclidean")
  distance = as.numeric(distance)
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  if (i==1) {
    if (pred >= 0.1) {
      x=x1
    } else {
      x=x
    } 
  } else {
    if (pred >= d6[i-1,23] & distance <=0.05 & all(x1>=0) & all(x1<=1)) {
      x=x1
      solution=pred
    } else {
      x=x
    }
  }
  end_time <- Sys.time()
  t=end_time - start_time
  time=as.numeric(t,units="secs")
  d6 = rbind(d6, data.frame(i,Seed,test7,pred,distance,time))
}
write.csv(d6,"protein6.csv") 

#plot solution for different iterations
d4=read.csv("optimization4.csv")
plot(x = d5$i,y = d5$time,
     xlab = "Iteration",
     ylab = "Objective function",
     pch=".",
     #xlim = c(1,41781),
     #ylim = c(min(d4$pred),max(d4$pred)),  
     main = "Objective function VS. iteration"
)
#plot time vs iteration for d4
d4=read.csv("protein4.csv")
d5=read.csv("protein6.csv")
plot(x = d5$i,y = d5$time,
     xlab = "Iteration",
     ylab = "Time",
     #pch=".",
     #xlim = c(1,41781),
     #ylim = c(min(d4$pred),max(d4$pred)),  
     main = "Time VS. iteration"
)

d5[1,27]=d5[1,26]
for (i in seq(2,100000,1)){
  d5[i,27]=d5[i,26]+d5[i-1,27]
}
colnames(d5)[27]<-"Accumulate_time"

plot(x = d5$i,y = d5$Accumulate_time,
     xlab = "Iteration",
     ylab = "Accumulate_Time(secs)",
     #pch=".",
     main = "Accumulate_Time VS. iteration"
)

#add 30 amino acids
d1=NULL 
for (i in seq(1,100,1)) {
  Seed=sample(1:100, size = 1)
  y=as.vector(t(RandVec(a=0, b=30, s=30, n=20, m=1, Seed=Seed)$RandVecOutput))
  x1=(x0*350+y)/380
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  amino<- data.frame(t(y))
  colnames(amino)<-clean8[1:20,1]
  d1 = rbind(d1, data.frame(i,Seed,amino,pred))
}
max_solubility=max(d1$pred)
View(d1)
max_solubility
#0.1431891

d1=NULL 
for (i in seq(1,100000,1)) {
  Seed=sample(1:100000, size = 1)
  y=as.vector(t(RandVec(a=0, b=30, s=30, n=20, m=1, Seed=Seed)$RandVecOutput))
  x1=(x0*350+y)/380
  test7 <- data.frame(t(x1)) 
  colnames(test7)<-clean8[1:20,1]
  pred <- predict(svm_model, newdata = test7)
  amino<- data.frame(t(y))
  colnames(amino)<-clean8[1:20,1]
  d1 = rbind(d1, data.frame(i,Seed,amino,pred))
}
write.csv(d1,"addtion.csv")
max(d1$pred)
#0.2084518

#find minimum mutated amino acids for 10 proteins
solution=NULL
for (i in seq(6,6,1)) {
  n=paste("protein",i,".csv",sep = "", collapse = "")
  d5=read.csv(n)
  d5[,1]<-NULL
  for (x in seq(3,22,1)){
    d5<- d5[which(d5[ ,x]>=0 & d5[ ,x]<=1),]
  }
  new=d5[order(d5$pred,decreasing = TRUE),]
  new2<- new[which(new$distance<=0.05),]
  x1=new2[1,1:24]
  solution = rbind(solution, data.frame(i,x1))
}
write.csv(solution,"solution_10.csv") 
qq=read.csv("solution_10.csv")
qq[,1]=NULL
qq[6,1:25]=solution[1,1:25]
write.csv(qq,"solution_10.csv") 
# For the solution recorded, it may not be the highest resonale solution because if, right 0.7, wrong 0.8, right 0.75, 
# it is not higher the last row but it is higher than the last right solution. So this screen result is more reasonable!
optimize=read.csv("new_optimization_sample.csv")
solution=read.csv("solution_10.csv")
sample=read.csv("new_sample_10.csv")
mutation=NULL
for (i in seq(6,6,1)){
  x2=as.vector(solution[i,5:24])
  x2=as.vector(t(x2))
  x0=as.vector(sample[i,2:21])
  x0=as.vector(t(x0))
  fun <- function(y) {
    length=optimize[i,26]
    mutation=sum(abs(x2*(length+y)-x0*length))
    return(-mutation)
  }
  GA <- ga(type = "real-valued", fitness = fun, lower = 0, upper = 100, popSize = 100, maxiter = 1000)
  objective=GA@fitnessValue
  length=GA@solution
  mutation=rbind(mutation, data.frame(i,length,objective))
}
write.csv(mutation,"mutation_10.csv")
original=read.csv("mutation_10.csv")
original[,1]=NULL
original[6,1:3]=mutation[1,1:3]
write.csv(original,"mutation_10.csv")
#length of protein changed for 10 proteins

#re run for sample 6, 9, 10
optimize=read.csv("new_optimization_sample.csv")
solution=read.csv("solution_10.csv")
sample=read.csv("sample_10.csv")
mutation=NULL
for (i in seq(1,10,1)){
  x2=as.vector(solution[i,5:24])
  x2=as.vector(t(x2))
  x0=as.vector(sample[i,2:21])
  x0=as.vector(t(x0))
  fun <- function(y) {
    length=optimize[i,26]
    mutation=sum(abs(x2*(length+y)-x0*length))
    return(-mutation)
  }
  GA <- ga(type = "real-valued", fitness = fun, lower = 0, upper = 100, popSize = 100, maxiter = 1000)
  objective=GA@fitnessValue
  length=GA@solution
  mutation=rbind(mutation, data.frame(i,length,objective))
}
write.csv(mutation,"mutation_10.csv")


# find the protein function
optimize=read.csv("optimization_sample.csv")
total=read.csv("solubility2.csv")

#another 6 samples
library("readxl")
traningdata <- read_excel("6_sample.xlsx")
x=data.frame(traningdata)
length(x[ ,2])
x[ ,2]<-as.character(x[ ,2])
library("protr")
x[ ,2] = x[ ,2][(sapply(x[ ,2], protcheck))]
length(x[ ,2])
x1 = t(sapply(x[ ,2], extractAAC))
x[1:6,3]=nchar(x[1:6,2])
x[ ,4:23] <- x1[ ,1:20]
write.csv(x,"sample_6.csv")
#
clean8=read.csv("composition_difference.csv")
clean7=read.csv("cleandata_1.csv")
clean7[,1]<-NULL
colnames(clean7)<-clean8[1:20,1]
colnames(clean7)[21]<-"solubility"
new=x[,4:23]
colnames(new)<-clean8[1:20,1]
library("e1071")
svm_model <- svm(solubility ~ ., clean7, method="eps-regression", cost=1.26, gamma=0.029, epsilon = 0.057)
pred <- predict(svm_model, newdata = new)
prediction<-data.frame(pred)

x$solubility=prediction[,1]
write.csv(x,"sample_6.csv")
matlab_5=x[1:3,4:23]
matlab_5[4,1:20]=x[5,4:23]
matlab_5[5,1:20]=clean7[385,2:21]
write.csv(matlab_5,"matlab_5.csv")

new=read.csv("cleandata_1.csv")
fold2 = new[which(new$solubility<=0.7 & new$solubility > 0.3),]



# screen the 10 proteins for experiments

