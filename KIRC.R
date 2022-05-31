path<-"D:/learning/KIRC"
source("D:/learning/methods/rfunction/mer.clin/mer.clin.R")
clinical<-mer.clin(path)
clinical<-as.data.frame(clinical)
daysToDeath <- clinical$days_to_death
nonComplt <- is.na(daysToDeath)
vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))


sumNA<-function(x){
  sum(is.na(x))
}

apply(Methylation,2,sumNA)
setwd("D:/learning/KIRC/methana")
#bladder methylation
rm(list=ls())

Methylation<-read.table("D:/learning/KIRC/TCGA-KIRC.methylation450.tsv/TCGA-KIRC.methylation450.tsv",sep="\t",header=T,row.names = 1,check.names = F)
Methylation[1:7,1:7]
source("D:/learning/methods/rfunction/rmNA/rmNA.R.txt")
Methylation<-rmNA(Methylation,0.9)
Methylation<-Methylation[,substr(colnames(Methylation), 14, 16) %in% c("01A","11A")]

###################################################################################################################
#different methylation analysis
###################################################################################################################
sample_id<-substr(colnames(Methylation), 14, 15)
#####鐢熸垚杩囨护浣嶇偣鐨刾d鏂囦???
Group=ifelse(sample_id=="01","C","T")
N=length(colnames(Methylation))
Array.matrix=combn(c(paste(0,1:9,sep=""),as.character(11:50)),2)
Array1=c(Array.matrix[1,],Array.matrix[2,])
Array2=c(Array.matrix[2,],Array.matrix[1,])
pd.data<-data.frame(Sample_Name=colnames(Methylation),Sample_Plate=rep(NA,N),
                    Sample_Group=Group,Pool_ID=rep(NA,N),Project=rep(NA,N),
                    Sample_Well=rep("E01",N),Slide=seq(7990895118,7990895117+N,1),
                    Array=paste(rep("R",N),Array1[1:N],rep("C",N),Array2[1:N],sep=""),stringsAsFactors=F)

library(ChAMP)
###
my.filter<-champ.filter(beta=Methylation,pd=pd.data,fixOutlier=F)
dim(my.filter)
Methylation<-my.filter$beta

Methylation<-rmNA(Methylation,0.9)
dim(Methylation)
################################################################鐢插熀鍖栫煩闃典笌涓村簥淇℃伅鏁寸???
####################################################################################################
#differentially analysis +UNICOX
#Differentially analysis

patient_id<-substr(colnames(Methylation), 1, 12)
paired<-which(duplicated(patient_id)==T)
paired_id<-colnames(Methylation)[patient_id %in% patient_id[paired]]
paired_id<-c(paired_id[substr(paired_id, 14, 15)=="11"],paired_id[substr(paired_id, 14, 15)=="01"])
Group<-ifelse(substr(paired_id,14,15)=="01","T","C")
paried.methylation<-Methylation[,paired_id]
tumor.methylation<-Methylation[,sample_id=="01"]
save(paried.methylation,Group,file="paired.meth-Group.Rdata")
#load("paired.meth-Group.Rdata")
####################################################################################################
#宸紓鐢插熀鍖栦綅鐐?
adj.pvalue=0.9999999999
cut.deltaBeta=0.3
myDMP<-champ.DMP(beta=paried.methylation,pheno=factor(Group),adjPVal=adj.pvalue)
all.methylation<-myDMP[[1]]


adj.pvalue=0.05
cut.deltaBeta=0.3
myDMP<-champ.DMP(beta=paried.methylation,pheno=factor(Group),adjPVal=adj.pvalue)

alldiff.methylation<-myDMP[[1]]
Diff.methylation<-subset(alldiff.methylation,abs(deltaBeta)>cut.deltaBeta)
dim(Diff.methylation)
dir.create("result")
write.csv(Diff.methylation,"./result/Diff.methylation.csv")
index<-intersect(row.names(Methylation),row.names(Diff.methylation))
Methylation.diff<-Methylation[index,]
write.csv(Methylation.diff,"./result/methylation.diff.csv")
alldiff.methylation<-alldiff.methylation[complete.cases(alldiff.methylation$deltaBeta),]
#methylation.diff<-read.csv("methylation.diff.csv",header=T,row.names=1)
##########################################################################################
#volcano plot
#rm(list=ls())
setwd("D:/learning/KIRC/methana")
library(ggplot2)
library(ggrepel)
all.methylation<-all.methylation[complete.cases(all.methylation$deltaBeta),]
for (i in 1:nrow(all.methylation)) {
  if ((abs(all.methylation[i,'deltaBeta']) >= 0.3)|all.methylation[i,'deltaBeta'] %in% NA) all.methylation[i,'select_change'] <- 'y' else all.methylation[i,'select_change'] <- 'n'
  if (abs(all.methylation[i,'adj.P.Val']) >= 0.05) all.methylation[i,'select_pvalue'] <- 'n' else all.methylation[i,'select_pvalue'] <- 'y'
  all.methylation[i,'select'] <- paste(all.methylation[i,'select_change'], all.methylation[i,'select_pvalue'], sep = '')
}
all.methylation$select <- factor(all.methylation$select, levels = c('nn', 'ny', 'yn', 'yy'), labels = c('Adjusted P >= 0.05, Delta Beta < 0.3', 'Adjusted P < 0.05, Delta Beta < 0.3', 'Adjusted P >= 0.05, Delta Beta >= 0.3', 'Adjusted P < 0.05, Delta Beta >= 0.3'))
#all.methylation$select <- factor(all.methylation$select, levels = c('nn', 'ny', 'yn', 'yy'), labels = c('p >= 0.05, FC < 2', 'p < 0.05, FC < 2', 'p >= 0.05, FC >= 2', 'p < 0.05, FC >= 2'))
#????Ϊ?????? p ֵ
volcano_plot_pvalue <- ggplot(all.methylation, aes(deltaBeta, -log(adj.P.Val, 10))) +
  geom_point(aes(color = select), alpha = 0.6) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.3, 0.3), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  labs(x = 'Delta Beta', y = '-log10 adjusted P value')+
  theme(legend.position=c(0.8,0.9))
save(all.methylation,volcano_plot_pvalue,file="volcano_plot_pvalue.Rdata")
pdf("./result/vocanoplot.pdf")
volcano_plot_pvalue
dev.off()
####################################################################################################
#
rm(list=ls())
clinical1<-read.table("D:/learning/KIRC/KIRC_clinicalMatrix/KIRC_clinicalMatrix",header=T,sep="\t",row.names=1,check.names=F)
clinical2<-read.csv("D:/learning/KIRC/TCGA-KIRC.GDC_phenotype/TCGA-KIRC.GDC_phenotype.tsv",header=T,sep="\t",row.names=1,check.names=F)
clinical3<-read.csv("D:/learning/KIRC/TCGA-KIRC.survival.tsv/TCGA-KIRC.survival.tsv",header=T,sep="\t",row.names=1,check.names=F)
head(row.names(clinical1))
head(row.names(clinical2))
head(row.names(clinical3))
keepsam<-function(x,location,mathch){
  x<-x[substr(row.names(x),14,location) %in% mathch,]
  if(location==16){
    row.names(x)<-substr(row.names(x),1,15)
  }
  return(x)
}

clinical2<-keepsam(clinical2,16, "01A")
clinical3<-keepsam(clinical3,16,"01A")
clinical1<-keepsam(clinical1,15,"01")

index<-intersect(intersect(row.names(clinical1),row.names(clinical2)),row.names(clinical3))
clinical1<-clinical1[index,]
clinical2<-clinical2[index,]
clinical3<-clinical3[index,]

clinical<-clinical2
clinical$OS.time<-clinical3$`_OS`
clinical$OS<-clinical3$`_OS_IND`
clinical$RFS.time<-clinical1$RFS.time
clinical$RFS<-clinical1$RFS
remove(clinical1,clinical2,clinical3)
row.names(clinical) <- gsub("-",".",row.names(clinical))
clinical<-clinical[complete.cases(clinical$OS.time)&clinical$OS.time!=0,]
clinical<-subset(clinical,select=c(OS.time,OS,RFS.time,RFS,additional_pharmaceutical_therapy,additional_radiation_therapy,
                                   additional_surgery_locoregional_procedure,additional_surgery_metastatic_procedure,age_at_initial_pathologic_diagnosis,
                                   eastern_cancer_oncology_group,karnofsky_performance_score,laterality,hemoglobin_result,neoplasm_histologic_grade,
                                   number_pack_years_smoked,pathologic_M,pathologic_N,pathologic_T,radiation_therapy,radiation_treatment_ongoing,radiation_type,
                                   stopped_smoking_year,gender.demographic,serum_calcium_result,race.demographic,morphology.diagnoses,tumor_grade.diagnoses,
                                   tumor_stage.diagnoses))
write.csv(clinical,"./result/clinical.csv")
####################################################################################################
setwd("D:/learning/KIRC/methana")
rm(list=ls())
Methylation.diff<-read.csv("./result/Methylation.diff.csv",header=T,row.names=1)
clinical<-read.csv("./result/clinical.csv",header=T,row.names=1)
names(Methylation.diff) <- gsub("-",".",names(Methylation.diff))
names(clinical)

clinical$tumor_stage.diagnoses<-as.character(clinical$tumor_stage.diagnoses)
clinical$tumor_stage.diagnoses[clinical$tumor_stage.diagnoses == "stage i"] <- 1
clinical$tumor_stage.diagnoses[clinical$tumor_stage.diagnoses == "stage ii"] <- 2
clinical$tumor_stage.diagnoses[clinical$tumor_stage.diagnoses == "stage iii"] <- 3
clinical$tumor_stage.diagnoses[clinical$tumor_stage.diagnoses == "stage iv"] <- 4
clinical$tumor_stage.diagnoses[clinical$tumor_stage.diagnoses == "not reported"] <- NA
clinical$tumor_stage.diagnoses<-as.numeric(clinical$tumor_stage.diagnoses)
clinical$neoplasm_histologic_grade<-as.character(clinical$neoplasm_histologic_grade)
clinical$neoplasm_histologic_grade[clinical$neoplasm_histologic_grade == "G1"] <- 1
clinical$neoplasm_histologic_grade[clinical$neoplasm_histologic_grade == "G2"] <- 2
clinical$neoplasm_histologic_grade[clinical$neoplasm_histologic_grade == "G3"] <- 3
clinical$neoplasm_histologic_grade[clinical$neoplasm_histologic_grade == "G4"] <- 4
clinical$neoplasm_histologic_grade[clinical$neoplasm_histologic_grade == ""] <- NA
clinical$neoplasm_histologic_grade<-as.numeric(clinical$neoplasm_histologic_grade)
names(Methylation.diff)<-substr(names(Methylation.diff),1,15)
index<-intersect(row.names(clinical),names(Methylation.diff))
clinical<-clinical[index,]
expres<-Methylation.diff[,index]
dim(Methylation.diff)
dim(clinical)
names(clinical)
################################################################################
#Uinvariate Cox analysis
################################################################################
source("D:/learning/methods/rfunction/unicox/unicoxph.R")
expres<-as.data.frame(expres)
unicox<-unicox(clinical$OS.time,clinical$OS,expres)
unicox<-as.data.frame(unicox)
#FDR
unicox$fdr<-p.adjust(unicox$pValue,method = "bonferroni")
# annotation
load("D:/learning/methods/rfunction/Annotation.450K.Rdata")
index<-intersect(row.names(annotation),row.names(unicox))
unicox$genesymbol<-annotation[row.names(unicox),]$gene
unicox<-unicox[unicox$genesymbol!="",]
dim(unicox)
write.csv(unicox,"./result/unicox.os.csv")
sigunicox<-unicox[which(unicox$fdr<0.01),]
#sigunicox<-unicox[which(unicox$pValue<0.00001),]
write.csv(sigunicox,"./result/sigunicox.csv")
can<-expres[row.names(sigunicox),]
write.csv(can,"./esult/can.csv")
write.csv(clinical,"./result/clinical.csv")
#####################################################################################################
#elastic net
#####################################################################################################
#####################################################################################################
rm(list=ls())
setwd("D:/learning/KIRC/methana")
can<-read.csv("./result/can.csv",header=T,row.names=1)
clinical<-read.csv("./result/clinical.csv",header=T,row.names=1)
can<-na.omit(can)
library(survivalROC)
library(ggplot2)
library(survminer)
library(glmnet)
library(c060)
library(caret)
seed=2357
#####################################################################################################
while(seed>0){
  set.seed(seed+767) #4513
  train.idx <- createDataPartition(clinical$OS, p = 0.6, list = FALSE)
  
  test.surv <- clinical[-train.idx,]
  test <- can[,-train.idx]
  train.surv <- clinical[train.idx,]
  y.train <- train.surv[,1:2]
  y.test<-test.surv[,1:2]
  
  train <- can[,train.idx]
  names(y.train) <- c("time","status")
  y.train<-as.matrix(y.train)
  train<-as.matrix(t(train))
  names(y.test)<-c("time","status")
  test<-as.matrix(t(test))
  bounds <- t(data.frame(alpha = c(0, 1)))
  colnames(bounds) <- c("lower", "upper")
  nfolds <- 10
  foldid <- balancedFolds(class.column.factor = y.train[, 2],cross.outer = nfolds)
  fitEN <- epsgo(Q.func = "tune.glmnet.interval", bounds = bounds,
                 parms.coding = "none", seed = seed, fminlower = -100, x = train, y = y.train,
                 family = "cox", foldid = foldid, type.min = "lambda.1se",
                 type.measure = "deviance")
  sumint <- summary(fitEN, verbose = TRUE)
  fit <- glmnet(x = train, y = y.train,family = "cox",standardize=F,lambda = sumint$opt.models$model$cvreg$lambda.1se, alpha=sumint$opt.alpha)
  risk.test <- predict(fit,test,type = "response")
  surv.test<- test.surv
  surv.test$risk <- risk.test
  years<-seq(from=0,to=max(surv.test$OS.time),by=365)
  years<-years[-1]
  sumROC<-list()
  for (i in 1:length(years)){
    sumROC[[i]] <- survivalROC(Stime = surv.test$OS.time,status = surv.test$OS,marker = surv.test$risk,
                               predict.time =years[i],method = "NNE",span = 0.25*nrow(surv.test)^(-0.20))
  }
  sumAUC<-list()
  for (i in 1:length(sumROC)){
    sumAUC[[i]]<-sumROC[[i]]$AUC
  }
  minAUC<-min(sumAUC[[1]],sumAUC[[3]],sumAUC[[5]],sumAUC[[7]])
  coef<- coef(fit, s = fit$lambda.lse)
  coef<-as.matrix(coef)
  coef<-as.data.frame(coef)
  coef<-subset(coef,s0!=0)
  if(minAUC>0.65 & nrow(coef)<18 & nrow(coef)>3){
    break
  }else {
    seed<-seed+785
  }
}
coef
#seed5497+767
########################################################################################################
save(coef,fit,test.surv,train,train.surv,test,y.test,y.train,file="./result/seed5497.Rdata")
#load("./result/seed5497.Rdata")
write.csv(test.surv,"./result/test.surv.csv")
write.csv(train.surv,"./result/train.surv.csv")
write.csv(test,"./result/test.csv")
write.csv(train,"./result/train.csv")
save(fit,file="./result/fit.Rdata")
save(fitEN,file="./result/fitEN.Rdata")
load("./result/fitEN.Rdata")
sumint <- summary(fitEN, verbose = TRUE)
sumint$opt.alpha
sumint$opt.lambda
#Optimal parameters found are: 
#alpha =  0.9932 	 lambda =  0.1187 deviance =  15.2505
pdf("./result/tuning parameters alpha and loglamda.pdf")
plot(sumint)
dev.off()
pdf("./result/The distribution of initial and visited points of the interval search plotted in chronological order.pdf")
plot(sumint, type = "points")
dev.off()
fit <- glmnet(x = train, y = y.train,family = "cox",standardize=F,lambda = sumint$opt.models$model$cvreg$lambda.1se, alpha=sumint$opt.alpha)
library(BootValidation)
#vboot<-vboot(fit, x = train, y = y.train, s=sumint$opt.models$model$cvreg$lambda.1se, nfolds = 5, B = 200,  cv_replicates = 100, lambda = TRUE, n_cores = max(1,parallel::detectCores() - 1))
#save(vboot,file="./result/vboot.Rdata")
#pdf("./result/vboot.pdf",width=32,height = 28)
#plot(vboot)
#dev.off()

#load("./result/seed4513.Rdata")
coef<- coef(fit, s = fit$lambda.lse)
coef<-as.matrix(coef)
coef<-as.data.frame(coef)
coef<-subset(coef,s0!=0)
coef
write.csv(coef,"./result/coef.csv")
#draw the heatmap for coef gene
load("paired.meth-Group.Rdata")
coef<-read.csv("./result/coef.csv",header=T,row.names=1)
library(ComplexHeatmap)
can.met<-paried.methylation[row.names(coef),]
library(ComplexHeatmap)
#heatdata<-as.data.frame(t(normalize(t(heatdata))))
pdf("./result/difmet.pdf",width=15,height=5)
Heatmap(as.matrix(can.met),border = F,column_title =c ("Normal renal tissue","Renal cell carcinoma"),
        column_title_gp = gpar(fill = c("red","blue"), col = "black",fontsize = 22),cluster_rows = F,
        #row_title =c (paste(row.names(coef))),
        cluster_columns = T,name = "", show_column_names = F,column_split = Group)
dev.off()
metpcr<-read.csv("metpcr.csv",header = T,row.names=1)
Group<-substr(names(metpcr),3,3)
pdf("./result/metpcr.pdf",width=15,height=5)
Heatmap(
  as.matrix(metpcr),
  border = F,
  #column_title = c ("Normal renal tissue", "Renal cell carcinoma"),
  #column_title_gp = gpar(
  #  fill = c("red", "blue"),
  #  col = "black",
  #  fontsize = 22
  #),
  cluster_rows = F,
  #row_title =c (paste(row.names(coef))),
  cluster_columns = F,
  name = "",
  show_column_names = F
  #column_split = Group
)
dev.off()
Group


ls()
rm(list=ls())
load("./result/seed5497.Rdata")
ls()
risk.train<- predict(fit,train,type = "response")
surv.train<- train.surv
surv.train$risk <- risk.train
write.csv(surv.train, file="./result/surv.train.csv")
##caclulate risk score in the test
head(row.names(test)) #row sample

risk.test <- predict(fit,test,type = "response")
dim(risk.test)
surv.test<- test.surv
dim(test.surv)
surv.test$risk <- risk.test
write.csv(surv.test, file="./result/surv.test.csv")
#get risk score for the validation set
#GSE44001<-read.csv("D:/learning/cervicalcancer/data/GSE44001_series_matrix.txt/GSE44001.csv",header=T,row.names=1)
#GSE44001clinical<-read.csv("D:/learning/cervicalcancer/data/GSE44001_series_matrix.txt/GSE44001clinical.csv",header=T,row.names=1)
#GSE44001<-GSE44001[,order(names(GSE44001))]
#GSE44001clinical<-GSE44001clinical[order(row.names(GSE44001clinical)),]
####################################################################################
#base risk score
setwd("D:/learning/KIRC/methana")
#bladder methylation
rm(list=ls())
Methylation<-read.table("D:/learning/KIRC/TCGA-KIRC.methylation450.tsv/TCGA-KIRC.methylation450.tsv",sep="\t",header=T,row.names = 1,check.names = F)
Methylation<-Methylation[,substr(colnames(Methylation), 14, 16) %in% c("01A","11A")]
coef<-read.csv("./result/coef.csv",header=T,row.names=1)
index<-intersect(row.names(coef),row.names(Methylation))
methylation<-Methylation[index,]

train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
#survdata<-survdata[order(survdata$risk),]
names(methylation)<-gsub("-",".",names(methylation))
names(methylation)<-substr(names(methylation),1,15)

index<-intersect(names(methylation),rownames(survdata))
methylation.risk<-methylation[,index]
survdata<-survdata[index,]

methylation.normal<-methylation[,substr(names(methylation),14,15) %in% "11"]
methylation.low.risk<-methylation.risk[,which(survdata$Group %in% "Low risk")]
methylation.high.risk<-methylation.risk[,which(survdata$Group %in% "High risk")]
methylation<-cbind(methylation.normal,methylation.low.risk,methylation.high.risk)
group<-c()
for(i in 1:length(names(methylation))){
  if(names(methylation)[i] %in% names(methylation.normal)){
    group[i]<-"Normal renal tissue"
  } else if (names(methylation)[i] %in% names(methylation.low.risk)){
    group[i]<-"Low risk"
  } else {
    group[i]<-"High risk"
  }
}

Group<-cbind(names(methylation),group)
Group<-as.data.frame(Group)
names(Group)<-c("sample","group")
row.names(Group)<-names(methylation)
library(ComplexHeatmap)
#heatdata<-as.data.frame(t(normalize(t(heatdata)))) #FC4E07","#00AFBB"
methylation<-as.matrix(methylation)
pdf("./result/heatmaprisk.pdf",width = 8,height=3.5)
Heatmap(methylation, cluster_rows = FALSE, cluster_columns = F,name = "",show_column_names = F,
        column_title = c (paste(levels(Group$group)[1], "group"), paste(levels(Group$group)[2], "group"),paste(levels(Group$group)[3], "group")),
        column_title_gp = gpar( fill = c("#FC4E07","#00AFBB","green" ),col = "black",fontsize = 12 ),
        column_split = Group$group
)
dev.off()



#global DNA picture

rm(list=ls())
#Methylation<-read.table("D:/learning/KIRC/TCGA-KIRC.methylation450.tsv/TCGA-KIRC.methylation450.tsv",sep="\t",header=T,row.names = 1,check.names = F)
Methylation<-Methylation[,substr(colnames(Methylation), 14, 16) %in% c("01A","11A")]
coef<-read.csv("./result/coef.csv",header=T,row.names=1)
index<-intersect(row.names(coef),row.names(Methylation))
methylation<-Methylation[index,]

train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
#survdata<-survdata[order(survdata$risk),]
names(Methylation)<-gsub("-",".",names(Methylation))
names(Methylation)<-substr(names(Methylation),1,15)

index<-intersect(names(Methylation),rownames(survdata))
Methylation.risk<-Methylation[,index]
survdata<-survdata[index,]

Methylation.normal<-Methylation[,substr(names(Methylation),14,15) %in% "11"]
Methylation.low.risk<-Methylation.risk[,which(survdata$Group %in% "Low risk")]
Methylation.high.risk<-Methylation.risk[,which(survdata$Group %in% "High risk")]
Methylation<-cbind(Methylation.normal,Methylation.low.risk,Methylation.high.risk)
group<-c()
for(i in 1:length(names(Methylation))){
  if(names(Methylation)[i] %in% names(Methylation.normal)){
    group[i]<-"Normal renal tissue"
  } else if (names(Methylation)[i] %in% names(Methylation.low.risk)){
    group[i]<-"Low risk"
  } else {
    group[i]<-"High risk"
  }
}

Group<-cbind(names(Methylation),group)
Group<-as.data.frame(Group)
names(Group)<-c("sample","group")
row.names(Group)<-names(Methylation)
library(ComplexHeatmap)
#heatdata<-as.data.frame(t(normalize(t(heatdata)))) #FC4E07","#00AFBB"
Methylation<-as.matrix(Methylation)
pdf("./result/heatmapriskall.pdf",width = 8,height=3.5)
Heatmap(Methylation, cluster_rows = FALSE, cluster_columns = F,name = " ",show_column_names = F,
        column_title = c (paste(levels(Group$group)[1], "group"), paste(levels(Group$group)[2], "group"),paste(levels(Group$group)[3], "group")),
        column_title_gp = gpar( fill = c("#FC4E07","#00AFBB","green" ),col = "black",fontsize = 12 ),
        column_split = Group$group
)
dev.off()

####################################################################################
#survival analysis in the test set
rm(list=ls())
setwd("D:/learning/KIRC/methana")
library(survivalROC)
library(ggplot2)
library(survminer)
surv.test<-read.csv("./result/surv.test.csv",header = T,row.names=1)
summary(surv.test$OS.time)
years<-seq(from=0,to=max(surv.test$OS.time),by=365)
years<-years[-1]
sumROC<-list()
for (i in 1:length(years)){
  sumROC[[i]] <- survivalROC(Stime = surv.test$OS.time,status = surv.test$OS,marker = surv.test$risk,
                             predict.time =years[i],method = "NNE",span = 0.25*nrow(surv.test)^(-0.20))
}
sumAUC<-list()
for (i in 1:length(sumROC)){
  sumAUC[[i]]<-sumROC[[i]]$AUC
}
which.max(unlist(sumAUC))
minAUC<-min(sumAUC[[1]],sumAUC[[3]],sumAUC[[5]],sumAUC[[7]],sumAUC[[10]])

summary(sumROC[[1]])
ROCdata<-c()
for(i in 1:length(sumROC)){
  predict.time<-sumROC[[i]]$predict.time
  TP<-sumROC[[1]]$TP
  FP<-sumROC[[i]]$FP
  auc<-sumROC[[i]]$AUC
  tmp<-c(predict.time,TP,FP,auc)
  ROCdata<-rbind(ROCdata,tmp)
}

survivalROC_helper <- function(t) {
  survivalROC(Stime = surv.test$OS.time,status = surv.test$OS,marker = surv.test$risk,
              predict.time =t,method = "NNE",span = 0.25*nrow(surv.test)^(-0.20))
}
library(tidyverse)
library(survivalROC)
survivalROC_data <- data_frame(t = c(years[1],years[3],years[5],years[7],years[10])) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)

survivalROC_data1<-mutate(survivalROC_data,auc =sprintf("%.3f",auc))
survivalROC_data1$years<-survivalROC_data1$t/365
survivalROC_data1<-unite(survivalROC_data1,year, years,auc,sep = " year AUC: " )
AUC =factor(survivalROC_data1$year)
sumAUC1<-list(sumAUC[[1]],sumAUC[[3]],sumAUC[[5]],sumAUC[[7]],sumAUC[[10]])
sumROC1<-list(sumROC[[1]],sumROC[[3]],sumROC[[5]],sumROC[[7]],sumROC[[10]])
ROC.1<-sumROC1[[which.max(sumAUC1)]]
#ROC.1<-sumROC[[which.max(sumAUC)]]
cutoff.imm <- ROC.1$cut.values[with(ROC.1, which.min((1-TP)^2+ FP^2))]
cutoff.imm
dot <- data.frame(TP = ROC.1$TP[with(ROC.1, which.min((1-TP)^2+ FP^2))],
                  FP = ROC.1$FP[with(ROC.1, which.min((1-TP)^2+ FP^2))])
dot <- rbind(c(1,0),dot)
ROC.test<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
  ylab("True positive rate") + 
  theme(legend.position = c(0.8,0.2))+
  geom_path(mapping = aes(x = FP,y = TP),data = dot)+
  annotate("text",x = dot$FP[2] + 0.05,y = dot$TP[2],label = paste0("Cutoff: ",round(cutoff.imm,3)))
pdf("./result/ROC.test.pdf",width=5.5,height=5.5)
ROC.test
dev.off()
save(ROC.test,AUC,file="ROC.test.Rdata")
##################################################################################################################
group<- ifelse(surv.test$risk<=cutoff.imm,"Low risk", "High risk")
group
surv.test$Group<-group
write.csv(surv.test, "./result/test.surv.risk.csv")
##survival in the test set
library(survminer)
library(survival)
##survival analysis for test set for OS
rm(list=ls()) 
surv.test<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
levels(surv.test$hemoglobin_result)
surv.test$hemoglobin_result<-factor(surv.test$hemoglobin_result, levels=c("Elevated","Low","Normal"), labels=c("Elevated","Low","Normal"))
levels(surv.test$serum_calcium_result)
surv.test$serum_calcium_result<-factor(surv.test$serum_calcium_result, levels=c("Elevated","Low","Normal"), labels=c("Elevated","Low","Normal"))
names(surv.test)
levels(surv.test$histological_type)
surv.test$Group
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(surv.test$OS.time, surv.test$OS) ~ Group, data = surv.test)
# Basic survival curves #FC4E07","#00AFBB
level<-levels(factor(surv.test$Group))
c<-ifelse(level[1]=="Low risk",paste("#00AFBB"),paste("#FC4E07"))
d<-ifelse(level[2]=="High risk",paste("#FC4E07"),paste("#00AFBB"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
p.test <- ggsurvplot(fit, data = surv.test, risk.table = TRUE,
                     risk.table.height = 0.20,
                     risk.table.y.text = FALSE,
                     risk.table.title ="",
                     main = "Survival curve",
                     palette =c(c, d),pval=T, 
                     risk.table.y.text.col = T,
                     legend = c(0.85, 0.90),
                     legend.title = "",
                     xlab="Overall survival(days)",
                     legend.labs=c(e,f))
pdf(file="./result/KM.test.OS.pdf")
p.test
dev.off()
save(p.test,file="p.test.OS.Rdata")
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#log-rank test for the survival analysis
survdiff(Surv(surv.test$OS.time, surv.test$OS) ~Group,data=surv.test,rho=0)
#univariate cox analysis
fitu<-coxph(Surv(surv.test$OS.time, surv.test$OS) ~risk,data=surv.test)
summary(fitu)
names(surv.test)
#mutivariable cox
surv.test$race.demographic<-as.character(surv.test$race.demographic)
surv.test$race.demographic[surv.test$race.demographic == "not reported"] <- NA
names(surv.test)
surv.test$Group<-as.vector(surv.test$Group)
fitm<-coxph(Surv(surv.test$OS.time, surv.test$OS)~Group+age_at_initial_pathologic_diagnosis+hemoglobin_result
            +serum_calcium_result+tumor_stage.diagnoses+gender.demographic,data=surv.test)
fitm<-coxph(Surv(surv.test$OS.time, surv.test$OS)~Group+age_at_initial_pathologic_diagnosis
            +tumor_stage.diagnoses,data=surv.test)

fitm<-coxph(Surv(surv.test$OS.time, surv.test$OS)~risk+age_at_initial_pathologic_diagnosis+hemoglobin_result
            +serum_calcium_result+tumor_stage.diagnoses+gender.demographic,data=surv.test)
summary(fitm)
#fitm<-coxph(Surv(surv.test$OS, surv.test$OS_status)~Group+AGE+B2M+CRP+LDH+ALB+HGB+BMPC+MRI,data=surv.test)
sumfitm<-summary(fitm)
a<-cbind(sumfitm$conf.int[,1:4],sumfitm$coefficients[,5])
colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
a<-as.data.frame(a)
write.csv(a,file="./result/mcox.test.os.csv")

#unicox
Surv<-Surv(surv.test$OS.time, surv.test$OS)
var<-c("risk","age_at_initial_pathologic_diagnosis","hemoglobin_result","serum_calcium_result","tumor_stage.diagnoses","gender.demographic")
var
unicox<-c()
for (i in 1:length(var)){
  fomu<-as.formula(paste("Surv~",var[i],sep=""))
  cox<-coxph(fomu,data=surv.test)
  cox<-summary(cox)
  conf.int<-cox$conf.int
  coef1<-cox$coef
  a<-cbind(conf.int,coef1[,5])
  #a<-c(as.vector(conf.int),as.vector(coef1)[5])
  colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
  unicox<-rbind(unicox,a)
}
unicox
write.csv(unicox,file="./result/unicox.test.os.csv")

############################################################################################################
##survival analysis for test set for RFS
rm(list=ls()) 
surv.test<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
names(surv.test)
levels(surv.test$histological_type)
surv.test$Group
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(surv.test$RFS.time, surv.test$RFS) ~ Group, data = surv.test)
# Basic survival curves
level<-levels(factor(surv.test$Group))
c<-ifelse(level[1]=="Low risk",paste("#00AFBB"),paste("#FC4E07"))
d<-ifelse(level[2]=="High risk",paste("#FC4E07"),paste("#00AFBB"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
p.test <- ggsurvplot(fit, data = surv.test, risk.table = TRUE,
                     risk.table.height = 0.20,
                     risk.table.y.text = FALSE,
                     risk.table.title ="",
                     main = "Survival curve",
                     palette =c(c, d),pval=T, 
                     risk.table.y.text.col = T,
                     legend = c(0.85, 0.90),
                     legend.title = "",
                     xlab="Recurrence free survival(day)",
                     legend.labs=c(e,f))
pdf(file="./result/KM.test.RFS.pdf")
p.test
dev.off()
save(p.test,file="p.test.RFS.Rdata")
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#log-rank test for the survival analysis
survdiff(Surv(surv.test$RFS.time, surv.test$RFS) ~Group,data=surv.test,rho=0)
#univariate cox analysis
fitu<-coxph(Surv(surv.test$RFS.time, surv.test$RFS) ~risk,data=surv.test)
summary(fitu)
names(surv.test)
#mutivariable cox
names(surv.test)
fitm<-coxph(Surv(surv.test$RFS.time, surv.test$RFS)~Group+age_at_initial_pathologic_diagnosis+hemoglobin_result
            +serum_calcium_result+tumor_stage.diagnoses,data=surv.test)
summary(fitm)
#fitm<-coxph(Surv(surv.test$RFS., surv.test$RFS._status)~Group+AGE+B2M+CRP+LDH+ALB+HGB+BMPC+MRI,data=surv.test)
sumfitm<-summary(fitm)
a<-cbind(sumfitm$conf.int[,1:4],sumfitm$coefficients[,5])
colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
a<-as.data.frame(a)
write.csv(a,file="./result/mcox.test.RFS.csv")
#unicox
Surv<-Surv(surv.test$RFS.time, surv.test$RFS)
var<-c("Group","age_at_initial_pathologic_diagnosis","hemoglobin_result","serum_calcium_result","tumor_stage.diagnoses")
var
unicox<-c()
for (i in 1:length(var)){
  fomu<-as.formula(paste("Surv~",var[i],sep=""))
  cox<-coxph(fomu,data=surv.test)
  cox<-summary(cox)
  conf.int<-cox$conf.int
  coef1<-cox$coef
  a<-cbind(conf.int,coef1[,5])
  #a<-c(as.vector(conf.int),as.vector(coef1)[5])
  colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
  unicox<-rbind(unicox,a)
}
unicox
write.csv(unicox,file="./result/unicox.test.RFS.csv")
####################################################################################################
#survival analysis for the training set OS
rm(list=ls())
library(survivalROC)
library(ggplot2)
library(survminer)
surv.train<-read.csv("./result/surv.train.csv",header = T,row.names=1)
summary(surv.train$OS.time)
years<-seq(from=0,to=max(surv.train$OS.time),by=365)
years<-years[-1]
sumROC<-list()
for (i in 1:length(years)){
  sumROC[[i]] <- survivalROC(Stime = surv.train$OS.time,status = surv.train$OS,marker = surv.train$risk,
                             predict.time =years[i],method = "NNE",span = 0.25*nrow(surv.train)^(-0.20))
}
sumAUC<-list()
for (i in 1:length(sumROC)){
  sumAUC[[i]]<-sumROC[[i]]$AUC
}
which.max(unlist(sumAUC))
summary(sumROC[[1]])
survivalROC_helper <- function(t) {
  survivalROC(Stime = surv.train$OS.time,status = surv.train$OS,marker = surv.train$risk,
              predict.time =t,method = "NNE",span = 0.25*nrow(surv.train)^(-0.20))
}
library(tidyverse)
library(survivalROC)
survivalROC_data <- data_frame(t = c(years[1],years[3],years[5],years[7],years[10])) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
survivalROC_data1<-mutate(survivalROC_data,auc =sprintf("%.3f",auc))
survivalROC_data1$years<-survivalROC_data1$t/365
survivalROC_data1<-unite(survivalROC_data1,year, years,auc,sep = " year AUC: " )
AUC =factor(survivalROC_data1$year)
save(AUC,file="AUC.train.Rdata")
sumAUC1<-list(sumAUC[[1]],sumAUC[[3]],sumAUC[[5]],sumAUC[[7]],sumAUC[[10]])
sumROC1<-list(sumROC[[1]],sumROC[[3]],sumROC[[5]],sumROC[[7]],sumROC[[10]])
ROC.1<-sumROC1[[which.max(sumAUC1)]]
cutoff.imm <- ROC.1$cut.values[with(ROC.1, which.min((1-TP)^2+ FP^2))]
cutoff.imm
dot <- data.frame(TP = ROC.1$TP[with(ROC.1, which.min((1-TP)^2+ FP^2))],
                  FP = ROC.1$FP[with(ROC.1, which.min((1-TP)^2+ FP^2))])
dot <- rbind(c(1,0),dot)
ROC.train<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
  ylab("True positive rate") + 
  theme(legend.position = c(0.8,0.2))+
  geom_path(mapping = aes(x = FP,y = TP),data = dot)+
  annotate("text",x = dot$FP[2] + 0.05,y = dot$TP[2],label = paste0("Cutoff: ",round(cutoff.imm,3)))
pdf("./result/ROC.train.pdf",width=5.5,height=5.5)
ROC.train
dev.off()
save(ROC.train,AUC,file="ROC.train.Rdata")
##################################################################################################################
group<- ifelse(surv.train$risk<=cutoff.imm,"Low risk", "High risk")
group
surv.train$Group<-group
write.csv(surv.train, "./result/train.surv.risk.csv")
##survival in the train set
library(survminer)
library(survival)
##survival analysis for train set for OS
rm(list=ls()) 
surv.train<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
levels(surv.train$hemoglobin_result)
surv.train$hemoglobin_result<-factor(surv.train$hemoglobin_result, levels=c("Elevated","Low","Normal"), labels=c("Elevated","Low","Normal"))
levels(surv.train$serum_calcium_result)
surv.train$serum_calcium_result<-factor(surv.train$serum_calcium_result, levels=c("Elevated","Low","Normal"), labels=c("Elevated","Low","Normal"))
names(surv.train)
names(surv.train)
surv.train$Group
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(surv.train$OS.time, surv.train$OS) ~ Group, data = surv.train)
# Basic survival curves
level<-levels(factor(surv.train$Group))
c<-ifelse(level[1]=="Low risk",paste("#00AFBB"),paste("#FC4E07"))
d<-ifelse(level[2]=="High risk",paste("#FC4E07"),paste("#00AFBB"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
p.train <- ggsurvplot(fit, data = surv.train, risk.table = TRUE,
                      risk.table.height = 0.20,
                      risk.table.y.text = FALSE,
                      risk.table.title ="",
                      main = "Survival curve",
                      palette =c(c, d),pval=T, 
                      risk.table.y.text.col = T,
                      legend = c(0.85, 0.90),
                      legend.title = "",
                      xlab="Overall survival(days)",
                      legend.labs=c(e,f))
pdf(file="./result/KM.train.OS.pdf")
p.train
dev.off()
save(p.train,file="p.train.OS.Rdata")
###############################################################################################
#log-rank train for the survival analysis
survdiff(Surv(surv.train$OS.time, surv.train$OS) ~Group,data=surv.train,rho=0)
#univariate cox analysis
fitu<-coxph(Surv(surv.train$OS.time, surv.train$OS) ~risk,data=surv.train)
summary(fitu)
names(surv.train)
#mutivariable cox
surv.train$race.demographic[surv.train$race.demographic == "not reported"] <- NA
names(surv.train)
fitm<-coxph(Surv(surv.train$OS.time, surv.train$OS)~risk+age_at_initial_pathologic_diagnosis+hemoglobin_result
            +serum_calcium_result+tumor_stage.diagnoses+gender.demographic,data=surv.train)
summary(fitm)
#fitm<-coxph(Surv(surv.train$OS, surv.train$OS_status)~Group+AGE+B2M+CRP+LDH+ALB+HGB+BMPC+MRI,data=surv.train)
sumfitm<-summary(fitm)
a<-cbind(sumfitm$conf.int[,1:4],sumfitm$coefficients[,5])
colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
a<-as.data.frame(a)
write.csv(a,file="./result/mcox.train.os.csv")
#unicox
Surv<-Surv(surv.train$OS.time, surv.train$OS)
var<-c("risk","age_at_initial_pathologic_diagnosis","hemoglobin_result",
       "serum_calcium_result","tumor_stage.diagnoses","gender.demographic")
var
unicox<-c()
for (i in 1:length(var)){
  fomu<-as.formula(paste("Surv~",var[i],sep=""))
  cox<-coxph(fomu,data=surv.train)
  cox<-summary(cox)
  conf.int<-cox$conf.int
  coef1<-cox$coef
  a<-cbind(conf.int,coef1[,5])
  #a<-c(as.vector(conf.int),as.vector(coef1)[5])
  colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
  unicox<-rbind(unicox,a)
}
unicox
write.csv(unicox,file="./result/unicox.train.os.csv")
save(list=ls(),file="survivaltrain.OS.Rdata")
############################################

#multiplot 
rm(list=ls())
# multiple plot on one plot
load("p.train.OS.Rdata")
load("ROC.train.Rdata")
library(grid)
pdf("./result/train.OS.surv.pdf",width=11,height=6)
grid.newpage()     
pushViewport(viewport(layout=grid.layout(4,36)))
print(p.train$plot, vp=viewport(layout.pos.row=1:3, layout.pos.col=19:35))
print(p.train$table, vp=viewport(layout.pos.row=4, layout.pos.col=19:35))
print(ROC.train, vp=viewport(layout.pos.row=1:4, layout.pos.col=1:18))
dev.off()
rm(list=ls())
# multiple plot on one plot
load("p.test.OS.Rdata")
load("ROC.test.Rdata")


library(grid)
pdf("./result/test.OS.surv.pdf",width=11,height=6)
grid.newpage()     
pushViewport(viewport(layout=grid.layout(4,36)))
print(p.test$plot, vp=viewport(layout.pos.row=1:3, layout.pos.col=19:35))
print(p.test$table, vp=viewport(layout.pos.row=4, layout.pos.col=19:35))
print(ROC.test, vp=viewport(layout.pos.row=1:4, layout.pos.col=1:18))
dev.off()

#####################################################################################################################
##survival analysis for train set for RFS
rm(list=ls()) 
surv.train<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
names(surv.train)
surv.train$Group
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(surv.train$RFS.time, surv.train$RFS) ~ Group, data = surv.train)
# Basic survival curves
level<-levels(factor(surv.train$Group))
c<-ifelse(level[1]=="Low risk",paste("#00AFBB"),paste("#FC4E07"))
d<-ifelse(level[2]=="High risk",paste("#FC4E07"),paste("#00AFBB"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
p.train <- ggsurvplot(fit, data = surv.train, risk.table = TRUE,
                      risk.table.height = 0.20,
                      risk.table.y.text = FALSE,
                      risk.table.title ="",
                      main = "Survival curve",
                      palette =c(c, d),pval=T, 
                      risk.table.y.text.col = T,
                      legend = c(0.85, 0.90),
                      legend.title = "",
                      xlab="Recurrence-free survival (days)",
                      legend.labs=c(e,f))
pdf(file="./result/KM.train.RFS.pdf")
p.train
dev.off()
save(p.train,file="p.train.RFS.Rdata")

###############################################################################################
#log-rank train for the survival analysis
survdiff(Surv(surv.train$RFS.time, surv.train$RFS) ~Group,data=surv.train,rho=0)
#univariate cox analysis
fitu<-coxph(Surv(surv.train$RFS.time, surv.train$RFS) ~risk,data=surv.train)
summary(fitu)
names(surv.train)
#mutivariable cox
fitm<-coxph(Surv(surv.train$RFS.time, surv.train$RFS)~risk+age_at_initial_pathologic_diagnosis+hemoglobin_result
            +serum_calcium_result+tumor_stage.diagnoses+gender.demographic,data=surv.train)
summary(fitm)
#fitm<-coxph(Surv(surv.train$OS, surv.train$OS_status)~Group+AGE+B2M+CRP+LDH+ALB+HGB+BMPC+MRI,data=surv.train)
sumfitm<-summary(fitm)
a<-cbind(sumfitm$conf.int[,1:4],sumfitm$coefficients[,5])
colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
a<-as.data.frame(a)
write.csv(a,file="./result/mcox.train.RFS.csv")

#unicox
Surv<-Surv(surv.train$RFS.time, surv.train$RFS)
var<-c("risk","age_at_initial_pathologic_diagnosis","hemoglobin_result",
       "serum_calcium_result","tumor_stage.diagnoses","gender.demographic")
var
unicox<-c()
for (i in 1:length(var)){
  fomu<-as.formula(paste("Surv~",var[i],sep=""))
  cox<-coxph(fomu,data=surv.train)
  cox<-summary(cox)
  conf.int<-cox$conf.int
  coef1<-cox$coef
  a<-cbind(conf.int,coef1[,5])
  #a<-c(as.vector(conf.int),as.vector(coef1)[5])
  colnames(a)<-c("HR","exp(-coef)","LCI","UCI","P value")
  unicox<-rbind(unicox,a)
}
unicox
write.csv(unicox,file="./result/unicox.train.RFS.csv")
save(list=ls(),file="survivaltrain.RFS.Rdata")
######################################################################################################################
######################################################################################################################
#correlation analysis with mRNA expression
rm(list=ls())
setwd("D:/learning/KIRC/methana")
expres.mRNA<-read.table("D:/learning/KIRC/TCGA-KIRC.htseq_fpkm.tsv/TCGA-KIRC.htseq_fpkm.tsv",header=T,row.names=1)
sigunicox<-read.csv("./result/sigunicox.csv",header=T,row.names=1)
coef<-read.csv("./result/coef.csv",header=T,row.names=1)
genecode<-read.table("D:/learning/methods/supportdata/gencode.v22.annotation.gene.txt",header=T)
names(genecode)
names(sigunicox)
expres.mRNA$id<-row.names(expres.mRNA)
expres.mRNA<-merge(genecode,expres.mRNA,by="id")
expres.mRNA<-subset(expres.mRNA,select=-c(id,chrom,chromStart,chromEnd,strand))
dim(expres.mRNA)

sigunicox.coef<-sigunicox[which(rownames(sigunicox) %in% rownames(coef)),]
write.csv(sigunicox.coef,"./result/sigunicox.coef.csv")
#replce the fifth vector of the genesymbol in sigunicox.coef
sigunicox.coef$genesymbol<-as.vector(sigunicox.coef$genesymbol)
setdiff(sigunicox.coef$genesymbol,expres.mRNA$gene)
#clinical$tumor_stage.diagnoses[clinical$tumor_stage.diagnoses == "DLX6AS"] <- "DLX6-AS1"
"DLX6-AS1" %in% expres.mRNA$gene
sigunicox.coef$genesymbol[sigunicox.coef$genesymbol=="DLX6AS"]<-"DLX6-AS1"
sigunicox.coef$genesymbol
expres<-expres.mRNA[which(expres.mRNA$gene %in% sigunicox.coef$genesymbol), ]
if(nrow(expres)<nrow(sigunicox.coef)){
  setdiff(sigunicox.coef$genesymbol,expres$gene)
}

"DLX6-AS1" %in% expres.mRNA$gene
#a<-setdiff(sigunicox.coef$genesymbol,expres$gene)
#sigunicox.coef<-sigunicox.coef[-which(sigunicox.coef$genesymbol %in% a),]

row.names(expres)<-expres$gene
genes<-row.names(expres)
clinical<-read.csv("./result/clinical.csv",header=T,row.names=1)
names(expres)<-substr(colnames(expres),1,15)
row.names(clinical)<-gsub("-",".",row.names(clinical))
index<-intersect(row.names(clinical),names(expres))
expres<-expres[,index]
row.names(expres)
# two methylation loci matched to the same LTBP1,
RUNX3<-expres[which(row.names(expres)=="RUNX3"),]
expres<-rbind(expres,RUNX3)
expres<-rbind(expres,RUNX3)
expres<-rbind(expres,RUNX3)


dim(expres)
expres<-expres[order(row.names(expres)),]

can<-read.csv("./result/can.csv",header=T,row.names=1)
can<-can[,index]
sigunicox.coef<-sigunicox.coef[order(sigunicox.coef$genesymbol),]

can<-can[row.names(sigunicox.coef),]

identical(names(can),names(expres))
dim(expres)
dim(can)
row.names(expres)
row.names(can)

expres<-as.data.frame(t(expres))
#names(expres)[5]<-"LTBP1"
can<-as.data.frame(t(can))
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header = T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
train.test.surv<-rbind(train.surv.risk,test.surv.risk)
index<-intersect(row.names(can),row.names(clinical))
can<-can[index,]
expres<-expres[index,]
clinical<-clinical[index,]
train.test.surv<-train.test.surv[index,]
save(can,expres,clinical,train.test.surv,file="formediationanalysis.Rdata")
head(row.names(clinical))
head(row.names(can))
head(row.names(expres))
head(row.names(train.test.surv))
summary(can)
normalize<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}

can<-apply(can,2,normalize)
summary(can)
expres<-apply(expres,2,normalize)
summary(expres)
can<-as.data.frame(can)
expres<-as.data.frame(expres)
########################################################################################
# correlation analysis between the methylation level and mRNA expression level
#make surv the expres and the can have the same dim and order 
library("ggpubr")
coplot<-list()
for(i in 1: ncol(expres)){
  data<-cbind(expres[,i],can[,i])
  data<-as.data.frame(data)
  names(data)<-c("mRNA","methylation")
  data$Group<-train.test.surv$Group
  coplot[[i]]<- ggscatter(data, x = "mRNA", y = "methylation", 
                          color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                          add = "reg.line", conf.int = TRUE,
                          legend = c(0.15, 0.15),
                          legend.title = "",
                          add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                          cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                          cor.coeff.args = list(label.x = 0.5,label.y=1, label.sep = "\n"),
                          xlab = paste(names(expres)[i], "expression", sep = " "), ylab = paste(names(can)[i],"methylation",sep=" "))
}
coplot[[2]]

save(coplot,file="./result/coplot.Rdata")

#drow the survival plot based on the mRNA levels in the TCGA-BLCA
dim(clinical)
dim(can)
index<-intersect(row.names(clinical),rownames(can))
clinical<-clinical[index,]
dim(expres)
# Beacase DLX6-AS1 has been recognized as "DLX6.AS1", thus, we replace it with the previously form "DLX6AS"
names(expres)[names(expres) == "DLX6-AS1"] <- "DLX6AS"
head(row.names(clinical))
head(row.names(expres))
data<-cbind(clinical,expres)
names(data)
library(survminer)
res.cut <- surv_cutpoint(data, time = "OS.time", event = "OS",
                         variables = names(expres))
res.cat <- surv_categorize(res.cut)
library("survival")

KMplot<-list()
for (i in 1:length(names(expres))){
  a<-paste(names(expres)[i],"low expression group",sep=" ")
  b<-paste(names(expres)[i],"high expression group",sep=" ")
  level<-levels(factor(res.cat[,i+2]))
  c<-ifelse(level[1]=="low",paste("blue"),paste("orange"))
  d<-ifelse(level[2]=="high",paste("orange"),paste("blue"))
  e<-ifelse(level[1]=="low",a,b)
  f<-ifelse(level[2]=="high",b,a)
  fit<- survfit(as.formula(paste("Surv(OS.time, OS)~",names(expres)[i],sep=""))
                , data = res.cat)
  KMplot[[i]]<-ggsurvplot(fit, data = res.cat,
                          pval=T, pval.size=4,pval.coord=c(0.15,0.05),
                          pval.method=T,pval.method.size=4,pval.method.coord=c(0.15,0.1),
                          legend = c(0.74, 0.9),
                          legend.title = "",
                          palette =c(c, d),
                          xlab="Overall survival (days)",
                          legend.labs=c(e,f)
  )
}
coplot[[5]]
KMplot[[6]]
save(KMplot,file="./result/KMplot.Rdata")

names(expres)
variable<-paste(names(expres),collapse="+")
Surv<-Surv(clinical$OS.time,clinical$OS)
fomu<-as.formula(paste("Surv~",variable,seq=""))

mcox.mRNA<-coxph(fomu,data=expres)
risk<-predict(mcox.mRNA,type="risk")
clinical$risk.mRNA<-risk

library(survivalROC)
library(ggplot2)
library(survminer)

summary(clinical$OS.time)
years<-seq(from=0,to=max(clinical$OS.time),by=365)
years<-years[-1]
sumROC<-list()
for (i in 1:length(years)){
  sumROC[[i]] <- survivalROC(Stime = clinical$OS.time,status = clinical$OS,marker = clinical$risk.mRNA,
                             predict.time =years[i],method = "NNE",span = 0.25*nrow(clinical)^(-0.20))
}
sumAUC<-list()
for (i in 1:length(sumROC)){
  sumAUC[[i]]<-sumROC[[i]]$AUC
}
which.max(unlist(sumAUC))
survivalROC_helper <- function(t) {
  survivalROC(Stime = clinical$OS.time,status = clinical$OS,marker = clinical$risk,
              predict.time =t,method = "NNE",span = 0.25*nrow(clinical)^(-0.20))
}
library(tidyverse)
library(survivalROC)
survivalROC_data <- data_frame(t = c(years[1],years[3],years[5],years[7],years[10])) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)

survivalROC_data1<-mutate(survivalROC_data,auc =sprintf("%.3f",auc))
survivalROC_data1$years<-survivalROC_data1$t/365
survivalROC_data1<-unite(survivalROC_data1,year, years,auc,sep = " year AUC: " )
AUC =factor(survivalROC_data1$year)
#save(AUC,file="AUC.mRNA.Rdata")
sumAUC1<-list(sumAUC[[1]],sumAUC[[3]],sumAUC[[5]],sumAUC[[7]],sumAUC[[10]])
sumROC1<-list(sumROC[[1]],sumROC[[3]],sumROC[[5]],sumROC[[7]],sumROC[[10]])
ROC.1<-sumROC1[[which.max(sumAUC1)]]
cutoff.imm <- ROC.1$cut.values[with(ROC.1, which.min((1-TP)^2+ FP^2))]
cutoff.imm
dot <- data.frame(TP = ROC.1$TP[with(ROC.1, which.min((1-TP)^2+ FP^2))],
                  FP = ROC.1$FP[with(ROC.1, which.min((1-TP)^2+ FP^2))])
dot <- rbind(c(1,0),dot)
ROC.train<-ggplot(survivalROC_data1,mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+xlab("False positive rate")+
  ylab("True positive rate") + 
  theme(legend.position = c(0.8,0.2))

#geom_path(mapping = aes(x = FP,y = TP),data = dot)
#annotate("text",x = dot$FP[2] + 0.05,y = dot$TP[2],label = paste0("Cutoff: ",round(cutoff.imm,3))
Diff.methylation<-read.csv("./result/Diff.methylation.csv",header=T,row.names=1)
sig.methylation<-Diff.methylation[row.names(sigunicox.coef),]
write.csv(sig.methylation,"./result/sig.methylation.csv")

pdf("./result/ROC.mRNA.pdf",width=5.5,height=5.5)
ROC.train
dev.off()
save(ROC.train,AUC,file="./result/ROC.mRNA.Rdata")
#multiplot
rm(list=ls())
# multiple plot on one plot
load("./result/KMplot.Rdata")
load("./result/coplot.Rdata")
load("./result/ROC.mRNA.Rdata")
library(grid)
pdf("./result/mRNA.methy.pdf",width=20,height=30)
grid.newpage()     
pushViewport(viewport(layout=grid.layout(30,20)))
print(coplot[[1]], vp=viewport(layout.pos.row=1:5, layout.pos.col=1:5))
print(coplot[[2]], vp=viewport(layout.pos.row=6:10, layout.pos.col=1:5))
print(coplot[[3]], vp=viewport(layout.pos.row=11:15, layout.pos.col=1:5))
print(coplot[[4]], vp=viewport(layout.pos.row=16:20, layout.pos.col=1:5))
print(coplot[[5]], vp=viewport(layout.pos.row=21:25, layout.pos.col=1:5))
print(coplot[[6]], vp=viewport(layout.pos.row=26:30, layout.pos.col=1:5))

print(KMplot[[1]]$plot, vp=viewport(layout.pos.row=1:5, layout.pos.col=6:10))
print(KMplot[[2]]$plot, vp=viewport(layout.pos.row=6:10, layout.pos.col=6:10))
print(KMplot[[3]]$plot, vp=viewport(layout.pos.row=11:15, layout.pos.col=6:10))
print(KMplot[[4]]$plot, vp=viewport(layout.pos.row=16:20, layout.pos.col=6:10))
print(KMplot[[5]]$plot, vp=viewport(layout.pos.row=21:25, layout.pos.col=6:10))
print(KMplot[[6]]$plot, vp=viewport(layout.pos.row=26:30, layout.pos.col=6:10))

print(coplot[[7]], vp=viewport(layout.pos.row=1:5, layout.pos.col=11:15))
print(coplot[[8]], vp=viewport(layout.pos.row=6:10, layout.pos.col=11:15))
print(coplot[[9]], vp=viewport(layout.pos.row=11:15, layout.pos.col=11:15))
print(coplot[[10]], vp=viewport(layout.pos.row=16:20, layout.pos.col=11:15))
print(coplot[[11]], vp=viewport(layout.pos.row=21:25, layout.pos.col=11:15))

print(KMplot[[7]]$plot, vp=viewport(layout.pos.row=1:5, layout.pos.col=16:20))
print(KMplot[[8]]$plot, vp=viewport(layout.pos.row=6:10, layout.pos.col=16:20))
print(KMplot[[9]]$plot, vp=viewport(layout.pos.row=11:15, layout.pos.col=16:20))
print(KMplot[[10]]$plot, vp=viewport(layout.pos.row=16:20, layout.pos.col=16:20))
print(KMplot[[11]]$plot, vp=viewport(layout.pos.row=21:25, layout.pos.col=16:20))
dev.off()

# valid in indedpendent dataset GSE113501
rm(list=ls())
setwd("D:/learning/KIRC/methana")
clinical<-read.csv("D:/learning/KIRC/GSE113501_series_matrix.txt/GSE113501_clinical.csv", header=T,row.names=1)
methylation<-read.csv("D:/learning/KIRC/GSE113501_series_matrix.txt/GSE113501_series_matrix.csv", header=T,row.names=1)
coef<-read.csv("./result/coef.csv",header=T,row.names=1)
methylation<-methylation[row.names(coef),]
dim(clinical)
dim(methylation)
clinical<-clinical[complete.cases(clinical$Age),]
clinical<-clinical[which(clinical$Metastatsis %in% c("M0-P","M0-PF","TNMIV")),]
clinical$Group<-ifelse(as.character(clinical$Metastatsis)=="TNMIV","Metastasis","Non-metastasis")
index<-intersect(row.names(clinical),names(methylation))
clinical<-clinical[index,]
methylation<-methylation[,index]
methylation<-na.omit(methylation)
dim(methylation)
index<-intersect(row.names(methylation),row.names(coef))
methylation<-methylation[index,]
coef<-coef[index,]
clinical$risk<-colSums(methylation*coef)
library(pROC)
library(ggplot2)
library(Rmisc)
SuSE<-summarySE(clinical, measurevar="risk", groupvars=c("Group"))
modelroc <- roc(clinical$Group,clinical$risk)
roc<-ggroc(modelroc)+geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
  theme_bw()+ annotate("text",x = 0.2,y = 0.2,label = paste0("AUC: ",round(modelroc$auc,3)))
pdf("./result/ROC.Meta.GSE113501.pdf")
roc
dev.off()
save(roc,file="ROC.Meta.GSE113501.Rdata")
p<-ggplot(data=SuSE,aes(x=Group,y=risk,fill=Group))+geom_bar(position=position_dodge(), stat="identity")
#进行坐标轴，图例等的微调
p1<-p+ylab("Methylation risk")+xlab("")

p2<-p1+theme_bw()+theme(legend.position='bottom')+labs(fill='')  #设置背景为空，存放图例位置及设置图例的名???
p3<-p2+theme(legend.text = element_text(colour="black", size = 12, face = "bold"))
p4<-p2+theme(axis.text.x=element_text(size=12,angle=0,color="Black"),
             axis.text.y=element_text(size=12,angle=0,color="Black")
)
p5<-p4+theme(axis.title.y = element_text(size = 12, color = "Black"))+
  geom_errorbar(aes(ymin=risk-se,ymax=risk+se),width=0.1,size=0.5,position = position_dodge(width = 0.9))+
  guides(fill=FALSE)
save(p5,file="GSE113501.barplot.Rdata")
pdf("./result/GSE113501.barplot.pdf")
p5
dev.off()
#boxplot
library(data.table)
library(ggplot2)
library(ggsignif)
library("ggpubr")
p <- ggboxplot(clinical, x="Group", y="risk",fill="Group",
               xlab = FALSE,ylab = "Methylation risk",font.y=12,font.xtickslab = 12,
               font.ytickslab = 12)+
  theme(legend.position='none')+
  stat_compare_means(method = "wilcox.test",label= "p.format",label.x.npc=0.45,hide.ns=T)
pdf("./result/boxplot.GSE113501.pdf")
p
dev.off()
save(p,file="boxplot.GSE113501.Rdata")
#validation immune cells in GSE113501
rm(list=ls())
setwd("D:/learning/KIRC/methana")
clinical<-read.csv("D:/learning/KIRC/GSE113501_series_matrix.txt/GSE113501_clinical.csv", header=T,row.names=1)
methylation<-read.csv("D:/learning/KIRC/GSE113501_series_matrix.txt/GSE113501_series_matrix.csv", header=T,row.names=1)

centEpiFibIC.m<-load("D:/learning/methods/rfunction/EpiDISH_2.0.2/EpiDISH/data/centEpiFibIC.m.rda")
centBloodSub.m<-load("D:/learning/methods/rfunction/EpiDISH_2.0.2/EpiDISH/data/centBloodSub.m.rda")
source("D:/learning/methods/rfunction/EpiDISH_2.0.2/EpiDISH/R/hepidish.R")
source("D:/learning/methods/rfunction/EpiDISH_2.0.2/EpiDISH/R/epidish.R")
frac.m <- hepidish(beta.m = methylation, ref1.m = centEpiFibIC.m,ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
#validation in GSE61441
rm(list=ls())
setwd("D:/learning/KIRC/methana")
clinical<-read.csv("D:/learning/KIRC/GSE61441_series_matrix.txt/GSE61441_clinical.csv",header=T,row.names=1)
methylation<-read.csv("D:/learning/KIRC/GSE61441_series_matrix.txt/GSE61441_series_matrix.csv",header=T,row.names=1)
dim(clinical)
dim(methylation)
coef<-read.csv("./result/coef.csv",header=T,row.names=1)
methylation<-methylation[row.names(coef),]
dim(methylation)
index<-intersect(row.names(clinical),names(methylation))
clinical<-clinical[index,]
methylation<-methylation[,index]
methylation<-na.omit(methylation)
dim(methylation)
index<-intersect(row.names(methylation),row.names(coef))
methylation<-methylation[index,]
coef<-coef[index,]
clinical$risk<-colSums(methylation*coef)

modelroc <- roc(clinical$Group,clinical$risk)
roc<-ggroc(modelroc)+geom_abline(intercept = 1, slope = 1, linetype = "dashed")+
  theme_bw()+ annotate("text",x = 0.2,y = 0.2,label = paste0("AUC: ",round(modelroc$auc,3)))
SuSE<-summarySE(clinical, measurevar="risk", groupvars=c("Group"))
pdf("./result/ROC.nomal.pdf")
roc
dev.off()
save(roc,file="ROC.nomal.Rdata")
p<-ggplot(data=SuSE,aes(x=Group,y=risk,fill=Group))+geom_bar(position=position_dodge(),stat="identity") #绘制柱状图，核心代码就一???
#进行坐标轴，图例等的微调
p1<-p+ylab("Methylation risk")+xlab("")

p2<-p1+theme_bw()+theme(legend.position='bottom')+labs(fill='')  #设置背景为空，存放图例位置及设置图例的名???
p3<-p2+theme(legend.text = element_text(colour="black", size = 12, face = "bold"))
p4<-p2+theme(axis.text.x=element_text(size=12,angle=0,color="Black"),
             axis.text.y=element_text(size=12,angle=0,color="Black")
)
p5<-p4+theme(axis.title.y = element_text(size = 12, color = "Black"))+
  geom_errorbar(aes(ymin=risk-se,ymax=risk+se),width=0.1,size=0.5,position = position_dodge(width = 0.9))+
  guides(fill=FALSE)
save(p5,file="./result/GSE61441.bar.Rdata")

#boxplot
library(data.table)
library(ggplot2)
library(ggsignif)
library("ggpubr")
p <- ggboxplot(clinical, x="Group", y="risk",fill="Group",
               xlab = FALSE,ylab = "Methylation risk",font.y=12,font.xtickslab = 12,
               font.ytickslab = 12)+
  theme(legend.position='none')+
  stat_compare_means(method = "wilcox.test",label= "p.format",label.x.npc=0.3,hide.ns=T)
save(p,file="./result/GSE61441.box.Rdata")
#multiplot
rm(list=ls())
load("./result/GSE61441.bar.Rdata")
GSE61441bar<-p5
load("ROC.nomal.Rdata")
GSE61441ROC<-roc
load("GSE113501.barplot.Rdata")
GSE113501bar<-p5
load("ROC.Meta.GSE113501")
GSE113501roc<-roc
pdf("./result/validation.pdf",height=6,width=7)
ggarrange(
  GSE61441bar,
  GSE61441ROC,
  GSE113501bar,
  GSE113501roc,
  ncol = 2,
  nrow = 2,
  labels = "AUTO",
  align = "hv",
  font.label=list(face = "plain",family="Times")
)
dev.off()
load("./result/seed5497.Rdata")
mode(clinical$Age)
####
#GSEAanalysis
setwd("D:/learning/KIRC/methana")
expres<-read.table("D:/learning/KIRC/TCGA-KIRC.htseq_fpkm.tsv/TCGA-KIRC.htseq_fpkm.tsv",header=T,row.names=1)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
names(expres)<-substr(names(expres),1,15)
dim(survdata)
names(survdata)
row.names(survdata)
names(expres)
genecode<-read.table("D:/learning/methods/supportdata/gencode.v22.annotation.gene.txt",header=T)
names(genecode)
names(sigunicox)
expres$id<-row.names(expres)
expres<-merge(genecode,expres,by="id")
expres<-subset(expres,select=-c(id,chrom,chromStart,chromEnd,strand))

index<-intersect(names(expres),row.names(train.surv.risk))
train.surv.risk<-train.surv.risk[index,]
expres1<-expres[,index]
expres<-cbind(expres$gene,expres1)

write.csv(expres,"D:/learning/KIRC/methana/GSEA/expres.csv")
write.csv(train.surv.risk,"D:/learning/KIRC/methana/GSEA/train.surv.risk.csv")
a<-rowSums(expres[,-1])
summary(a)
a<-as.data.frame(a)
expres$a<-a
names(expres)
expres<-expres[which(expres$a!=0),]

#CIBERSORT
setwd("D:/learning/KIRC/methana")
rm(list=ls())
expres<-read.table("D:/learning/KIRC/TCGA_KIRC_tpm.tsv/TCGA_KIRC_tpm.tsv",header=T,row.names=1)
dim(expres)
row.names(expres)
na<-strsplit(row.names(expres), "|", fixed=TRUE)   
expres<-cbind(na[[6]],expres)
source("D:/learning/methods/rfunction/metaDE.match/Match.gene.txt")
expres<-Match.gene(expres,"IQR")
dim(expres)

length(names(expres))
expres<-as.data.frame(expres)
names(expres) <- gsub("X","",names(expres))
metaTCGA<-read.csv("D:/learning/KIRC/TCGA_Metadata.csv/TCGA_Metadata.csv",header=T)
metaTCGA$a_CGHubAnalysisID <- gsub("-",".",metaTCGA$a_CGHubAnalysisID)
index<-intersect(names(expres),metaTCGA$a_CGHubAnalysisID)
KiRCmeta<-metaTCGA[which(metaTCGA$a_CGHubAnalysisID %in% index),]
identical(KiRCmeta$a_CGHubAnalysisID,names(expres))
head(KiRCmeta$a_CGHubAnalysisID)
head(names(expres))
names(expres)<-substr(KiRCmeta$a_AliquotBarcode,1,15)
source("D:/learning/methods/immuneinfilltration/cibersort.R")
sig_matrix<-read.table("D:/learning/methods/immuneinfilltration/ref.new.txt",header=T,sep="\t",row.names=1,check.names=F)
KIRCTME<-CIBERSORT(expres,perm=1000,QN=T)
#save(KIRCTME,file="./result/KIRCTME.Rdata")
rm(list=ls())
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
load("./result/KIRCTME.Rdata")
KIRCTME<-as.data.frame(KIRCTME)
KIRCTME<-KIRCTME[KIRCTME$`P-value`<0.05,]

index<-intersect(row.names(KIRCTME),row.names(survdata))
KIRCTME<-KIRCTME[index,]
survdata<-survdata[index,]
KIRCTME<-KIRCTME[,1:22]
normalize<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}

KIRCTME<-apply(KIRCTME,2,normalize)
KIRCTME<-as.data.frame(KIRCTME)

library("ggpubr")
coplot<-list()
for(i in 1: ncol(KIRCTME)){
  data<-cbind(survdata$risk,KIRCTME[,i])
  data<-as.data.frame(data)
  names(data)<-c("Methylation","Immune")
  data$Group<-survdata$Group
  coplot[[i]]<- ggscatter(data, y = "Methylation", x = "Immune", 
                          color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                          add = "reg.line", conf.int = TRUE,
                          legend = c(0.15, 0.15),
                          legend.title = "",
                          add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                          cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                          cor.coeff.args = list(label.x = 0.5,label.y=1.02, label.sep = "\n"),
                          xlab =paste(names(KIRCTME)[i]) , ylab = "Methylation risk")}

coplot[[22]]



CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] ) }

data<-cbind(survdata$risk,KIRCTME)
library(Hmisc)
res <- rcorr(as.matrix(data))
res<-CorMatrix (res$r, res$P)
res<-res[which(res$row=="survdata$risk"),]
res<-na.omit(res)
res$Adjusted.P<-p.adjust(res$p,method="bonferroni")
res<-res[res$Adjusted.P<0.05,]

library(ggplot2)
correplot<-ggplot(data=res)+geom_bar(aes(x=column,y=cor, fill=p), stat='identity')+
  coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+ 
  ylab("Spearman rho")+
  theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="black", size=48))+
  scale_y_continuous(expand=c(0, 0))+ 
  scale_x_discrete(expand=c(0,0))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
save(correplot,file="./result/correplot.bar.immune.Rdata")
correplot<-ggplot(data=res,aes(x=cor,y=column))+geom_point(aes(col=cor,size=p))+
  scale_color_gradient(low="blue",high="red")+
  labs(color="R square", size="P value")+theme_bw()+
  ylab("")+xlab("R square")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
        axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))
save(correplot,file="./result/correplot.immnue.point.Rdata")

pdf("./result/correplot.pdf",width=6,height=5)
correplot
dev.off()

library(stats)
library(CCA)
library(mvstats)
index<-intersect(res$column,names(KIRCTME))
TME<-KIRCTME[,index]
x<-TME
y<-survdata$risk
ca<-cancor(scale(x),scale(y))
cancor.test(x,y,plot=T)

U.TME <- as.matrix(scale(x)) %*% ca$xcoef
V.TME <- as.matrix(scale(y)) %*% ca$ycoef
TMEdata<-data.frame(U.TME[,1],V.TME,survdata$Group)
names(TMEdata)<-c("signature","risk", "Group")

library(ggpubr)
#if the P is not equal to can.test, try pearson test
coplot.TME.CCA<- ggscatter(data=TMEdata, y = "risk",x = "signature", 
                           color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                           add = "reg.line", conf.int = TRUE,
                           legend = c(0.15, 0.98),
                           legend.title = "",
                           add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                           cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                           cor.coeff.args = list(label.x =0.16,label.y=0.28, label.sep = "\n"),
                           xlab ="Immune cell infiltration signature" , ylab = "Methylation risk")
save(coplot.TME.CCA,file="coplot.TME.CCA.Rdata")

#survival analysis of the total effect of cell TME gene on the survival KIRC patients
U.TME<-as.data.frame(U.TME)
if(identical(row.names(U.TME),row.names(survdata))) {
  survdata$U.TME.group<-ifelse(as.numeric(U.TME[,1])<=median(as.numeric(U.TME[,1])),"Low risk","High risk")
}
#
library(survminer)
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(survdata$OS.time, survdata$OS) ~ U.TME.group, data = survdata)
survdiff(Surv(survdata$OS.time, survdata$OS) ~ U.TME.group, data = survdata)
summary(coxph(Surv(survdata$OS.time, survdata$OS) ~ U.TME.group, data = survdata))
# Basic survival curves
level<-levels(factor(survdata$U.TME.group))
c<-ifelse(level[1]=="Low risk",paste("blue"),paste("red"))
d<-ifelse(level[2]=="High risk",paste("red"),paste("blue"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
KMTME <- ggsurvplot(fit, data = survdata, risk.table = F,
                    risk.table.height = 0.20,
                    risk.table.y.text = FALSE,
                    risk.table.title ="",
                    main = "Survival curve",
                    palette =c(c, d),pval=T, 
                    risk.table.y.text.col = T,
                    legend = c(0.85, 0.90),
                    legend.title = "",
                    xlab="Overall survival (days)",
                    legend.labs=c(e,f))
save(KMTME,file="KMcellTME.Rdata")

#estimate analysis
setwd("D:/learning/KIRC/methana")
rm(list=ls())
library(estimate)
help(package="estimate")

#MACH analyis
setwd("D:/learning/KIRC/methana")
rm(list=ls())
library(maftools)
maf = read.maf(maf="D:/learning/KIRC/CGA.KIRC.varscan.ee42d944-ccbf-4406-9e7c-ffe1a0a4a1d7.DR-10.0.somatic.maf/TCGA.KIRC.varscan.somatic.maf")
getFields(maf)
getClinicalData(maf) 
getSampleSummary(maf)
plotmafSummary(maf = maf, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)#
oncoplot(maf = maf, top = 20, fontSize = 12)

maf@data$VAF <- maf@data$t_alt_count/maf@data$t_depth
head(maf@data$Tumor_Sample_Barcode)

barcode<-levels(maf@data$Tumor_Sample_Barcode)
rem.sam<- which(barcode=="TCGA-B8-A8YJ-01A-13D-A38X-10")#remove bad sample with no t_vaf
barcode<-barcode[-rem.sam]
rem.sam<- which(barcode=="TCGA-DV-5576-01A-01D-1534-10")
barcode<-barcode[-rem.sam]
MATH.score<-c()
for(i in 1:length(barcode)){
  res <- inferHeterogeneity(maf = maf, tsb = barcode[i], vafCol = 'VAF')
  MATH.score[i]<-unique(res$clusterData$MATH)
}
MATH<-data.frame(barcode,MATH.score)
MATH$barcode<-gsub("-",".",MATH$barcode)
row.names(MATH)<-substr(MATH$barcode,1,15)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
index<-intersect(row.names(MATH),row.names(survdata))
length(index)
survdata<-survdata[index,]
MATH<-MATH[index,]
MATH$risk<-survdata$risk
MATH$Group<-survdata$Group
#correlation analysis
library("ggpubr")
coplot.MATH<- ggscatter(data=MATH, y = "risk", x = "MATH.score", 
                        color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                        add = "reg.line", conf.int = TRUE,
                        legend = c(0.15, 0.98),
                        legend.title = "",
                        add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                        cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                        cor.coeff.args = list(label.x = 40,label.y=1.016, label.sep = "\n"),
                        xlab ="MATH score" , ylab = "Methylation risk")

save(coplot.MATH,file="coplot.MATH.Rdata")

#MSI correlation analysis
setwd("D:/learning/KIRC/methana")
rm(list=ls())
library(BiocOncoTK)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
dingMSI
row.names(dingMSI)<-gsub("-",".",row.names(dingMSI))
row.names(survdata)<-substr(row.names(survdata),1,12)
index<-intersect(row.names(dingMSI),row.names(survdata))
length(index)
survdata<-survdata[index,]
MSI<-dingMSI[index,]
MSI$risk<-survdata$risk
MSI$Group<-survdata$Group
#MSI<-MSI[MSI$MSIsensor.score!=0,]
coplot.MSI<- ggscatter(data=MSI, y = "risk",x = "MSIsensor.score", 
                       color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                       add = "reg.line", conf.int = TRUE,
                       legend = c(0.15, 0.15),
                       legend.title = "",
                       add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                       cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                       cor.coeff.args = list(label.x =15,label.y=1.015, label.sep = "\n"),
                       xlab ="MSIsensor.score" , ylab = "Methylation risk")

MSI$MSIgroup<-ifelse(MSI$MSIsensor.score>=4,"MSI-high","MSI-low")
table(MSI$MSIgroup)
#SCNA analysis
#cellcycle_genes and cytokine_genes
rm(list=ls())
setwd("D:/learning/KIRC/methana")
library(KEGG.db)
library(org.Hs.eg.db)
#entrez2symbol<-as.data.frame(org.Hs.egALIAS2EG)
#names(entrez2symbol)<-c("entrezID","symbol")
#write.table(entrez2symbol,file="D:/learning/methods/supportdata/entrz2symbol.txt",sep = "\t",row.names=F)
#entrez2ensembl<-as.data.frame(org.Hs.egENSEMBL)
#names(entrez2ensembl)<-c("entrezID","ensembl")
#write.table(entrez2ensembl,file="D:/learning/methods/supportdata/entrz2ensembl.txt",sep = "\t",row.names=F)
#entrez_symbol_ensembl<-merge(entrez2ensembl,entrez2symbol,by="entrezID")
#write.table(entrez_symbol_ensembl,file="D:/learning/methods/supportdata/entrez_symbol_ensembl.txt",sep = "\t",row.names=F)

ls("package:KEGG.db")
cellcycle_genes=KEGGPATHID2EXTID[['hsa04110']]
cytokine<-KEGGPATHID2EXTID[['hsa04060']]
expres.mRNA<-read.table("D:/learning/KIRC/TCGA-KIRC.htseq_fpkm.tsv/TCGA-KIRC.htseq_fpkm.tsv",header=T,row.names=1)
expres.mRNA<-expres.mRNA[,substr(colnames(expres.mRNA), 14, 16) %in% c("01A","11A")]
genecode<-read.table("D:/learning/methods/supportdata/gencode.v22.annotation.gene.txt",header=T)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
expres.mRNA$id<-row.names(expres.mRNA)
expres.mRNA<-merge(genecode,expres.mRNA,by="id")
expres.mRNA<-subset(expres.mRNA,select=-c(id,chrom,chromStart,chromEnd,strand))
entrez2symbol<-read.table("D:/learning/methods/supportdata/entrz2symbol.txt",header=T,sep="\t")
names(expres.mRNA)<-substr(names(expres.mRNA),1,15)
names(expres.mRNA)[1]<-"symbol"
expres.mRNA<-merge(entrez2symbol,expres.mRNA,by="symbol")

######################################################################################################
expres.mRNA.cellcycle<-expres.mRNA[expres.mRNA$entrezID %in% cellcycle_genes,]
expres.mRNA.cellcycle<-subset(expres.mRNA.cellcycle,select=-entrezID)
source("D:/learning/methods/rfunction/metaDE.match/Match.gene.txt")
expres.mRNA.cellcycle<-Match.gene(expres.mRNA.cellcycle,"IQR")
expres.mRNA.cellcycle<-as.data.frame(expres.mRNA.cellcycle)

index<-intersect(names(expres.mRNA.cellcycle),rownames(survdata))
expres.mRNA.cellcycle<-expres.mRNA.cellcycle[,index]
survdata<-survdata[index,]

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] ) }

data<-cbind(survdata$risk,as.data.frame(t(expres.mRNA.cellcycle)))
library(Hmisc)
res <- rcorr(as.matrix(data),type="spearman")
res<-CorMatrix (res$r, res$P)
res<-res[which(res$row=="survdata$risk"),]
res<-na.omit(res)
#res<-res[res$p<0.05,]
res$Adjusted.P<-p.adjust(res$p,method = "bonferroni")
res<-res[res$Adjusted.P<0.05,]
write.csv(res,"./result/res.cellcycle.csv")
library(ggplot2)
correplot.cycle.bar<-ggplot(data=res)+geom_bar(aes(x=column,y=cor, fill=p), stat='identity')+
  coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+ 
  ylab("Correlation coefficient")+
  theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="black", size=48))+
  scale_y_continuous(expand=c(0, 0))+ 
  scale_x_discrete(expand=c(0,0))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
save(correplot.cycle.bar,file="correplot.cycle.bar.Rdata")

correplot.cycle.point<-ggplot(data=res,aes(x=cor,y=column))+geom_point(aes(col=cor,size=p))+
  scale_color_gradient(low="blue",high="red")+
  labs(color="R square", size="P value")+theme_bw()+
  ylab("")+xlab("R square")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
        axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))

library(stats)
library(CCA)
library(mvstats)
index<-intersect(res$column,row.names(expres.mRNA.cellcycle))
expres.mRNA.cellcycle<-expres.mRNA.cellcycle[index,]
ca<-cancor(scale(t(expres.mRNA.cellcycle)),scale(survdata$risk))
cancor.test(t(expres.mRNA.cellcycle),survdata$risk,plot=T)

corcoef.test <- function(r, n, p, q, alpha = 0.1) {
  #r为相关系??? n为样本个??? 且n>p+q
  m <- length(r)
  Q <- rep(0, m)
  lambda <- 1
  for (k in m:1) {
    lambda <- lambda * (1 - r[k] ^ 2)
    #检验统计量
    Q[k] <- -log(lambda)   #检验统计量取对???
  }
  s <- 0
  i <- m
  for (k in 1:m) {
    Q[k] <- (n - k + 1 - 1 / 2 * (p + q + 3) + s) * Q[k]  #统计???
    chi <- 1 - pchisq(Q[k], (p - k + 1) * (q - k + 1))
    if (chi > alpha) {
      i <- k - 1
      break
    }
    s <- s + 1 / r[k] ^ 2
  }
  i  #显示输出结果 选用第几对典型变???
}

n<-corcoef.test(r=ca$cor,n=ncol(expres.mRNA.cellcycle),p=nrow(expres.mRNA.cellcycle),q=1)
U.cycle <- as.matrix(scale(t(expres.mRNA.cellcycle))) %*% ca$xcoef
V.cycle <- as.matrix(scale(survdata$risk)) %*% ca$ycoef
cycledata<-data.frame(U.cycle[,n],V.cycle,survdata$Group)
names(cycledata)<-c("signature","risk", "Group")

library(ggpubr)
coplot.cycle.CCA<- ggscatter(data=cycledata, y = "risk",x = "signature", 
                             color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                             add = "reg.line", conf.int = TRUE,
                             legend = c(0.15, 0.98),
                             legend.title = "",
                             add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                             cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                             cor.coeff.args = list(label.x =0,label.y=0.15, label.sep = "\n"),
                             xlab ="Cell cycle gene signature" , ylab = "Methylation risk")
save(coplot.cycle.CCA,file="coplot.cycle.CCA.Rdata")
#survival analysis of the total effect of cell cycle gene on the survival KIRC patients
U.cycle<-as.data.frame(U.cycle)
if(identical(row.names(U.cycle),row.names(survdata))) {
  survdata$U.cycle.group<-ifelse(as.numeric(U.cycle[,n])<=median(as.numeric(U.cycle[,n])),"Low risk","High risk")
}
#
library(survminer)
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(survdata$OS.time, survdata$OS) ~ U.cycle.group, data = survdata)
survdiff(Surv(survdata$OS.time, survdata$OS) ~ U.cycle.group, data = survdata)
summary(coxph(Surv(survdata$OS.time, survdata$OS) ~ U.cycle.group, data = survdata))
# Basic survival curves
level<-levels(factor(survdata$U.cycle.group))
c<-ifelse(level[1]=="Low risk",paste("blue"),paste("red"))
d<-ifelse(level[2]=="High risk",paste("red"),paste("blue"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
cellcycle <- ggsurvplot(fit, data = survdata, risk.table = F,
                        risk.table.height = 0.20,
                        risk.table.y.text = FALSE,
                        risk.table.title ="",
                        main = "Survival curve",
                        palette =c(c, d),pval=T, 
                        risk.table.y.text.col = T,
                        legend = c(0.85, 0.90),
                        legend.title = "",
                        xlab="Overall survival (days)",
                        legend.labs=c(e,f))
save(cellcycle,file="KMcellcycle.Rdata")
#
##################################################################################
#cytokine correlation analysis
expres.mRNA.cytokine<-expres.mRNA[expres.mRNA$entrezID %in% cytokine,]
expres.mRNA.cytokine<-as.data.frame(expres.mRNA.cytokine)
expres.mRNA.cytokine<-subset(expres.mRNA.cytokine,select=-entrezID)
source("D:/learning/methods/rfunction/metaDE.match/Match.gene.txt")
expres.mRNA.cytokine<-Match.gene(expres.mRNA.cytokine,"IQR")
expres.mRNA.cytokine<-as.data.frame(expres.mRNA.cytokine)


index<-intersect(names(expres.mRNA.cytokine),rownames(survdata))
expres.mRNA.cytokine<-expres.mRNA.cytokine[,index]
survdata<-survdata[index,]

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] ) }

data<-cbind(survdata$risk,as.data.frame(t(expres.mRNA.cytokine)))
library(Hmisc)
res <- rcorr(as.matrix(data),type="spearman")
res<-CorMatrix (res$r, res$P)
res<-res[which(res$row=="survdata$risk"),]
res<-na.omit(res)
#res<-res[res$p<0.05,]
res$Adjusted.P<-p.adjust(res$p,method = "bonferroni")
res<-res[res$Adjusted.P<0.05,]
write.csv(res,"./result/res.cytokine.csv")
library(ggplot2)
correplot.cytokine.bar<-ggplot(data=res)+geom_bar(aes(x=column,y=cor, fill=p), stat='identity')+
  coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+ 
  ylab("Correlation coefficient")+
  theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="black", size=48))+
  scale_y_continuous(expand=c(0, 0))+ 
  scale_x_discrete(expand=c(0,0))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
save(correplot.cytokine.bar,file="correplot.cytokine.bar.Rdata")

correplot.cytokine.point<-ggplot(data=res,aes(x=cor,y=column))+geom_point(aes(col=cor,size=p))+
  scale_color_gradient(low="blue",high="red")+
  labs(color="R square", size="P value")+theme_bw()+
  ylab("")+xlab("R square")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
        axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))
save(correplot.cytokine.point,file="correplot.cytokine.point")
library(stats)
library(CCA)
library(mvstats)
index<-intersect(res$column,row.names(expres.mRNA.cytokine))
expres.mRNA.cytokine<-expres.mRNA.cytokine[index,]
ca<-cancor(scale(t(expres.mRNA.cytokine)),scale(survdata$risk))
cancor.test(t(expres.mRNA.cytokine),survdata$risk,plot=T)

corcoef.test <- function(r, n, p, q, alpha = 0.1) {
  #r为相关系??? n为样本个??? 且n>p+q
  m <- length(r)
  Q <- rep(0, m)
  lambda <- 1
  for (k in m:1) {
    lambda <- lambda * (1 - r[k] ^ 2)
    #检验统计量
    Q[k] <- -log(lambda)   #检验统计量取对???
  }
  s <- 0
  i <- m
  for (k in 1:m) {
    Q[k] <- (n - k + 1 - 1 / 2 * (p + q + 3) + s) * Q[k]  #统计???
    chi <- 1 - pchisq(Q[k], (p - k + 1) * (q - k + 1))
    if (chi > alpha) {
      i <- k - 1
      break
    }
    s <- s + 1 / r[k] ^ 2
  }
  i  #显示输出结果 选用第几对典型变???
}

n<-corcoef.test(r=ca$cor,n=ncol(expres.mRNA.cytokine),p=nrow(expres.mRNA.cytokine),q=1)
U.cytokine <- as.matrix(scale(t(expres.mRNA.cytokine))) %*% ca$xcoef
V.cytokine <- as.matrix(scale(survdata$risk)) %*% ca$ycoef
cytokinedata<-data.frame(U.cytokine[,n],V.cytokine,survdata$Group)
names(cytokinedata)<-c("signature","risk", "Group")

library(ggpubr)
coplot.cytokine.CCA<- ggscatter(data=cytokinedata, y = "risk",x = "signature", 
                                color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                                add = "reg.line", conf.int = TRUE,
                                legend = c(0.15, 0.98),
                                legend.title = "",
                                add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                                cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                                cor.coeff.args = list(label.x =0,label.y=0.15, label.sep = "\n"),
                                xlab ="cytokine gene signature" , ylab = "Methylation risk")
save(coplot.cytokine.CCA,file="coplot.cytokine.CCA.Rdata")
#survival analysis of the total effect of cell cytokine gene on the survival KIRC patients
U.cytokine<-as.data.frame(U.cytokine)
if(identical(row.names(U.cytokine),row.names(survdata))) {
  survdata$U.cytokine.group<-ifelse(as.numeric(U.cytokine[,n])<=median(as.numeric(U.cytokine[,n])),"Low risk","High risk")
}
#
library(survminer)
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(survdata$OS.time, survdata$OS) ~ U.cytokine.group, data = survdata)
survdiff(Surv(survdata$OS.time, survdata$OS) ~ U.cytokine.group, data = survdata)
summary(coxph(Surv(survdata$OS.time, survdata$OS) ~ U.cytokine.group, data = survdata))
# Basic survival curves
level<-levels(factor(survdata$U.cytokine.group))
c<-ifelse(level[1]=="Low risk",paste("blue"),paste("red"))
d<-ifelse(level[2]=="High risk",paste("red"),paste("blue"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
cytokine <- ggsurvplot(fit, data = survdata, risk.table = F,
                       risk.table.height = 0.20,
                       risk.table.y.text = FALSE,
                       risk.table.title ="",
                       main = "Survival curve",
                       palette =c(c, d),pval=T, 
                       risk.table.y.text.col = T,
                       legend = c(0.85, 0.90),
                       legend.title = "",
                       xlab="Overall survival (days)",
                       legend.labs=c(e,f))
save(cytokine,file="KMcytokine.Rdata")

#################################################################################################
#apotosis correlation analysis
apotosis=KEGGPATHID2EXTID[['hsa04210']]
expres.mRNA.apotosis<-expres.mRNA[expres.mRNA$entrezID %in% apotosis,]
expres.mRNA.apotosis<-as.data.frame(expres.mRNA.apotosis)
expres.mRNA.apotosis<-subset(expres.mRNA.apotosis,select=-entrezID)
source("D:/learning/methods/rfunction/metaDE.match/Match.gene.txt")
expres.mRNA.apotosis<-Match.gene(expres.mRNA.apotosis,"IQR")
expres.mRNA.apotosis<-as.data.frame(expres.mRNA.apotosis)


index<-intersect(names(expres.mRNA.apotosis),rownames(survdata))
expres.mRNA.apotosis<-expres.mRNA.apotosis[,index]
survdata<-survdata[index,]

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] ) }

data<-cbind(survdata$risk,as.data.frame(t(expres.mRNA.apotosis)))
library(Hmisc)
res <- rcorr(as.matrix(data),type="spearman")
res<-CorMatrix (res$r, res$P)
res<-res[which(res$row=="survdata$risk"),]
res<-na.omit(res)
#res<-res[res$p<0.05,]
res$Adjusted.P<-p.adjust(res$p,method = "bonferroni")
res<-res[res$Adjusted.P<0.05,]
write.csv(res,"./result/res.apotosis.csv")
library(ggplot2)
correplot.apotosis.bar<-ggplot(data=res)+geom_bar(aes(x=column,y=cor, fill=p), stat='identity')+
  coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+ 
  ylab("Correlation coefficient")+
  theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="black", size=48))+
  scale_y_continuous(expand=c(0, 0))+ 
  scale_x_discrete(expand=c(0,0))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
save(correplot.apotosis.bar,file="correplot.apotosis.bar.Rdata")

correplot.apotosis.point<-ggplot(data=res,aes(x=cor,y=column))+geom_point(aes(col=cor,size=p))+
  scale_color_gradient(low="blue",high="red")+
  labs(color="R square", size="P value")+theme_bw()+
  ylab("")+xlab("R square")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
        axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))
save(correplot.apotosis.point,file="correplot.apotosis.point")
library(stats)
library(CCA)
library(mvstats)
index<-intersect(res$column,row.names(expres.mRNA.apotosis))
expres.mRNA.apotosis<-expres.mRNA.apotosis[index,]
ca<-cancor(scale(t(expres.mRNA.apotosis)),scale(survdata$risk))
cancor.test(t(expres.mRNA.apotosis),survdata$risk,plot=T)

corcoef.test <- function(r, n, p, q, alpha = 0.1) {
  #r为相关系??? n为样本个??? 且n>p+q
  m <- length(r)
  Q <- rep(0, m)
  lambda <- 1
  for (k in m:1) {
    lambda <- lambda * (1 - r[k] ^ 2)
    #检验统计量
    Q[k] <- -log(lambda)   #检验统计量取对???
  }
  s <- 0
  i <- m
  for (k in 1:m) {
    Q[k] <- (n - k + 1 - 1 / 2 * (p + q + 3) + s) * Q[k]  #统计???
    chi <- 1 - pchisq(Q[k], (p - k + 1) * (q - k + 1))
    if (chi > alpha) {
      i <- k - 1
      break
    }
    s <- s + 1 / r[k] ^ 2
  }
  i  #显示输出结果 选用第几对典型变???
}

n<-corcoef.test(r=ca$cor,n=ncol(expres.mRNA.apotosis),p=nrow(expres.mRNA.apotosis),q=1)
U.apotosis <- as.matrix(scale(t(expres.mRNA.apotosis))) %*% ca$xcoef
V.apotosis <- as.matrix(scale(survdata$risk)) %*% ca$ycoef
apotosisdata<-data.frame(U.apotosis[,n],V.apotosis,survdata$Group)
names(apotosisdata)<-c("signature","risk", "Group")

library(ggpubr)
coplot.apotosis.CCA<- ggscatter(data=apotosisdata, y = "risk",x = "signature", 
                                color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                                add = "reg.line", conf.int = TRUE,
                                legend = c(0.15, 0.98),
                                legend.title = "",
                                add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                                cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                                cor.coeff.args = list(label.x =0,label.y=0.15, label.sep = "\n"),
                                xlab ="apotosis gene signature" , ylab = "Methylation risk")
save(coplot.apotosis.CCA,file="coplot.apotosis.CCA.Rdata")
#survival analysis of the total effect of cell apotosis gene on the survival KIRC patients
U.apotosis<-as.data.frame(U.apotosis)
if(identical(row.names(U.apotosis),row.names(survdata))) {
  survdata$U.apotosis.group<-ifelse(as.numeric(U.apotosis[,n])<=median(as.numeric(U.apotosis[,n])),"Low risk","High risk")
}
#
library(survminer)
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(survdata$OS.time, survdata$OS) ~ U.apotosis.group, data = survdata)
survdiff(Surv(survdata$OS.time, survdata$OS) ~ U.apotosis.group, data = survdata)
summary(coxph(Surv(survdata$OS.time, survdata$OS) ~ U.apotosis.group, data = survdata))
# Basic survival curves
level<-levels(factor(survdata$U.apotosis.group))
c<-ifelse(level[1]=="Low risk",paste("blue"),paste("red"))
d<-ifelse(level[2]=="High risk",paste("red"),paste("blue"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
apotosis <- ggsurvplot(fit, data = survdata, risk.table = F,
                       risk.table.height = 0.20,
                       risk.table.y.text = FALSE,
                       risk.table.title ="",
                       main = "Survival curve",
                       palette =c(c, d),pval=T, 
                       risk.table.y.text.col = T,
                       legend = c(0.85, 0.90),
                       legend.title = "",
                       xlab="Overall survival (days)",
                       legend.labs=c(e,f))
save(apotosis,file="KMapotosis.Rdata")

# SCNA correlation analysis
rm(list=ls())
setwd("D:/learning/KIRC/methana")
SCNA<-read.csv("D:/learning/methods/SCNA/SCNA.csv",header=T,row.names=1)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
row.names(survdata)<-substr(row.names(survdata),1,12)
row.names(SCNA)<-gsub("-",".",rownames(SCNA))
index<-intersect(row.names(SCNA),rownames(survdata))
length(index)
SCNA<-SCNA[index,]
survdata<-survdata[index,]

CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] ) }
risk<-survdata$risk
data<-data.frame(risk,as.data.frame(SCNA))
data<-as.data.frame(data)
dim(data)
data1<-subset(data,select=c(risk,SCNA.Level))
data<-subset(data,select=c(risk,Chrom.SCNA.Level,Arm.SCNA.Level,Focal.SCNA.Level))
library(Hmisc)
res <- rcorr(as.matrix(data),type="spearman")
res<-CorMatrix (res$r, res$P)
res<-res[which(res$row=="risk"),]
res<-na.omit(res)
#res<-res[res$p<0.00001,]
res$Adjusted.P<-p.adjust(res$p,method = "bonferroni")
res<-res[res$Adjusted.P<0.05,]
write.csv(res,"./result/res.SCNA.csv")
library(ggplot2)
correplot.SCNA.bar<-ggplot(data=res)+geom_bar(aes(x=column,y=cor, fill=Adjusted.P), stat='identity')+
  coord_flip() +scale_fill_gradient(low="red", high = "blue") +xlab("")+ 
  ylab("Correlation coefficient")+
  theme(axis.text.x=element_text(color="black",size=48),axis.text.y=element_text(color="green", size=48),axis.title.y = element_text(size = 15,color = "green"))+
  scale_y_continuous(expand=c(0, 0))+ 
  scale_x_discrete(expand=c(0,0))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
save(correplot.SCNA.bar,file="correplot.SCNA.bar.Rdata")
correplot.SNA.point<-ggplot(data=res,aes(x=cor,y=column))+geom_point(aes(col=cor,size=Adjusted.P))+
  scale_color_gradient(low="blue",high="red")+
  labs(color="R square", size="Adjusted P")+theme_bw()+
  ylab("")+xlab("R square")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
        axis.text=element_text(size=12,color="black"),legend.text=element_text(size=10))
save(correplot.SNA.point,file="correplot.SCNA.bar.Rdata")
library(stats)
library(CCA)
library(mvstats)
x<-data[,2:4]
y<-data[,1]
ca<-cancor(scale(x),scale(y))
cancor.test(x,y,plot=T)

corcoef.test <- function(r, n, p, q, alpha = 0.1) {
  m <- length(r)
  Q <- rep(0, m)
  lambda <- 1
  for (k in m:1) {
    lambda <- lambda * (1 - r[k] ^ 2)
    #检验统计量
    Q[k] <- -log(lambda)   
  }
  s <- 0
  i <- m
  for (k in 1:m) {
    Q[k] <- (n - k + 1 - 1 / 2 * (p + q + 3) + s) * Q[k]  
    chi <- 1 - pchisq(Q[k], (p - k + 1) * (q - k + 1))
    if (chi > alpha) {
      i <- k - 1
      break
    }
    s <- s + 1 / r[k] ^ 2
  }
  i  
}

n<-corcoef.test(r=ca$cor,p=ncol(x),n=nrow(x),q=1)
U.SCNA <- as.matrix(scale(x)) %*% ca$xcoef
V.SCNA <- as.matrix(scale(y)) %*% ca$ycoef
SCNAdata<-data.frame(U.SCNA[,1],V.SCNA,survdata$Group)
names(SCNAdata)<-c("signature","risk", "Group")
library(ggpubr)
coplot.SCNA.CCA<- ggscatter(data=SCNAdata, y = "risk",x = "signature", 
                            color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                            add = "reg.line", conf.int = TRUE,
                            legend = c(0.15, 0.98),
                            legend.title = "",
                            add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                            cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                            cor.coeff.args = list(label.x =-0.2,label.y=0.15, label.sep = "\n"),
                            xlab ="SCNA" , ylab = "Methylation risk")
save(coplot.SCNA.CCA,file="coplot.SCNA.CCA.Rdata")
#survival analysis of the total effect of cell cycle gene on the survival KIRC patients
U.SCNA<-as.data.frame(U.SCNA)
row.names(U.SCNA)<-row.names(survdata)
if(identical(row.names(U.SCNA),row.names(survdata))) {
  survdata$U.SCNA.group<-ifelse(as.numeric(U.SCNA[,n])<=median(as.numeric(U.SCNA[,n])),"Low risk","High risk")
}
#
library(survminer)
#surv$Gender<-factor(surv$Gender,levels=c("1","2"), labels=c("1","2"))
fit<- survfit(Surv(survdata$OS.time, survdata$OS) ~ U.SCNA.group, data = survdata)
survdiff(Surv(survdata$OS.time, survdata$OS) ~ U.SCNA.group, data = survdata)
summary(coxph(Surv(survdata$OS.time, survdata$OS) ~ U.SCNA.group, data = survdata))
# Basic survival curves
level<-levels(factor(survdata$U.SCNA.group))
c<-ifelse(level[1]=="Low risk",paste("blue"),paste("red"))
d<-ifelse(level[2]=="High risk",paste("red"),paste("blue"))
e<-ifelse(level[1]== "Low risk",paste("Low risk group"),paste("High risk group"))
f<-ifelse(level[2]== "Low risk",paste("Low risk group"),paste("High risk group"))
SCNA<- ggsurvplot(fit, data = survdata, risk.table = F,
                  risk.table.height = 0.20,
                  risk.table.y.text = FALSE,
                  risk.table.title ="",
                  main = "Survival curve",
                  palette =c(c, d),pval=T, 
                  risk.table.y.text.col = T,
                  legend = c(0.85, 0.90),
                  legend.title = "",
                  xlab="Overall survival (days)",
                  legend.labs=c(e,f))
save(SCNA,file="KMSCNA.Rdata")
#total SCNA levels correlation
library(ggpubr)
data1<-as.data.frame(scale(data1))
data1$Group<-survdata$Group
coplot.SCNA.total<- ggscatter(data=data1, y = "risk",x = "SCNA.Level", 
                              color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                              add = "reg.line", conf.int = TRUE,
                              legend = c(0.15, 0.98),
                              legend.title = "",
                              add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                              cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                              cor.coeff.args = list(label.x =4,label.y=3.5, label.sep = "\n"),
                              xlab ="SCNA level" , ylab = "Methylation risk")
save(coplot.SCNA.total,file="coplot.SCNA.total.Rdata")
pdf("./result/coplot.SCNA.total.pdf")
coplot.SCNA.total
dev.off()
################################################################################################################
#Differentially analysis
rm(list=ls())
setwd("D:/learning/KIRC/methana")
library(edgeR)
expres<-read.table("D:/learning/KIRC/KIRC-HTSeq - Counts.merge.txt",header=T,row.names=1)
genecode<-read.table("D:/learning/methods/supportdata/gencode.v22.annotation.gene.txt",header=T)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
expres<-expres[,substr(colnames(expres), 14, 15) %in% c("01")]
index<-intersect(row.names(survdata),names(expres))
length(index)
expres<-expres[,index]
survdata<-survdata[index,]
identical(row.names(survdata),names(expres))
survdata<-survdata[order(survdata$risk,decreasing = T),]
expres<-expres[,row.names(survdata)]
identical(row.names(survdata),names(expres))
#differentially expression analysis using edgeR
group=survdata$Group
y <- DGEList(counts=expres, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#To perform quasi-likelihood F-tests:
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
summary(decideTests(lrt))

save(lrt,qlf,file="DEfit.Rdata")
DEG<-topTags(lrt,sort.by="logFC",n=nrow(expres))
DEG<-DEG$table
DEG$id<-row.names(DEG)
DEG.ann<-merge(genecode,DEG,by="id")
DEG.ann<-subset(DEG.ann,select=-c(chrom, chromStart,  chromEnd, strand))
imm.gene<-read.csv("D:/learning/methods/supportdata/immunegenescience.csv",header=T)
CD4.mature<-imm.gene$CD4.mature
imm.dif <- function(x) {
  levels <- levels(x)[-1]
  x <- factor(x, levels = levels)
  x <- x[complete.cases(x)]
  DEG.ann.x <- DEG.ann[which(DEG.ann$gene %in% x), ]
  DEG.ann.x <-DEG.ann.x[DEG.ann.x$FDR < 0.05 & abs(DEG.ann.x$logFC)>=0.5, ]
  return(DEG.ann.x)
}
& abs(DEG.ann$logFC)>=0.5
summary(DEG.ann$logFC)

CD4.mature.diff<-imm.dif(imm.gene$CD4.mature)
CD8.effector.diff<-imm.dif(imm.gene$CD8.effector)
NK.cells.diff<-imm.dif(imm.gene$NK.cells)
B.cells.diff<-imm.dif(imm.gene$B.cells)
T.reg.diff<-imm.dif(imm.gene$T.reg)
Dendritic.diff<-imm.dif(imm.gene$Dendritic)
CD8.effector.NK.cells.diff<-imm.dif(imm.gene$CD8.effector.NK.cells)
Macrophages.diff<-imm.dif(imm.gene$Macrophages)
Macrophages.M2.diff<-imm.dif(imm.gene$Macrophages.M2)
Macrophages.M1.diff<-imm.dif(imm.gene$Macrophages.M1)

expres$id<-row.names(expres)
expres<-merge(genecode,expres,by="id")
expres<-subset(expres,select=-c(id,chrom,chromStart,chromEnd,strand))
imm.exp <- function(x, y) {
  #x=CD4.mature.diff;y=expres
  imm.exp <- merge(x, y, by = "gene")
  row.names(imm.exp) <- imm.exp$gene
  imm.exp <- subset(imm.exp, select = -c(id, logFC, logCPM, LR, PValue, gene))
  return(imm.exp)
}

CD4.mature.exp<-imm.exp(CD4.mature.diff,expres)
CD8.effector.exp<-imm.exp(CD8.effector.diff,expres)
NK.cells.exp<-imm.exp(NK.cells.diff,expres)
B.cells.exp<-imm.exp(B.cells.diff,expres)
T.reg.exp<-imm.exp(T.reg.diff,expres)
Dendritic.exp<-imm.exp(Dendritic.diff,expres)
CD8.effector.NK.cells.exp<-imm.exp(CD8.effector.NK.cells.diff,expres)
Macrophages.exp<-imm.exp(Macrophages.diff,expres)
Macrophages.M2.exp<-imm.exp(Macrophages.M2.diff,expres)
Macrophages.M1.exp<-imm.exp(Macrophages.M1.diff,expres)

CD4.mature.exp$cell<-rep("CD4.mature",nrow(CD4.mature.exp))
CD8.effector.exp$cell<-rep("CD8.effector",nrow(CD8.effector.exp))
NK.cells.exp$cell<-rep("NK.cells",nrow(NK.cells.exp))
B.cells.exp$cell<-rep("B.cells",nrow(B.cells.exp))
T.reg.exp$cell<-rep("T.reg",nrow(T.reg.exp))
Dendritic.exp$cell<-rep("Dendritic",nrow(Dendritic.exp))
CD8.effector.NK.cells.exp$cell<-rep("CD8.effector.NK.cells",nrow(CD8.effector.NK.cells.exp))
Macrophages.exp$cell<-rep("Macrophages",nrow(Macrophages.exp))
Macrophages.M2.exp$cell<-rep("Macrophages.M2",nrow(Macrophages.M2.exp))
Macrophages.M1.exp$cell<-rep("Macrophages.M1",nrow(Macrophages.M1.exp))
heat<-rbind(CD4.mature.exp,CD8.effector.exp,NK.cells.exp,B.cells.exp,T.reg.exp,Dendritic.exp,CD8.effector.NK.cells.exp,
            Macrophages.exp,Macrophages.M2.exp,Macrophages.M1.exp)
heat<-heat[order(heat$cell),]
annot<-cbind(heat$FDR,heat$cell)
annot<-as.data.frame(annot)
names(annot)<-c("FDR","Cell")
row.names(annot)<-row.names(heatdata)
annot<-annot[order(annot$Cell),]
FDR1=log10(as.numeric(as.vector(annot$FDR)))
heatdata<-subset(heat,select=-c(FDR,cell))

identical(row.names(survdata),names(heatdata))
column_split = survdata$Group
normalize<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}
heatdata<-as.data.frame(t(scale(t(heatdata),center =T,scale=T)))
library(ComplexHeatmap)
#heatdata<-as.data.frame(t(normalize(t(heatdata))))
pdf("./result/heatmap2.pdf",width=15,height=16.5)
Heatmap(heatdata,border = F,column_title =c (paste(levels(survdata$Group)[1],"group"),paste(levels(survdata$Group)[2],"group")),
        column_title_gp = gpar(fill = c("red","blue"), col = "black",fontsize = 22),cluster_rows = FALSE,
        row_title =c (paste(levels(annot$Cell))),
        cluster_columns = F,name = "", show_column_names = F,column_split = survdata$Group,
        row_split = annot$Cell,right_annotation = rowAnnotation(log10FDR = anno_barplot(FDR1))
)
dev.off()
# ratios between different type of immune cells
rm(list=ls())
expres<-read.table("D:/learning/KIRC/TCGA-KIRC.htseq_fpkm.tsv/TCGA-KIRC.htseq_fpkm.tsv",header=T,row.names=1)
source("D:/learning/methods/rfunction/metaDE.match/Match.gene.txt")
genecode<-read.table("D:/learning/methods/supportdata/gencode.v22.annotation.gene.txt",header=T)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
expres<-expres[,substr(colnames(expres), 14, 16) %in% c("01A")]
names(expres)<-substr(names(expres),1,15)
index<-intersect(row.names(survdata),names(expres))
length(index)
#312
expres<-expres[,index]
survdata<-survdata[index,]
identical(row.names(survdata),names(expres))
survdata<-survdata[order(survdata$risk,decreasing = T),]
expres<-expres[,row.names(survdata)]
identical(row.names(survdata),names(expres))
imm.gene<-read.csv("D:/learning/methods/supportdata/immunegenescience.csv",header=T)
expres$id<-row.names(expres)
expres<-merge(genecode,expres,by="id")
expres<-subset(expres,select=-c(id,chrom,chromStart,chromEnd,strand))


imm.cell.gene <- function(x) {
  levels <- levels(x)[-1]
  x <- factor(x, levels = levels)
  x <- x[complete.cases(x)]
  DEG.ann.x <- expres[which(expres$gene %in% x), ]
  #DEG.ann.x <-DEG.ann.x[DEG.ann.x$FDR < 0.05 & abs(DEG.ann.x$logFC)>=0.5, ]
  DEG.ann.x<-Match.gene(DEG.ann.x,"IQR")
  return(DEG.ann.x)
}

CD4.mature<-imm.cell.gene(imm.gene$CD4.mature)
CD8.effector<-imm.cell.gene(imm.gene$CD8.effector)
NK.cells<-imm.cell.gene(imm.gene$NK.cells)
B.cells<-imm.cell.gene(imm.gene$B.cells)
T.reg<-imm.cell.gene(imm.gene$T.reg)
Dendritic<-imm.cell.gene(imm.gene$Dendritic)
CD8.effector.NK.cells<-imm.cell.gene(imm.gene$CD8.effector.NK.cells)
Macrophages<-imm.cell.gene(imm.gene$Macrophages)
Macrophages.M2<-imm.cell.gene(imm.gene$Macrophages.M2)
Macrophages.M1<-imm.cell.gene(imm.gene$Macrophages.M1)


CD4.mature.m<-colMeans(CD4.mature)
CD8.effector.m<-colMeans(CD8.effector)
NK.cells.m<-colMeans(NK.cells)
B.cells.m<-colMeans(B.cells)
T.reg.m<-colMeans(T.reg)
Dendritic.m<-colMeans(Dendritic)
CD8.effector.NK.cells.m<-colMeans(CD8.effector.NK.cells)
Macrophages.m<-colMeans(Macrophages)
Macrophages.M2.m<-colMeans(Macrophages.M2)
Macrophages.M1.m<-colMeans(Macrophages.M1)

library(ggpubr)
wt.ratio <- function( x,y) {
  #x<-CD8.effector.m;y<-T.reg.m
  require(ggpubr)
  ratio <- as.data.frame(x / y)
  ratio$group <- survdata$Group
  p <-
    ggboxplot(
      ratio.CD8.effector.Treg,
      x = "group",
      y = "ratio",
      fill = "group",
      xlab = FALSE,
      ylab = "Ratio between two cell types",
      font.y = 12,
      font.xtickslab = 12,
      font.ytickslab = 12
    ) + theme(legend.position = 'none') +
    stat_compare_means(
      method = "wilcox.test",
      #method = "t.test",
      label = "p.format",
      label.x.npc = 0.3,
      hide.ns = T
    )
  return(p)
}

t.ratio <- function(x, y) {
  #x<-CD8.effector.m;y<-T.reg.m
  require(ggpubr)
  ratio <- as.data.frame(x / y)
  ratio$group <- survdata$Group
  p <-
    ggboxplot(
      ratio.CD8.effector.Treg,
      x = "group",
      y = "ratio",
      fill = "group",
      xlab = FALSE,
      ylab = "Ratio between two cell types",
      font.y = 12,
      font.xtickslab = 12,
      font.ytickslab = 12
    ) + theme(legend.position = 'none') +
    stat_compare_means(
      #method = "wilcox.test",
      method = "t.test",
      label = "p.format",
      label.x.npc = 0.3,
      hide.ns = T
    )
  return(p)
}

wt.ratio(CD8.effector.m,T.reg.m)
wt.ratio(Macrophages.M1.m,Macrophages.M2.m)
wt.ratio(Macrophages.M1.m,T.reg.m)
t.ratio(CD8.effector.m,T.reg.m)
t.ratio(Macrophages.M1.m,Macrophages.M2.m)

ratio.CD8.effector.Treg<-as.data.frame(CD8.effector.m/T.reg.m)
ratio.CD8.effector.Treg$group<-survdata$Group
names(ratio.CD8.effector.Treg)<-c("ratio","group")
wilcox.test(ratio ~ group, data = ratio.CD8.effector.Treg)


p <- ggboxplot(ratio.CD8.effector.Treg, x="group", y="ratio",fill="group",
               xlab = FALSE,ylab = "ratio.CD8.effector.Treg",font.y=12,font.xtickslab = 12,
               font.ytickslab = 12)+
  theme(legend.position='none')+
  stat_compare_means(method = "wilcox.test",label= "p.format",label.x.npc=0.3,hide.ns=T)


#total mutation
rm(list=ls())
setwd("D:/learning/KIRC/methana")
library(maftools)
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
sample.barcode<-c()
for(i in 1:length(row.names(survdata))){
  sample.barcode[i]<-paste(row.names(survdata)[i],"A",sep = "")
}
survdata$Tumor_Sample_Barcode<-sample.barcode
surv<-survdata
row.names(surv)<-surv$Tumor_Sample_Barcode
maf<-read.maf(maf="D:/learning/KIRC/CGA.KIRC.varscan.ee42d944-ccbf-4406-9e7c-ffe1a0a4a1d7.DR-10.0.somatic.maf/TCGA.KIRC.varscan.somatic.maf",
              clinicalData=survdata)
total.mutation<-getSampleSummary(maf)
total.mutation<-as.matrix(total.mutation)
total.mutation<-as.data.frame(total.mutation)
row.names(total.mutation)<-substr(total.mutation$Tumor_Sample_Barcode,1,16)
row.names(total.mutation)<-gsub("-",".",row.names(total.mutation))
index<-intersect(row.names(total.mutation),row.names(surv))
length(index)
#[1] 255
surv<-as.data.frame(surv)
surv<-surv[index,]
total.mutation<-as.data.frame(total.mutation)
total.mutation<-total.mutation[index,]
#corplot
mutationdata<-as.data.frame(cbind(surv$risk,total.mutation$total))
mutationdata$Group<-surv$Group
names(mutationdata)<-c("risk","mutation","Group")
library(ggpubr)
coplot.mutation<- ggscatter(data=mutationdata, y = "risk",x = "mutation", 
                            color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                            add = "reg.line", conf.int = TRUE,
                            legend = c(0.15, 0.98),
                            legend.title = "",
                            add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                            cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                            cor.coeff.args = list(label.x =50,label.y=1.016, label.sep = "\n"),
                            xlab ="Number of mutation" , ylab = "Methylation risk")
save(coplot.mutation,file="coplot.mutation.Rdata")
pdf("./result/coplot.mutation.pdf")
coplot.mutation
dev.off()
write.mafSummary(maf = maf, basename = 'maf')
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = maf, top = 20, fontSize = 12)
oncostrip(maf = maf, top = 20)
maf.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = maf.titv)
#lollipopPlot(maf = maf, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)
#rainfallPlot(maf = maf, detectChangePoints = TRUE, pointSize = 0.6)
plotVaf(maf = maf)
geneCloud(input = maf, minMut = 3)
laml.sig = oncodrive(maf = maf, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
fab.ce = clinicalEnrichment(maf = maf, clinicalFeature = 'Group')
head(laml.sig)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)

expres$id<-row.names(expres)
expres<-merge(genecode,expres,by="id")
expres<-subset(expres,select=-c(id,chrom,chromStart,chromEnd,strand))
entrez2symbol<-read.table("D:/learning/methods/supportdata/entrz2symbol.txt",header=T,sep="\t")
names(expres)[1]<-"symbol"
expres.m<-merge(expres,entrez2symbol,by="symbol")
#####################################################################################################
#multiple plot genomic metrics
library(ggplot2)
load("coplot.mutation.Rdata")
load("coplot.MATH.Rdata")
load("Coplot.SCNA.total.Rdata")
pdf("./result/coplot.genomicmetrix.pdf",height = 16.5,width = 5.5)
ggarrange(
  coplot.mutation,
  coplot.SCNA.total,
  coplot.MATH,
  ncol = 1,
  nrow = 3,
  labels = "AUTO",
  align = "hv",
  font.label = list(face = "plain",family="Times")
)
dev.off()
#####################################################################################################
#Immune infltration
#####################################################################################################
#calculation score using ssGSEA  science article
rm(list=ls())
setwd("D:/learning/KIRC/methana")
expres<-read.table("D:/learning/KIRC/TCGA-KIRC.htseq_fpkm.tsv/TCGA-KIRC.htseq_fpkm.tsv",header=T,row.names=1)
expres<-expres[,substr(names(expres),14,16) %in% "01A"]
genecode<-read.table("D:/learning/methods/supportdata/gencode.v22.annotation.gene.txt",header=T)
entrez_symbol_ensembl<-read.table("D:/learning/methods/supportdata/entrez_symbol_ensembl.txt",header=T)
genes_Science<-read.csv("D:/learning/methods/supportdata/immunegenescience.csv",header=T)
names(genecode)
expres$id<-row.names(expres)
expres<-merge(genecode,expres,by="id")
expres<-subset(expres,select=-c(id,chrom,chromStart,chromEnd,strand))
source("D:/learning/methods/rfunction/metaDE.match/Match.gene.txt")
expres<-Match.gene(expres,"IQR")
save(expres,file="expres.Match.gene.Rdata")
#load("expres.Match.gene.Rdata")
names(genes_Science)

CD4.mature<-levels(genes_Science$CD4.mature)[-1]
CD4.mature %in% row.names(expres)

CD8.effector<-levels(genes_Science$CD8.effector)[-1]
CD8.effector %in% row.names(expres)

B.cells<-levels(genes_Science$B.cells)[-1]
B.cells %in% row.names(expres)

NK.cells<-levels(genes_Science$NK.cells)[-1]
NK.cells %in% row.names(expres)

T.reg<-levels(genes_Science$T.reg)[-1]
T.reg %in% row.names(expres)

Dendritic<-levels(genes_Science$Dendritic)[-1]
Dendritic %in% row.names(expres)

Macrophages<-levels(genes_Science$Macrophages)[-1]
Macrophages %in% row.names(expres)

CD8.effector.NK.cells<-levels(genes_Science$CD8.effector.NK.cells)[-1]
CD8.effector.NK.cells %in% row.names(expres)

Macrophages.M2<-levels(genes_Science$Macrophages.M2)[-1]
Macrophages.M2 %in% row.names(expres)

Macrophages.M1<-levels(genes_Science$Macrophages.M1)[-1]
Macrophages.M1 %in% row.names(expres)

IFN_Gamapathy<-levels(genes_Science$IFN_Gamapathy)[-1]
IFN_Gamapathy %in% row.names(expres)

Immue_response_cytokine_rich_microenvironment<-levels(genes_Science$Immue_response_cytokine_rich_microenvironment)[-1]
Immue_response_cytokine_rich_microenvironment %in% row.names(expres)

Proinflamatory<-levels(genes_Science$Proinflamatory)[-1]
Proinflamatory %in% row.names(expres)

Anti_inflamatory<-levels(genes_Science$Anti_inflamatory)[-1]
Anti_inflamatory %in% row.names(expres)

Angiogenesis<-levels(genes_Science$Angiogenesis)
Angiogenesis %in% row.names(expres)

PD1molecular<-levels(genes_Science$PD1molecular)[-1]
PD1molecular %in% row.names(expres)

Cytotoxic<-levels(genes_Science$Cytotoxic)[-1]
Cytotoxic %in% row.names(expres)

Immunosuppressive<-levels(genes_Science$Immunosuppressive)[-1]
Immunosuppressive %in% row.names(expres)

cytotoxic_T_lymphocyte<-levels(genes_Science$cytotoxic_T_lymphocyte)[-1]
cytotoxic_T_lymphocyte %in% row.names(expres)

Interferon_stimulated_chemokines<-levels(genes_Science$Interferon_stimulated_chemokines)[-1]
Interferon_stimulated_chemokines %in% row.names(expres)

names(genes_Science)
paste(names(genes_Science),collapse=",")
genelist <- list(
  CD4.mature,
  CD8.effector,
  NK.cells,
  B.cells,
  T.reg,
  Dendritic,
  CD8.effector.NK.cells,
  Macrophages,
  Macrophages.M2,
  Macrophages.M1,
  IFN_Gamapathy,
  Immue_response_cytokine_rich_microenvironment,
  Proinflamatory,
  Anti_inflamatory,
  Angiogenesis,
  PD1molecular,
  Cytotoxic,
  Immunosuppressive,
  cytotoxic_T_lymphocyte,
  Interferon_stimulated_chemokines
)

names(genelist)<-names(genes_Science)
library(GSVA)
#source("D:/learning/methods/rfunction/metaDE.match/Match.gene.txt")
#expres<-Match.gene(expres,"IQR")

#expres<-as.matrix(expres)# make sure the mode of expres is matrix
GSEAscore<-gsva(expres,genelist, annotation,method = "ssgsea", parallel.sz = 4, verbose = TRUE)
GSEAscore<-as.data.frame(GSEAscore)
write.csv(GSEAscore,"./result/GSEAscore.sci.csv")
#GSEAscore<-read.csv("./result/GSEAscore.sci.csv",header=T,row.names=1)
#correlation between methylation risk and angiogenesis
train.surv.risk<-read.csv("./result/train.surv.risk.csv",header=T,row.names=1)
test.surv.risk<-read.csv("./result/test.surv.risk.csv",header=T,row.names=1)
survdata<-rbind(train.surv.risk,test.surv.risk)
GSEAscore<-as.data.frame(GSEAscore)
names(GSEAscore)<-substr(names(GSEAscore),1,15)
index<-intersect(names(GSEAscore),row.names(survdata))
GSEAscore<-GSEAscore[,index]
survdata<-survdata[index,]
GSEAscore<-as.data.frame(t(GSEAscore))

coplot.ssGSEA.sci<-list()
library("ggpubr")
for(i in 1: ncol(GSEAscore)){
  data<-cbind(GSEAscore[,i],survdata$risk)
  data<-as.data.frame(data)
  names(data)<-c("GSEAscore","methylation")
  data$Group<-survdata$Group
  coplot.ssGSEA.sci[[i]]<- ggscatter(data, x = "GSEAscore", y = "methylation", 
                                     color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                                     add = "reg.line", conf.int = TRUE,
                                     legend = c(0.15, 0.15),
                                     legend.title = "",
                                     add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                                     cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                                     cor.coeff.args = list(label.x = 0.5,label.y=1, label.sep = "\n"),
                                     xlab = paste(names(GSEAscore)[i]), ylab = "Methylation risk")
}

coplot.ssGSEA.sci[[5]]
save(coplot.ssGSEA.sci,file="coplot.ssGSEA.science.Rdata")
#
##GSEA result in the science
GSEAscore.sci<-read.csv("./result/GSEAscore.sci.csv",header = T,row.names=1)
GSEAscore.sci<-as.data.frame(t(GSEAscore.sci))
row.names(GSEAscore.sci)<-substr(row.names(GSEAscore.sci),1,15)
index<-intersect(row.names(survdata),row.names(GSEAscore.sci))
GSEAscore.sci<-GSEAscore.sci[index,]
survdata<-survdata[index,]
proinf_antiprof<-(GSEAscore.sci$Proinflamatory-GSEAscore.sci$Anti_inflamatory)

data<-cbind(proinf_antiprof,survdata$risk)
data<-as.data.frame(data)
data$Group<-survdata$Group
names(data)<-c("signature","risk", "Group")

library(ggpubr)
coplot.proinf_antiprof.GSEA<- ggscatter(data=data, y = "risk",x = "signature", 
                                        color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                                        add = "reg.line", conf.int = TRUE,
                                        legend = c(0.15, 0.12),
                                        legend.title = "",
                                        add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                                        cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                                        cor.coeff.args = list(label.x =-0.15,label.y=1.016, label.sep = "\n"),
                                        xlab ="Proinflammatory cytokines / anti-inflammatory cytokines" , ylab = "Methylation risk")
save(coplot.proinf_antiprof.GSEA,file="coplot.proinf_antiprof.GSEA.Rdata")
pdf("./result/coplot.proinf_antiprof.GSEA.pdf",width=5.1,height=5.1)
coplot.proinf_antiprof.GSEA
dev.off()
# M1 vs M2
M1_M2<-(GSEAscore.sci$Macrophages.M1-GSEAscore.sci$Macrophages.M2)

data<-cbind(M1_M2,survdata$risk)
data<-as.data.frame(data)
data$Group<-survdata$Group
names(data)<-c("signature","risk", "Group")

library(ggpubr)
coplot.M1_M2.GSEA<- ggscatter(data=data, y = "risk",x = "signature", 
                              color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                              add = "reg.line", conf.int = TRUE,
                              legend = c(0.15, 0.12),
                              legend.title = "",
                              add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                              cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                              cor.coeff.args = list(label.x =-0.15,label.y=1.016, label.sep = "\n"),
                              xlab ="M1 macrophage / M2 macrophage" , ylab = "Methylation risk")
save(coplot.M1_M2.GSEA,file="coplot.M1_M2.GSEA.Rdata")
pdf("./result/coplot.M1_M2.GSEA.pdf",width=5.1,height=5.1)
coplot.M1_M2.GSEA
dev.off()
##########################################################################################################

Cytolytic.Activity<-GSEAscore.sci$Cytotoxic

data<-cbind(Cytolytic.Activity,survdata$risk)
data<-as.data.frame(data)
data$Group<-survdata$Group
names(data)<-c("signature","risk", "Group")

library(ggpubr)
coplot.Cytolytic.Activity.GSEA<- ggscatter(data=data, y = "risk",x = "signature", 
                                           color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                                           add = "reg.line", conf.int = TRUE,
                                           legend = c(0.15, 0.12),
                                           legend.title = "",
                                           add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                                           cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                                           cor.coeff.args = list(label.x =-0.15,label.y=1.016, label.sep = "\n"),
                                           xlab ="Cytolytic Activity" , ylab = "Methylation risk")
save(coplot.Cytolytic.Activity.GSEA,file="coplot.Cytolytic.Activity.GSEA.Rdata")
pdf("./result/coplot.Cytolytic.Activity.GSEA.pdf",width=5.1,height=5.1)
coplot.Cytolytic.Activity.GSEA
dev.off()

##
PD1<-GSEAscore.sci$PD1molecular

data<-cbind(PD1,survdata$risk)
data<-as.data.frame(data)
data$Group<-survdata$Group
names(data)<-c("signature","risk", "Group")

library(ggpubr)
coplot.PD1.GSEA<- ggscatter(data=data, y = "risk",x = "signature", 
                            color = "Group", size = 1,palette = c("#FC4E07","#00AFBB" ), # Points color, shape and size              
                            add = "reg.line", conf.int = TRUE,
                            legend = c(0.15, 0.12),
                            legend.title = "",
                            add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                            cor.coef = TRUE, cor.method = "spearman",cor.coef.size=4.5,
                            cor.coeff.args = list(label.x =0.45,label.y=1.016, label.sep = "\n"),
                            xlab ="PD1" , ylab = "Methylation risk")
save(coplot.PD1.GSEA,file="coplot.PD1.GSEA.Rdata")
pdf("./result/coplot.PD1.GSEA.pdf",width=5.1,height=5.1)
coplot.PD1.GSEA
dev.off()