##Load optmatch
library(optmatch)
library(nearfar)
library(RItools)
library(party)
source("/Users/jrigdon/Box sync/Rigdon/Useful Functions/Tables.R")

##Load SPRINT
outc = read.csv("/Users/jrigdon/Box sync/Rigdon/Sanjay/SPRINT_Challenge/outcomes.csv", header=TRUE)
bl = read.csv("/Users/jrigdon/Box sync/Rigdon/Sanjay/SPRINT_Challenge/baseline.csv", header=TRUE)

##Only keep data we need for this illustration
table(outc$EVENT_PRIMARYORDEATH, exclude=NULL)

dta = merge(bl, outc[, names(outc) %in% c("MASKID", "EVENT_PRIMARYORDEATH")], by="MASKID", all.x=TRUE)

covs = c("AGE", "FEMALE", "RACE4", "BMI", "TRR", "HDL", "CHR", "GLUR", "UMALCR", "SCREAT", "SBP", "DBP", "SMOKE_3CAT", "ASPIRIN", "STATIN")

dta2 = dta[, names(dta) %in% c("MASKID", "EVENT_PRIMARYORDEATH", "INTENSIVE", covs)]
head(dta2)

tab1 = mktab(data=dta2, var.names=covs, ind.cat=c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1), group.name="INTENSIVE", cfn=describeMean, miss="always", pval=TRUE, tot="last", digit=1)

##Recode race
table(dta2$RACE4, exclude=NULL)
dta2$race = NA
dta2$race[dta2$RACE4=="WHITE"] = 1
dta2$race[dta2$RACE4=="BLACK"] = 2
dta2$race[dta2$RACE4=="HISPANIC"] = 3
dta2$race[dta2$RACE4=="OTHER"] = 4
table(dta2$RACE4, dta2$race, exclude=NULL)

##Only keep variables we need and complete cases
dta3 = dta2[, names(dta2)!="RACE4"]
dim(dta3) #9361
dta4 = dta3[complete.cases(dta3), ]
dim(dta4) #8764
table(dta4$INTENSIVE)

##Create rank-based Mahalanobis distance matrix
dta4 = dta4[order(dta4$INTENSIVE, decreasing=TRUE), ]
dta4$ID = paste(dta4$INTENSIVE, dta4$MASKID, sep="_")
dta5 = dta4[ !names(dta4) %in% c("INTENSIVE", "EVENT_PRIMARYORDEATH", "MASKID", "ID")]

dist = smahal(X=dta5)
dist[1:5, 1:5]
summary(as.numeric(dist))
hist(as.numeric(dist)) #bell-shaped/Normal?

rownames(dist) = dta4$ID
colnames(dist) = dta4$ID

dist2 = dist[substr(rownames(dist), 1, 1)==1, substr(colnames(dist), 1, 1)==0]
dist2[1:5, 1:5]
summary(as.numeric(dist2))
hist(as.numeric(dist2)) #still bell-shaped?

##Execute pair-match
dim(dist2)
dist3 = t(dist2) #if more rows than columns, transpose
dist3[1:5, 1:5]
dim(dist3)

options("optmatch_max_problem_size" = Inf)
pm = pairmatch(dist3)
summary(pm) #4364 matched pairs
org = data.frame(ID=names(pm), match=pm)
dim(org)
table(org$match, exclude=NULL)
org2 = org[org$match), ]
dim(org2) #4364*2 = 8728

org3 = org2[order(org2$match), ]

ind = data.frame(IDt=rep(NA, dim(org3)[1]/2), IDc=rep(NA, dim(org3)[1]/2))

for (i in 1:dim(ind)[1]) {
    ind$IDt[i] = as.character(org3$ID[2*i])
    ind$IDc[i] = as.character(org3$ID[2*i-1])
}

##Now, remake data frame for what we need later
##Data frame of treatment and control IDs
trt = dta4
trt = dta4[match(ind$IDt, dta4$ID), ]
trt2 = trt[, names(trt) %in% covs2]

con = dta4[match(ind$IDc, dta4$ID), ]
con2 = con[, names(con) %in% covs2]

##Categorical variables are SMOKE_3CAT, ASPIRIN, FEMALE, STATIN, and race
##Only keep treatment and control data where exact match on categorical variables
keep = which(trt2$SMOKE_3CAT==con2$SMOKE_3CAT & trt2$ASPIRIN==con2$ASPIRIN & trt2$FEMALE==con2$FEMALE & trt2$STATIN==con2$STATIN & trt2$race==con2$race) #still 2452 matches

trt0 = trt[keep, ]
con0 = con[keep, ]

inf = data.frame(delta=trt0$EVENT_PRIMARYORDEATH-con0$EVENT_PRIMARYORDEATH, AGE=(trt0$AGE+con0$AGE)/2, FEMALE=(trt0$FEMALE+con0$FEMALE)/2, race=(trt0$race+con0$race)/2, BMI=(trt0$BMI+con0$BMI)/2, TRR=(trt0$TRR+con0$TRR)/2, HDL=(trt0$HDL+con0$HDL)/2, CHR=(trt0$CHR+con0$CHR)/2, GLUR=(trt0$GLUR+con0$GLUR)/2, UMALCR=(trt0$UMALCR+con0$UMALCR)/2, SCREAT=(trt0$SCREAT+con0$SCREAT)/2, SBP=(trt0$SBP+con0$SBP)/2, DBP=(trt0$DBP+con0$DBP)/2, SMOKE_3CAT=(trt0$SMOKE_3CAT+con0$SMOKE_3CAT)/2, ASPIRIN=(trt0$ASPIRIN+con0$ASPIRIN)/2, STATIN=(trt0$STATIN+con0$STATIN)/2)

table(inf$delta)


tab2 = mktab(data=inf, var.names=c("AGE", "FEMALE", "race", "BMI", "TRR", "HDL", "CHR", "GLUR", "UMALCR", "SCREAT", "SBP", "DBP", "SMOKE_3CAT", "ASPIRIN", "STATIN"), ind.cat=c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1), group.name="delta", cfn=describeMean, miss="always", pval=FALSE, tot=FALSE, digit=1)

##Now use conditional trees or multinomial regression for inference
#library(rpart)
#install.packages("party")
inf$FEMALE = as.factor(inf$FEMALE)
inf$race = as.factor(inf$race)
inf$SMOKE_3CAT = as.factor(inf$SMOKE_3CAT)
inf$ASPIRIN = as.factor(inf$ASPIRIN)
inf$STATIN = as.factor(inf$STATIN)

fit = ctree(as.factor(delta) ~ AGE + FEMALE + race + BMI + TRR + HDL + CHR + GLUR + UMALCR + SCREAT + SBP + DBP + SMOKE_3CAT + ASPIRIN + STATIN,  data=inf)

#fit = ctree(as.factor(delta) ~ AGE + as.factor(FEMALE) + as.factor(race) + BMI + TRR + HDL + CHR + GLUR + UMALCR + SCREAT + SBP + DBP + as.factor(SMOKE_3CAT) + as.factor(ASPIRIN) + as.factor(STATIN),  data=inf) #5 nodes in decision tree
pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/SPRINT_HTE.pdf", width=8, height=8)
plot(fit)
dev.off()

plot(fit, type="simple")


## predict delta for original data
tab3 = mktab(data=dta4, var.names=c("AGE", "FEMALE", "race", "BMI", "TRR", "HDL", "CHR", "GLUR", "UMALCR", "SCREAT", "SBP", "DBP", "SMOKE_3CAT", "ASPIRIN", "STATIN"), ind.cat=c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1), group.name="EVENT_PRIMARYORDEATH", cfn=describeMean, miss="always", pval=FALSE, tot=FALSE, digit=1)


#sel = c("AGE", "race", "SMOKE_3CAT", "SCREAT")
dt_test = dta4[dta4$SMOKE_3CAT!=4, names(dta4) %in% covs2]
dt_test$FEMALE = as.factor(dt_test$FEMALE)
dt_test$race = as.factor(dt_test$race)
dt_test$SMOKE_3CAT = as.factor(dt_test$SMOKE_3CAT)
#levels(dt_test$SMOKE_3CAT) = union(
#levels(inf$SMOKE_3CAT)
dt_test$ASPIRIN = as.factor(dt_test$ASPIRIN)
dt_test$STATIN = as.factor(dt_test$STATIN)
mm = match(names(inf), names(dt_test))
mm2 = mm[!is.na(mm)]
dt_test2 = dt_test[, mm2]
dt_test2$AGE = as.numeric(dt_test2$AGE)
dt_test2$TRR = as.numeric(dt_test2$TRR)
dt_test2$HDL = as.numeric(dt_test2$HDL)
dt_test2$CHR = as.numeric(dt_test2$CHR)
dt_test2$GLUR = as.numeric(dt_test2$GLUR)
dt_test2$SBP = as.numeric(dt_test2$SBP)
dt_test2$DBP = as.numeric(dt_test2$DBP)

str(inf)
str(dt_test2)

dt_test2$predClass = predict(fit, newdata=dt_test2, type="response")    # obtain the class (0/1) #won't work because all levels of factor don't appear in inf? #All most likely to be 0s - no change
dt_test2$predNeg1 = sapply(predict(fit, newdata=dt_test2, type="prob"),'[[',1)  # obtain probability of class 1 (second element from the lists)
dt_test2$predZero = sapply(predict(fit, newdata=dt_test2, type="prob"),'[[',2)
dt_test2$predPos1 = sapply(predict(fit, newdata=dt_test2, type="prob"),'[[',3)
dt_test2$predNode = predict(fit, newdata=dt_test2, type="node")   # obtain the predicted node (in case you need it)

##Save tables
word.doc(obj.list=list(tab1, tab3, tab2), obj.title=c("Table 1: Demographic and risk factors by treatment group in SPRINT", "Table 2: Demographic and risk factors by outcome for matched participants", "Table 3: Demographic and risk factors by estimated treatment effect for combined SPRINT participants"), dest="/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Tables_2017-07-27.docx", ftype="Arial", col.odd="white")


















