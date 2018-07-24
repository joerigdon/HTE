##Read in code
#install.packages('survival')
#library(survival)
#install.packages('modeltools')
#library(modeltools)
#install.packages('mvtnorm')
#library(mvtnorm)
#library(multcomp)
#library(vcd)
#install.packages('coin', dependencies=TRUE)
#library(coin)
#install.packages('party', dependencies=TRUE)
#library(party)
#install.packages('randomForest')
#library(randomForest)

##Simple function
NC = function(x) {
    as.numeric(as.character(x))
}

source('/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/HTE_v1.R')
source('/Users/jrigdon/Box sync/Rigdon/Useful Functions/Tables.R')
source('/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/TablesASD.R')

#Read in true potential outcomes,
true = read.csv('/Users/jrigdon/Box Sync/Rigdon/Sanjay/HTE/DATA1B.csv', header=TRUE, colClasses=c("numeric", rep("factor", 8), "character"))

##Name covs and categories
covs=c('age', 'stage4', 'site', 'prevTrt', 'ecog', 'diseaseFree')
ind.cat=c(0, 1, 1, 1, 1, 1)

#Take random sample of potential outcomes to get experiment,
set.seed(12)
ind.trt = sample(200, 100)
true$Z = 0
true$Z[ind.trt] = 1
true$Y = NA
true$Y[true$Z==1] = true$Y1[true$Z==1]
true$Y[true$Z==0] = true$Y0[true$Z==0]
true$Y = true$Y-1

##Look at ATE
r1 = mktabASD(data=true, var.names=c('age', 'stage4', 'site', 'prevTrt', 'ecog', 'diseaseFree'), ind.cat=c(0, 1, 1, 1, 1, 1), group.name="Z", miss='always', digit=2)

#word.doc(obj.list=list(r1), obj.title='Table 1: Data at randomization', dest='/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Baseline_tab_2.docx', ftype='Arial', col.odd='white')

table(true$Y, true$Z)
mean(true$Y[true$Z==1])-mean(true$Y[true$Z==0]) #-0.27 (close to true ATE of -0.2)
prop.test(x=c(29, 56), n=c(100, 100))

#Now use logistic regression w/ interaction term,
lr = glm(Y ~ (age + stage4 + site + prevTrt + ecog + diseaseFree)*Z , family='binomial', data=true)
lr.aic = step(lr, trace = F, direction=c("backward"))
df1 = true[, names(true) %in% c('age', 'stage4', 'site', 'prevTrt', 'ecog', 'diseaseFree')]
df1$Z = 1
p1 = predict(lr.aic, df1, type='response')
df0 = true[, names(true) %in% c('age', 'stage4', 'site', 'prevTrt', 'ecog', 'diseaseFree')]
df0$Z = 0
p0 = predict(lr.aic, df0, type='response')
summary(p1-p0) #min -0.49; max 0.21
true$predLR1 = p1-p0
true$predLR2 = round(p1-p0)
t1 = table(true$predLR2, true$delta)


#GRF
library(grf)
true2 = true
true3 = apply(as.matrix(true2[, names(true2) %in% c('age', 'stage4', 'site', 'prevTrt', 'ecog', 'diseaseFree')]), 2, function(x) as.numeric(x))
set.seed(42)
tau.forest = causal_forest(X=true3, Y=true2$Y, W=true2$Z)
tau.hat = predict(tau.forest, true3)
table(round(tau.hat$predictions))
true$predGRF1 = tau.hat$predictions
true$predGRF2 = round(tau.hat$predictions)
t2 = table(true$predGRF2, true$delta)

#Regular RF
gg = true[, names(true) %in% c('age', 'stage4', 'site', 'prevTrt', 'ecog', 'diseaseFree', 'Z', 'Y')]
gg$Y = factor(gg$Y)
library(randomForest)
set.seed(3)
rforest = randomForest(Y~., data=gg)
true$predRF = NC(predict(rforest, df1))-NC(predict(rforest, df0))
t4 = table(true$predRF, true$delta)


#Execute match,
dd2 = matchd(dta=true, id='UNIQID', trt='Z', outc='Y', covs=covs, ind.cat=ind.cat)


#mCART,
##Make sure levels line up
str(dd2)
newdata1 = true[, names(true) %in% covs]
str(newdata1)
levels(dd2$stage4) = union(levels(dd2$stage4), levels(newdata1$stage4))
str(dd2)
table(dd2$stage4)

jj = inf_delta(delt=dd2, newd=newdata1)
#pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Fig_sim4_v2.pdf")
#plot(jj$fit)
#dev.off()
true$predmCART = jj$dta$predDelta
t6 = table(true$predmCART, true$delta)


#Evaluate all
all2 = round(rbind(evalT(t1), evalT(t2), evalT(t4), evalT(t6)), 3)
rownames(all2) = c('LR (AIC)', 'GRF', 'RF', 'mCART')
all2


##Look at the simple balance tables in the identified subgroups
#pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Fig_sim8_v2.pdf")
#plot(jj$fit) #4 terminal nodes
#dev.off()
fitLR1 = ctree(formulize(outc="predLR1", covs=covs),  data=true)
pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/LR_1B.pdf")
plot(fitLR1, type="simple") #4 terminal nodes
dev.off()

#fitLR2 = ctree(formulize(outc="predLR2", covs=covs),  data=true)
#plot(fitLR2) #1 terminal node

fitGRF1 = ctree(formulize(outc="predGRF1", covs=covs),  data=true)
pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/GRF_1B.pdf")
plot(fitGRF1, type="simple") #10 terminal nodes
dev.off()

#fitGRF2 = ctree(formulize(outc="predGRF2", covs=covs),  data=true)
#plot(fitGRF2) #12 terminal nodes

fitRF = ctree(formulize(outc="predRF", covs=covs),  data=true)
pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/RF_1B.pdf")
plot(fitRF, type="simple") #33 terminal nodes
dev.off()

pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/mCART_1B.pdf")
plot(jj$fit, type="simple") #How to control this plot??
dev.off()

##Get nodes and look at ASD

true$nodeLR1 = predict(fitLR1, newdata=newdata1, type="node")
#true$nodeLR2 = predict(fitLR2, newdata=newdata1, type="node")
true$nodeGRF1 = predict(fitGRF1, newdata=newdata1, type="node")
#true$nodeGRF2 = predict(fitGRF2, newdata=newdata1, type="node")
true$nodeRF = predict(fitRF, newdata=newdata1, type="node")
true$nodemC = jj$dta$predNode



##How many nodes within each method
length(table(true$nodeLR1)) #4
#length(table(true$nodeLR2)) #1
length(table(true$nodeRF)) #6
length(table(true$nodeGRF1)) #10
#length(table(true$nodeGRF2)) #2
length(table(true$nodemC)) #1

##mCART
te1 = double()
l1 = double()
u1 = double()
asd1 = double()
tt1 = double()
true1 = double()
for (i in 1:length(table(true$nodemC))) {
dta = true[true$nodemC==names(table(true$nodemC))[i], ]
t1m = mktabASD(data=dta, var.names=covs, ind.cat=ind.cat, group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te1 = c(te1, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)
l1 = c(l1, jr$conf.int[1])
u1 = c(u1, jr$conf.int[2])

asd1 = c(asd1, max(as.numeric(t1m[, 3]) , na.rm=TRUE))

tt1 = rbind(tt1, t1m) #table

true1 = c(true1, mean(NC(dta$delta)))
}
tt1
te1
l1
u1
asd1
true1
#plot(te1, asd1)


##LR
te2 = double()
l2 = double()
u2 = double()
asd2 = double()
tt2 = double()
true2 = double()
for (i in 1:length(table(true$nodeLR1))) {
dta = true[true$nodeLR1==names(table(true$nodeLR1))[i], ]
t1m = mktabASD(data=dta, var.names=covs, ind.cat=ind.cat, group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te2 = c(te2, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)
l2 = c(l2, jr$conf.int[1])
u2 = c(u2, jr$conf.int[2])

asd2 = c(asd2, max(as.numeric(t1m[, 3]) , na.rm=TRUE))

tt2 = rbind(tt2, t1m) #table

true2 = c(true2, mean(NC(dta$delta)))
}
tt2
te2
l2
u2
asd2
true2


##RF
te3 = double()
l3 = double()
u3 = double()
asd3 = double()
tt3 = double()
true3 = double()
for (i in 1:length(table(true$nodeRF))) {
dta = true[true$nodeRF==names(table(true$nodeRF))[i], ]
t1m = mktabASD(data=dta, var.names=covs, ind.cat=ind.cat, group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te3 = c(te3, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)
l3 = c(l3, jr$conf.int[1])
u3 = c(u3, jr$conf.int[2])

asd3 = c(asd3, max(as.numeric(t1m[, 3]) , na.rm=TRUE))

tt3 = rbind(tt3, t1m) #table

true3 = c(true3, mean(NC(dta$delta)))
}
tt3
te3
l3
u3
asd3
true3


##GRF
te4 = double()
l4 = double()
u4 = double()
asd4 = double()
tt4 = double()
true4 = double()
for (i in 1:length(table(true$nodeGRF1))) {
dta = true[true$nodeGRF1==names(table(true$nodeGRF1))[i], ]
t1m = mktabASD(data=dta, var.names=covs, ind.cat=ind.cat, group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)

if (dim(tt)[1]==2 & dim(tt)[2]==2) {
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te4 = c(te4, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)
l4 = c(l4, jr$conf.int[1])
u4 = c(u4, jr$conf.int[2])

} else if (dim(tt)[1]!=2 | dim(tt)[2]!=2) {
te4 = c(te4, NA)
l4 = c(l4, NA)
u4 = c(u4, NA)
}

asd4 = c(asd4, max(as.numeric(t1m[, 3]) , na.rm=TRUE))
tt4 = rbind(tt4, t1m) #table
true4 = c(true4, mean(NC(dta$delta)))

}
tt4
te4
l4
u4
asd4
true4

##Look at treatment effects for each method
te.mCART = data.frame(eff=paste(paste(paste(paste(round(te1, 2), "(", sep=" "), round(l1, 2), sep=""), round(u1, 2), sep=", "), ")", sep=""))
te.LR = data.frame(eff=paste(paste(paste(paste(round(te2, 2), "(", sep=" "), round(l2, 2), sep=""), round(u2, 2), sep=", "), ")", sep=""))
te.RF = data.frame(eff=paste(paste(paste(paste(round(te3, 2), "(", sep=" "), round(l3, 2), sep=""), round(u3, 2), sep=", "), ")", sep=""))
te.GRF = data.frame(eff=paste(paste(paste(paste(round(te4, 2), "(", sep=" "), round(l4, 2), sep=""), round(u4, 2), sep=", "), ")", sep=""))

te.mCART
te.LR
te.RF
te.GRF

##Record absolute biases
b1 = abs(te1-true1)
b2 = abs(te2-true2)
b3 = abs(te3-true3)
b4 = abs(te4-true4)


##Make figure of TEs vs ASD for each method
pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Plot_TE_ASD_1B.pdf")
plot(asd1, b1, xlim=c(0, 2.2), ylim=c(0, 0.45), col=1, pch=19, xlab="Max of ASD across all covariates within subgroup", ylab="Absolute value of bias of estimated subgroup treatment effect", )
points(asd2, b2, col=1, pch=2)
points(asd3, b3, col=1, pch=0)
points(asd4, b4, col=1, pch=4)
abline(v=0.2, lty=2)
legend('topright', legend=c('mCART', 'LR', 'RF', 'Gradient RF'), col=c(1, 1, 1, 1), pch=c(19, 2, 0, 4), bty='n', inset=0.03)
dev.off()

##Save tables
word.doc(obj.list=list(tt1, tt2, tt3, tt4), obj.title=c('Table 1: mCART', 'Table 2: LR', 'Table 3: RF', 'Table 4: GRF'), dest='/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Tab_HTE_1B.docx', ftype='Arial', col.odd='white')




