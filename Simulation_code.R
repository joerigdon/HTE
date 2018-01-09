#Read in code,
source('/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/HTE_v1.R')
source('/Users/jrigdon/Box sync/Rigdon/Useful Functions/Tables.R')
source('/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/TablesASD.R')

#Read in true potential outcomes,
true = read.csv('/Users/jrigdon/Box Sync/Rigdon/Sanjay/HTE/DATA4.csv', header=TRUE)


#Take random sample of potential outcomes to get experiment,
set.seed(12)
ind.trt = sample(6000, 3000)
true$Z = 0
true$Z[ind.trt] = 1
true$Y = NA
true$Y[true$Z==1] = true$Y1[true$Z==1]
true$Y[true$Z==0] = true$Y0[true$Z==0]

##Look at ATE
r1 = mktabASD(data=true, var.names=c('black', 'eGFR', 'age', 'sex'), ind.cat=c(1, 0, 0, 1), group.name="Z", miss='always', digit=2)

#word.doc(obj.list=list(r1), obj.title='Table 1: Data at randomization', dest='/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Baseline_tab_2.docx', ftype='Arial', col.odd='white')

table(true$Y, true$Z)
mean(true$Y[true$Z==1])-mean(true$Y[true$Z==0]) #-0.11 (close to true ATE of -0.12)
prop.test(x=c(1317, 1648), n=c(3000, 3000))

#Now use logistic regression w/ interaction term,
lr = glm(Y ~ (black + eGFR + age + sex)*Z , family='binomial', data=true)
lr.aic = step(lr, trace = F, direction=c("backward"))
df1 = true[, c(1:4)]
df1$Z = 1
p1 = predict(lr.aic, df1, type='response')
df0 = true[, c(1:4)]
df0$Z = 0
p0 = predict(lr.aic, df0, type='response')
summary(p1-p0) #min -0.49; max 0.21
true$predLR = round(p1-p0)
t1 = table(true$predLR, true$delta)


#GRF
library(grf)
true2 = true
true3 = as.matrix(true2[, 1:4])
set.seed(42)
tau.forest = causal_forest(X=true3, Y=true2$Y, W=true2$Z)
tau.hat = predict(tau.forest, true3)
table(round(tau.hat$predictions))
true$predGRF = round(tau.hat$predictions)
t2 = table(true$predGRF, true$delta)


#Regular RF
gg = true[, names(true) %in% c("black", "eGFR", "age", "sex", "Z", "Y")]
gg$Y = factor(gg$Y)
library(randomForest)
set.seed(3)
rforest <- randomForest(Y~., data=gg)
true$predRF = as.numeric(as.character(predict(rforest, df1)))-as.numeric(as.character(predict(rforest, df0)))
t4 = table(true$predRF, true$delta)


#Execute match,
dd2 = matchd(dta=true, id='UNIQID', trt='Z', outc='Y', covs=c('black', 'eGFR', 'age', 'sex'), ind.cat=c(1, 0, 0, 1))


#mCART,
jj = inf_delta(delt=dd2, newd=true[, !names(true) %in% ('delta')])
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

covs = c('black', 'eGFR', 'age', 'sex')

fitLR = ctree(formulize(outc="predLR", covs=covs),  data=true)
plot(fitLR) #1 terminal node

fitGRF = ctree(formulize(outc="predGRF", covs=covs),  data=true)
plot(fitGRF) #12 terminal nodes

fitRF = ctree(formulize(outc="predRF", covs=covs),  data=true)
plot(fitRF) #33 terminal nodes


##Get nodes and look at ASD
true$nodeLR = predict(fitLR, newdata=true[, 1:4], type="node")
true$nodeGRF = predict(fitGRF, newdata=true[, 1:4], type="node")
true$nodeRF = predict(fitRF, newdata=true[, 1:4], type="node")

true$black = as.factor(true$black)
true$sex = as.factor(true$sex)
true$nodemC = predict(jj$fit, newdata=true[, 1:4], type="node")

##How many nodes within each method
table(true$nodeLR) #1
length(table(true$nodeRF)) #33
length(table(true$nodeGRF)) #12
table(true$nodemC) #4

##mCART
te1 = double()
asd1 = double()
tt1 = double()
true1 = double()
for (i in 1:length(table(true$nodemC))) {
dta = true[true$nodemC==names(table(true$nodemC))[i], ]
t1m = mktabASD(data=dta, var.names=c('black', 'eGFR', 'age', 'sex'), ind.cat=c(1, 0, 0, 1), group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te1 = c(te1, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)

asd1 = c(asd1, max(as.numeric(t1m[, 3]) , na.rm=TRUE))

tt1 = rbind(tt1, t1m) #table

true1 = c(true1, mean(dta$delta))
}
tt1
te1
asd1
true1
#plot(te1, asd1)


##LR
te2 = double()
asd2 = double()
tt2 = double()
true2 = double()
for (i in 1:length(table(true$nodeLR))) {
dta = true[true$nodeLR==names(table(true$nodeLR))[i], ]
t1m = mktabASD(data=dta, var.names=c('black', 'eGFR', 'age', 'sex'), ind.cat=c(1, 0, 0, 1), group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te2 = c(te2, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)

asd2 = c(asd2, max(as.numeric(t1m[, 3]) , na.rm=TRUE))

tt2 = rbind(tt2, t1m) #table

true2 = c(true2, mean(dta$delta))
}
tt2
te2
asd2
true2


##RF
te3 = double()
asd3 = double()
tt3 = double()
true3 = double()
for (i in 1:length(table(true$nodeRF))) {
dta = true[true$nodeRF==names(table(true$nodeRF))[i], ]
t1m = mktabASD(data=dta, var.names=c('black', 'eGFR', 'age', 'sex'), ind.cat=c(1, 0, 0, 1), group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te3 = c(te3, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)

asd3 = c(asd3, max(as.numeric(t1m[, 3]) , na.rm=TRUE))

tt3 = rbind(tt3, t1m) #table

true3 = c(true3, mean(dta$delta))
}
tt3
te3
asd3
true3


##GRF
te4 = double()
asd4 = double()
tt4 = double()
true4 = double()
for (i in 1:length(table(true$nodeGRF))) {
dta = true[true$nodeGRF==names(table(true$nodeGRF))[i], ]
t1m = mktabASD(data=dta, var.names=c('black', 'eGFR', 'age', 'sex'), ind.cat=c(1, 0, 0, 1), group.name="Z", miss="always", digit=1) #table
tt = table(dta$Y, dta$Z)
jr = prop.test(x=c(tt[2, 2], tt[2, 1]), n=c(apply(tt, 2, sum)[2:1]))
te4 = c(te4, jr$estimate[1]-jr$estimate[2]) #treatment effect (jr$conf.int[1] is conf int)

asd4 = c(asd4, max(as.numeric(t1m[, 3]) , na.rm=TRUE))

tt4 = rbind(tt4, t1m) #table

true4 = c(true4, mean(dta$delta))
}
tt4
te4
asd4
true4


##Record absolute biases
b1 = abs(te1-true1)
b2 = abs(te2-true2)
b3 = abs(te3-true3)
b4 = abs(te4-true4)


##Make figure of TEs vs ASD for each method
pdf("/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Plot_TE_ASD_v4.pdf")
plot(asd1, b1, xlim=c(0, 2.2), ylim=c(0, 0.45), col=1, pch=19, xlab="Max of ASD across all covariates within subgroup", ylab="Absolute value of bias of estimated subgroup treatment effect", )
points(asd2, b2, col=1, pch=2)
points(asd3, b3, col=1, pch=0)
points(asd4, b4, col=1, pch=4)
abline(v=0.2, lty=2)
legend('topright', legend=c('mCART', 'LR', 'RF', 'Gradient RF'), col=c(1, 1, 1, 1), pch=c(19, 2, 0, 4), bty='n', inset=0.03)
dev.off()

##Save tables
word.doc(obj.list=list(tt1, tt2, tt3, tt4), obj.title=c('Table 1: mCART', 'Table 2: LR', 'Table 3: RF', 'Table 4: GRF'), dest='/Users/jrigdon/Box sync/Rigdon/Sanjay/HTE/Tab_HTE_v6.docx', ftype='Arial', col.odd='white')

