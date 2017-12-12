##Load packages
library(optmatch)
#library(RItools)
library(party)
library(MASS)
#source("/Users/jrigdon/Box sync/Rigdon/Useful Functions/Tables.R")
#source("/Users/jrigdon/Box sync/Rigdon/Useful Functions/Functions.R")

##Must start with complete case data set

##Include functions smahal and formulize
smahal = function(X){
  X = as.matrix(X)
  n = dim(X)[1]
  k = dim(X)[2]
  for (j in 1:k) X[, j]=rank(X[, j]) #compute on ranks
  cv = cov(X)
  vuntied = var(1:n)
  rat = sqrt(vuntied/diag(cv))
  cv = diag(rat) %*% cv %*% diag(rat)
  out = matrix(NA, n, n)
  icov = ginv(cv)
  for (i in 1:n) out[i, ] = mahalanobis(X, X[i, ],icov, inverted=T)
  out
}

formulize = function(outc, covs) {
  string = paste( c(outc, paste(covs, collapse=" + ") ), collapse=" ~ " )
  return(as.formula(string))
}


##Execute pair-match with RCT data and create "combined" data set
matchd = function(dta, id, trt, outc, covs, ind.cat) {
    dta2 = dta[order(dta[, names(dta)==trt], decreasing=TRUE), ]
    trt2 = dta2[, names(dta2)==trt] #get trt variable
    id2 = dta2[, names(dta2)==id] #get ID variable
    dta2$tID = paste(trt2, id2, sep="_") #ID for match

    dta3 = dta2[, names(dta2) %in% covs] #data for match
    dist = smahal(X=dta3)
    rownames(dist) = dta2$tID
    colnames(dist) = dta2$tID
    dist2 = dist[substr(rownames(dist), 1, 1)==1, substr(colnames(dist), 1, 1)==0] #get distance matrix of treated (rows) by control (columns)
    if (dim(dist2)[1]>dim(dist2)[2]) {dist2 = t(dist2)} #transpose if more treated than controls

    options("optmatch_max_problem_size" = Inf) #accommodate larger trials
    pm = pairmatch(dist2) #execute pair-match
    org = data.frame(ID=names(pm), match=pm)

    ##Make into treated-control ID data frame
    org2 = org[order(org$match), ]
    org3 = org2[!is.na(org2$match), ]
    nmatch = dim(org3)[1]/2
    jj = data.frame(m1=rep(NA, nmatch), m2=rep(NA, nmatch))
    for (i in 1:nmatch) {
        jj$m1[i] = as.character(org3$ID[2*i-1])
        jj$m2[i] = as.character(org3$ID[2*i])
    }

    if (substr(jj[1,1], 1, 1)=="0") {
        names(jj) = c("con", "trt")
    } else if (substr(jj[1,1], 1, 1)=="1") {
        names(jj) = c("trt", "con")
    }

    ##Make delta data set
    trt = dta2[match(jj$trt, dta2$tID), ]
    t1 = trt[, names(trt) %in% covs]
    con = dta2[match(jj$con, dta2$tID), ]
    c1 = con[, names(con) %in% covs]
    aa = list(t1, c1)
    dd = data.frame(delta=trt[, names(trt)==outc]-con[, names(con)==outc], Reduce(`+`, aa) / length(aa))

    ##Remove the categorical variables where trt != con
    covs2 = covs[which(ind.cat==1)]
    t1cat = t1[, names(t1) %in% covs2]
    c1cat = c1[, names(c1) %in% covs2]
    elim = double()
    if (length(covs2)==1) {
        elim = as.numeric(t1cat!=c1cat)
    } else if (length(covs2)>1) {
        elim = apply(t1cat!=c1cat, 1, sum)
    }

    ##Return data frame and make sure categorical variables are factors
    dd2 = dd[elim==0, ]
    for (i in 1:length(covs2)) {
        dd2[, names(dd2)==covs2[i]] = as.factor(dd2[, names(dd2)==covs2[i]])
    }
    dd2
}



##Function for inference using ctree
##Make sure newd contains informative row names!  Only way to ID individuals
inf_delta = function(delt, newd) {
    ##Bind data to get all in same class
    newd$delta = NA
    newd2 = newd[, match(names(delt), names(newd))]
    all = rbind(delt, newd2)

    ##Make Table 1 of differences
#    all0 = all
#    all0$orig = ifelse(is.na(all0$delta), 1, 0)
#    all0 = all0[, names(all0)!="delta"]
#    disp = names(all0)[names(all0)!="orig"]
#    ind = double()
#    for (i in 1:length(disp)) {
#        ind = c(ind, class(all0[, names(all0)==disp[i]])=="factor"#)
#    }
#    tab1 = mktab(data=all0, var.names=disp, ind.cat=ind, group.name="orig", cfn=describeMean, miss="always", pval=FALSE, tot=FALSE, digit=1)
#    colnames(tab1) = c("Combined", "Original")

    ##Split into train and test
    train = all[!is.na(all$delta), ]
    train$delta = as.factor(train$delta) #make sure factor
    test = all[is.na(all$delta), ]
    test2 = test[, names(test)!="delta"]

    ##Estimate model on combined data set
    covs = names(test2)
    fit = ctree(formulize(outc="delta", covs=covs),  data=train)

    ##Inference for delta in initial model
    test2$predDelta = predict(fit, newdata=test2, type="response")
    test2$predNeg1 = sapply(predict(fit, newdata=test2, type="prob"),'[[',1)
    test2$predZero = sapply(predict(fit, newdata=test2, type="prob"),'[[',2)
    test2$predPos1 = sapply(predict(fit, newdata=test2, type="prob"),'[[',3)
    test2$predNode = predict(fit, newdata=test2, type="node")

    ##Return what we want
    #output.all = list(tab1=tab1, fit=fit, dta=test2)
    output.all = list(fit=fit, dta=test2)
    return(output.all)
}

##Evaluate all of the tables
evalT = function(tt) {

    ##Supplement with row of 0s for predictions that aren't there
    rr = rownames(tt)
    if (length(rr)==1 & rr[1]=="-1") {
        tt = rbind(tt, c(0,0,0), c(0,0,0))
    }    else if (length(rr)==1 & rr[1]=="0") {
        tt = rbind(c(0,0,0), tt, c(0,0,0))
    }    else if (length(rr)==1 & rr[1]=="1") {
        tt = rbind(c(0,0,0), c(0,0,0), tt)
    }    else if (length(rr)==2 & rr[1]=="-1" & rr[2]=="0") {
        tt = rbind(tt, c(0,0,0))
    }    else if (length(rr)==2 & rr[1]=="-1" & rr[2]=="1") {
        tt = rbind(tt[1, ], c(0,0,0), tt[2, ])
    }    else if (length(rr)==2 & rr[1]=="0" & rr[2]=="1") {
        tt = rbind(c(0,0,0), tt)
    }

    ##Do calculations
    rowS = apply(tt, 1, sum)
    colS = apply(tt, 2, sum)
    sens = tt[1,1] / colS[1]
    ppv = tt[1,1] / rowS[1]

    df = data.frame(true = (tt[1,1] + tt[2,2] + tt[3,3]) / sum(tt),
    sens = tt[1,1] / colS[1],
    spec = tt[3,3] / colS[3],
    ppv = tt[1,1] / rowS[1],
    npv = tt[3,3] / rowS[3],
    f1 = 2*(ppv*sens) / (ppv+sens),
    e1 = tt[1,3] / colS[3],
    e2 = tt[3,1] / rowS[3],
    e3 = tt[3,1] / colS[1],
    e4 = tt[1,3] / rowS[1]
                    )
    return(df)
}

##Send this file to Sherlock (must connect first)
#system(command="rsync -a --stats '/Users/jrigdon/Box Sync/Rigdon/Sanjay/HTE/HTE_v1.R' jrigdon@sherlock-dtn.stanford.edu:/home/jrigdon/HTE/HTE_v5.R")
