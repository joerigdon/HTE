library(Hmisc)

rmulti = function(prob, n) {
 as.numeric(rMultinom(probs=matrix(prob, 1, length(prob)), m=n))
}


##Set up data where B~bern(0.3), eGFR~norm(60, 20), age~norm(60,10), sex~bern(0.5)
true = data.frame(age=rnorm(6000, 68, 10), black=factor(rbinom(6000, 1, 0.3)), sbp=rnorm(6000, 140, 15), dbp=rnorm(6000, 78, 12), scr=rnorm(6000, 1.07, 0.34), eGFR=rnorm(6000, 72, 20), statin=factor(rbinom(6000, 1, 0.43)), aspirin=factor(rbinom(6000, 1, 0.51)), fram=rnorm(6000, 25, 12), smok=factor(rmulti(c(0.44, 0.42, 0.14), 6000)))


##Do the eGFR/aspirin 2-variable interaction
true$delta = NA #simulate delta
length(true$delta[true$aspirin==1 & true$eGFR<=72])
length(true$delta[true$aspirin==1 & true$eGFR>72])
length(true$delta[true$aspirin==0 & true$eGFR<=72])
length(true$delta[true$aspirin==0 & true$eGFR>72])


##(i) asp, eGFR<=72, ATE=-3D = 0.048, d~M(0.036, 0.88, 0.084)
true$delta[true$aspirin==1 & true$eGFR<=72] = rmulti(c(0.036, 0.88, 0.084), length(true$delta[true$aspirin==1 & true$eGFR<=72]))-2
mean(true$delta[true$aspirin==1 & true$eGFR<=72]) #0.049

##(ii) asp, eGFR>72, ATE=-7D = 0.112, d~M(0.004, 0.88, 0.116)
true$delta[true$aspirin==1 & true$eGFR>72] = rmulti(c(0.004, 0.88, 0.116), length(true$delta[true$aspirin==1 & true$eGFR>72]))-2
mean(true$delta[true$aspirin==1 & true$eGFR>72]) #0.117

##(iii) no asp, eGFR<=72, ATE=7D = -0.112, d~M(0.116, 0.88, 0.004)
true$delta[true$aspirin==0 & true$eGFR<=72] = rmulti(c(0.116, 0.88, 0.004), length(true$delta[true$aspirin==0 & true$eGFR<=72]))-2
mean(true$delta[true$aspirin==0 & true$eGFR<=72]) #-0.12

##(iv) no asp, eGFR>72, ATE=3D = -0.048, d~M(0.084, 0.88, 0.036)
true$delta[true$aspirin==0 & true$eGFR>72] = rmulti(c(0.084, 0.88, 0.036), length(true$delta[true$aspirin==0 & true$eGFR>72]))-2
mean(true$delta[true$aspirin==0 & true$eGFR>72]) #-0.049

##Look within groups
mean(true$delta) #-0.0008

true$Y1 = NA
true$Y0 = NA

true$Y1[true$delta==-1] = 0
true$Y0[true$delta==-1] = 1

true$Y1[true$delta==1] = 1
true$Y0[true$delta==1] = 0

true$Y1[true$delta==0] = 0
true$Y0[true$delta==0] = 0

sum((true$Y1-true$Y0)!=true$delta) #they all equal delta
mean(true$Y1-true$Y0) #-0.165 as specified above for ATE
mean(true$delta)
mean(true$Y1)
mean(true$Y0)


##Add in rowname
true$UNIQID = paste("A", rownames(true), sep="")

##Change types where necessary
true$delta = factor(true$delta)
true$Y1 = factor(true$Y1)
true$Y0 = factor(true$Y0)

##Save locally
write.csv(true, "/Users/jrigdon/Box Sync/Rigdon/Sanjay/HTE/DATA3A.csv", row.names=FALSE)


##Save potential outcome matrix to Sherlock
#system(command="rsync -a --stats '/Users/jrigdon/Box Sync/Rigdon/Sanjay/HTE/DATA1.csv' jrigdon@sherlock-dtn.stanford.edu:/home/jrigdon/HTE/DATA1.csv")







