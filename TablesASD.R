library(Gmisc, verbose=FALSE)
#install.packages("Hmisc")
library(Hmisc)
#install.packages("ReporteRs")
library(ReporteRs)

asd = function(v, Z, ind.cat) {
        v = as.numeric(as.character(v))
        trt = v[Z==1]
        ctrl = v[Z==0]
        m.trt = mean(trt, na.rm = TRUE)
        m.ctrl = mean(ctrl, na.rm = TRUE)
        sd2 = sqrt(var(trt, na.rm = TRUE)/2 + var(ctrl, na.rm = TRUE)/2)
        asd = abs(m.trt - m.ctrl)/sd2
        if (ind.cat==1) {
        asd = double()
        tt = as.numeric(names(table(v)))
        for (i in 1:length(tt)) {
         m.trt = sum(trt==tt[i])/length(trt)
         m.ctrl = sum(ctrl==tt[i])/length(ctrl)
         sd2 = sqrt(m.trt * (1 - m.trt)/2 + m.ctrl * (1 - m.ctrl)/2)
         asd = c(asd, abs(m.trt - m.ctrl)/sd2)
        }
    }
        asd
    }

mktabASD = function(data, var.names, ind.cat, group.name, miss, digit) {
    n = names(data)
    cols = data[, which(n==group.name)]
    r = table(data[, which(n==group.name)])
    j = c(apply(r, 1, function(x) paste("n=",x,sep="")))
    for (i in 1:length(var.names)) {
        if (ind.cat[i]==0) {
            outc = data[, which(n==var.names[i])]
            label(outc) = var.names[i]
            dd2 = getDescriptionStatsBy(outc, cols, html=TRUE, useNA=miss, statistics=FALSE, add_total_col=FALSE, continuous_fn=describeMean, digits=digit)
            dd2 = cbind(dd2, round(asd(outc, cols, 0), 2))
            rownames(dd2)[1] = var.names[i]
        }
        else if (ind.cat[i]==1) {
            outc = data[, which(n==var.names[i])]
            label(outc) = var.names[i]
            dd1 = getDescriptionStatsBy(factor(outc), cols, html=TRUE, useNA=miss, statistics=FALSE, add_total_col=FALSE, digits=digit)
            dd2 = rbind(rep("", dim(dd1)[2]), dd1)
            dd2 = cbind(dd2, c("", round(asd(outc, cols, 1), 2)))
            rownames(dd2)[1] = var.names[i]
        }
      j = rbind(j,dd2)
    }
    #Clean up a few things
    k = t(apply(j, 1, function(x) gsub("&plusmn;", "\\±", x)))
    k2 = t(apply(k, 1, function(x) gsub("&lt; ", "<", x)))
    rownames(k2)[rownames(k2)=="j"] = ""
    rownames(k2) = lapply(rownames(k2),function(x) gsub("a[.]", " ", x))
    rownames(k2) = lapply(rownames(k2),function(x) gsub("b[.]", " ", x))
    rownames(k2) = lapply(rownames(k2),function(x) gsub("c[.]", " ", x))
    rownames(k2) = lapply(rownames(k2),function(x) gsub("d[.]", " ", x))
    rownames(k2) = lapply(rownames(k2),function(x) gsub("e[.]", " ", x))
    rownames(k2) = lapply(rownames(k2),function(x) gsub("f[.]", " ", x))
    colnames(k2) = lapply(colnames(k2),function(x) gsub("a[.]", " ", x))
    colnames(k2) = lapply(colnames(k2),function(x) gsub("b[.]", " ", x))
    colnames(k2) = lapply(colnames(k2),function(x) gsub("c[.]", " ", x))
    colnames(k2) = lapply(colnames(k2),function(x) gsub("d[.]", " ", x))
    colnames(k2) = lapply(colnames(k2),function(x) gsub("e[.]", " ", x))
    colnames(k2) = lapply(colnames(k2),function(x) gsub("f[.]", " ", x))
    #Remove uninformative missing values if missing="ifany"
    if (miss=="always") {n0 = apply(k2, 1, function(x) sum(x=="0 (0%)" | x=="0 (0.0%)"))
    rmv = which(rownames(k2)=="Missing" & n0==2)
    if (length(rmv)>0) {k2 = k2[-as.numeric(rmv), ]}}
    k2[1, 3] = ""
    colnames(k2)[3] = "ASD"
    k2
}

#d = mtcars

#j1 = mktabASD(data=d, var.names=c("hp", "am", "gear"), ind.cat=c(0, 1, 1), group.name="vs", miss="always", digit=1)
#j1[1, ]


#asd(v=d$hp, Z=d$vs, ind.cat=0)
#asd(v=d$am, Z=d$vs, ind.cat=1)
#asd(v=d$gear, Z=d$vs, ind.cat=1) #Fix this one
