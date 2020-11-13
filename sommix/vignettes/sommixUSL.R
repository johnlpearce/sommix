## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load package-------------------------------------------------------------
library(sommix)

## ----dataset, echo=TRUE-------------------------------------------------------
data(dataset1)
summary(dataset1)

## ----datatrn,  echo=FALSE-----------------------------------------------------
data.trn<-dataset1[,3:9] 
summary(data.trn)

## ----datascale, echo=TRUE-----------------------------------------------------
data.trn.sc<-scale(data.trn, center=TRUE, scale=TRUE)
summary(data.trn.sc)

## ----vareval, echo=TRUE-------------------------------------------------------
vareval(data.trn)

## ----coreval, echo=TRUE-------------------------------------------------------
coreval(data.trn)

## ----pcaeval, echo=TRUE-------------------------------------------------------
pcaeval(data.trn.sc)

## ----grpeval, echo=TRUE-------------------------------------------------------
grpeval(data.trn.sc, kmx=30) 

## ----mapsize, echo=TRUE-------------------------------------------------------
mapsize(data.trn.sc, kmx=10)

## ----sommix, echo=TRUE--------------------------------------------------------
sommod<-sommix(data.trn.sc, somx=3, somy=3, seedopt="Y", nstarts=10) 
print(summary(sommod))

## ----somsumm, echo=TRUE-------------------------------------------------------
summ<-sommix_summ(som_obj=sommod)
names(summ)
head(summ$trn_data)

## ----somplots, echo=TRUE------------------------------------------------------
sommix_bar(som_obj=sommod, varnames=NULL, colormod=NULL,
                     nodelab=TRUE, labtype="XYs", labsize=1,
                     addFreq=TRUE, freqtype="frq", freqsize=0.75,
                     legsymsize=2, leglabsize=1, legtxtlas=2,
                     barstat="MED")

## ----somfitplots, echo=TRUE---------------------------------------------------
sommix_eval(sommod, labtype = "XYs", legcex=0.75)

## ----compeval, echo=TRUE------------------------------------------------------
compeval(summ)

