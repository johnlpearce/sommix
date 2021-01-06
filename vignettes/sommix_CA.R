## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load package-------------------------------------------------------------
library(sommix)

## ----dataset, echo=TRUE-------------------------------------------------------
data(dataset3)
summary(dataset3)

## ----datatrn,  echo=FALSE-----------------------------------------------------
data.x<-dataset3 
summary(data.x)

## ----vareval, echo=TRUE, fig.height=10, fig.width=7.5-------------------------
vareval(data.x)

## ----coreval, echo=TRUE, fig.height=10, fig.width=7.5-------------------------
coreval(data.x)

## ----pcaeval, echo=TRUE, fig.height=10, fig.width=7.5-------------------------
pcaeval(scale(data.x))

## ----datascale, echo=TRUE-----------------------------------------------------
data.x.ln<-log(data.x)
data.x.sc<-scale(data.x.ln, center=TRUE, scale=TRUE)
summary(data.x.sc)

## ----grpeval, echo=TRUE, fig.height=10, fig.width=7.5-------------------------
grpeval(data.x.sc, kmx=50) 

## ----mapsize, echo=TRUE,fig.height=10, fig.width=7.5--------------------------
mapsize(data.x.sc, kmx=50)

## ----sommix, echo=TRUE--------------------------------------------------------
sommod<-sommix(data.x.sc, somx=4, somy=3, seedopt="Y", nstarts=5, maptopo="hexagonal") 

## ----somsumm, echo=TRUE-------------------------------------------------------
summ<-sommix_summ(som_obj=sommod)
names(summ)
#Frequency table
summ$som_freq

#Class assignments with label IDs/XYs, map coordinates, and assignment error
head(summ$som_class)

## ----somplots, echo=TRUE, fig.height=8, fig.width=7.5-------------------------
sommix_bar(som_obj=sommod, varnames=NULL, colormod=NULL,
                     nodelab=TRUE, labtype="IDs", labsize=1,
                     addFreq=TRUE, freqtype="frq", freqsize=0.75,
                     legsymsize=2, leglabsize=1, legtxtlas=2,
                     barstat="MED")

## ----somfitplots, echo=TRUE, fig.height=10, fig.width=7.5---------------------
sommix_eval(sommod, labtype = "IDs", legcex=0.75)

## ----compeval, echo=TRUE, fig.height=10, fig.width=7.5------------------------
compeval(sommod)

