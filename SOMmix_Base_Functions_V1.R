#SOMmix
# A Collection of Functions to support SOM analyses of complex environmental exposures
#Date: 31AUG2020
#Author: John Pearce

############################################################################################
#Install necessary packages
install.packages('kohonen')
install.packages('e1071')
install.packages('fpc')
install.packages('colorspace')
install.packages('spdep')
install.packages('corrplot')


#Load necesary packages
library(kohonen)#self-organizing map algorithim
library(MASS)#Sammons Mapping
library(e1071) #k-nearest neighbor class assignments
library(fpc) #Cluster validation statistics
library(cluster)
library(colorspace) #Color options
library(spdep) #Spatial Statistics
library(corrplot) #Correllogram
############################################################################################

############################################################################################
#Functions to explore data characteristics
############################################################################################
#Explore variables independently
var.eval<-function(data, lab.cex=1, sym.cex=2){
  #Set up data
  data=data
  vars=names(data)
  
  #Set plot specifics
  par(mfrow=c(3,2), mar=c(5,4,2,1), family="serif", xaxs="r", ask=FALSE)
  var.labs=ifelse(nchar(vars)<8,vars,substr(vars, 0,8))
  lab.cex=lab.cex
  
  #Evaluate Variable completeness
  NAs<-apply(apply(data,2,is.na),2,sum)
  NAs.per<-round(100*(NAs/dim(data)[1]),2)
  N<-apply(apply(data,2,complete.cases),2,sum)
  N.per<-round(100*(N/dim(data)[1]),2)
  
  ZERO<-function (x){length(which(x <= 0))}
  ZEROs<-round(apply(data,2,FUN=ZERO),2)
  ZEROs.per<-round(100*(ZEROs/dim(data)[1]),2)
  
  AVG<-round(apply(data,2,FUN=mean, na.rm=TRUE),5)
  MED<-round(apply(data,2,FUN=median, na.rm=TRUE),5)
  SD<-round(apply(data,2,FUN=sd, na.rm=TRUE),5)
  CV<-SD/AVG
  CVrank=rank(-CV, ties.method = "random")
  MN<-round(apply(data,2,FUN=min, na.rm=TRUE),5)
  Q1<-round(apply(data,2,FUN=quantile, na.rm=TRUE)[2,],5)
  Q3<-round(apply(data,2,FUN=quantile, na.rm=TRUE)[4,],5)
  IR=round(apply(data,2,FUN=IQR, na.rm=TRUE),5)
  MX<-round(apply(data,2,FUN=max, na.rm=TRUE),5)
  
  
  summ.tab<-cbind(N, N.per, NAs, NAs.per, ZEROs, ZEROs.per, AVG, MED, SD, CV, CVrank, MN, Q1, MED, Q3, MX, IR)
  
  #Evaluate variable distribution
  SK<-NULL
  KR<-NULL
  
  for (i in 1:dim(data)[2]){
     #Set variable 
     V1<-as.numeric(na.omit(data[,i]))
     
      #Generate Distributional Statistics
      sk<-skewness(V1)
      kr<-kurtosis(V1, type=2)
     
    #Populate Summaries  
    SK<-c(SK,sk) 
    KR<-c(KR,kr) 
    
  }
  
  summ.tab<-data.frame(VARIABLE=vars,summ.tab, Skewness=SK, Kurtosis=KR)
  write.table(x=summ.tab, file="var.eval.summ.tab.csv", sep=",", row.names=FALSE)
  
  #############################################################
  #Evaluation Plots
  cat("*******************************************************************", 
      "\n")
  
  plot(summ.tab$Skewness, pch=20, col="darkblue", ylab="Skewness", 
       xlab="", main="b) Skewness", xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(0), col=c("lightgrey"), lwd=2, lty=2)
  
  plot(summ.tab$Kurtosis, pch=20, col="darkblue", ylab="Kurtosis",
       xlab="", main="b) Kurtosis", xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=0, col="lightgrey", lwd=2, lty=2) #Values below 0.05 not normal
  
  
  plot(summ.tab$SD, pch=20, col="darkblue", ylab="Standard Deviation",
       xlab="", main="c) Standard Deviations", xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(0), col=c("lightgrey"), lwd=2, lty=2)
  
  plot(summ.tab$CV, pch=20, col="darkblue", ylab="mean/sd", 
       xlab="", main="d) Coefficient of Variation (CV)", 
       xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(1), col=c("lightgrey"), lwd=2, lty=2)
  
  plot(summ.tab$NAs.per, pch=20, col="darkblue", ylab="Percentage", 
       xlab="", main="e) Missing Values (%)", ylim=c(0,100), xaxt="n",
       cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(90), col=c("lightgrey"), lwd=2, lty=2)
  
  plot(summ.tab$ZEROs.per, pch=20, col="darkblue", ylab="Percentage", 
       xlab="", main="f) Values less than or equal to zero (%)", ylim=c(0,100), xaxt="n", 
       cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(90), col=c("lightgrey"), lwd=2, lty=2)
  
  print.tab=summ.tab[,c("VARIABLE", "N", "AVG", "SD", "MN", "Q1", "MED","Q3", "MX", "IR",  "CV", "CVrank")]
  print(print.tab)
  cat("*******************************************************************", 
      "\n")
}
############################################################################################

############################################################################################
#Explore pairwise correlations  
cor.eval<-function(data, cor.method="pearson", lab.cex=1){
  
  #Set up data
  data=ex.data
  vars=names(data)
  
  #Set plot specifics

  var.labs=ifelse(nchar(vars)<8,vars,substr(vars, 0,8))
  lab.cex=lab.cex
  
  par(mfrow=c(1,1), mar=c(5,4,2,1), family='serif', cex=lab.cex, ask=FALSE, pty="m")
  
  #Examine group structure via distribution of pairwise distances between objects
  #Identy pairwise correlations
  cor.method="pearson"
  cor.mat<-cor(data, method=cor.method, use="pairwise.complete.obs")
  
  write.table(x=cor.mat, file="cor.summ.tab.csv", sep=",")
  cat("*******************************************************************", 
      "\n")
  cor.summ=round(data.frame(cor.mat),2)
  cor.summ[upper.tri(cor.mat, diag=TRUE)]<-NA
  print(cor.summ)
  cat("*******************************************************************", 
      "\n")
  
  col.r <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(cor.mat, method="color", col=rev(col.r(200)),type="lower",
           order="original", tl.cex=lab.cex,
           tl.col="black", tl.srt=45, 
           diag=FALSE)
}
############################################################################################

############################################################################################
#Explore variance-covariance structure with Principal Component Analysis
pca.eval<-function(data, lab.cex=1){
  data.df=as.data.frame(data)
  vars=names(data.df)
  var.labs=ifelse(nchar(vars)<8,vars,substr(vars, 0,8))
  
  par(mfrow=c(2,2), mar=c(5,4,2,1), ask=FALSE, family='serif')
  pca.mix<-prcomp(data, center=FALSE, scale=FALSE)
  cat("*******************************************************************", 
      "\n")
  print(summary(pca.mix))
  
  pca.loadings<-round(pca.mix$rotation,2)
  
  cat("*******************************************************************", 
      "\n")
  cat("Component Loadings",
      "\n")
  print(pca.loadings)
  cat("*******************************************************************", 
      "\n")
  write.table(x=pca.loadings, file="pca.loadings.summ.csv", sep=",")
  
  #Explore the variance explained by each component
  pr.var=pca.mix$sdev^2
  #Calculate the proportion of variance explained
  pve=round(pr.var/sum(pr.var),2)

  #Plot the component Standard deviations to see which are useful 
  plot(pca.mix$sdev^2, pch=20, ylab="Eigenvalue", xlab="Principal Component",
     cex=lab.cex, type="b", main="a) PCA Scree Plot", col="darkblue")
  abline(h=1, col="darkgrey", lty=2)
  box(which="figure")

  plot(pve, xlab="Principal Component", ylab="Proportion of Variance", ylim=c(0,1), type='b', 
       pch=20, col="darkblue", main="b) Variance Explained", cex=lab.cex)
  points(cumsum(pve), type='b',pch=20, col="darkred")
  abline(h=cumsum(pve)[2], col="darkgrey", lwd=2, lty=2)
  legend("bottomright", pch=20, col=c("darkblue", "darkred"), legend=c("Individual", "Cumulative"))
  box(which="figure")

  PC1<-pca.loadings[,1]
  barplot(PC1, main="c) PC1 Loadings", ylim=c(-1,1), las=2, 
          names.arg=var.labs, cex.names=lab.cex, col="lightgrey")
  box()
  box(which="figure")

  PC2<-pca.loadings[,2]
  barplot(PC2, main="d) PC2 Loadings", ylim=c(-1,1), las=2, 
          names.arg=var.labs, cex.names=lab.cex, col="lightgrey")
  box()
  box(which="figure")
}
############################################################################################

#Measures of Model Fit
somSSS<-function(som.obj){
  x=data.frame(som.obj$data)
  k=as.factor(som.obj$unit.classif)
  
  #Set range of functions to evaluate sum-of-squares 
  #Total Sum of Squares
  TSS <- sum(scale(x, scale = FALSE)^2)
  
  #Within-class sum of squares
  k.n<-length(levels(k))
  WCSS=NULL
  
  for (i in 1:k.n){
    x1=subset(x, k==i)
    wss<-sum(scale(x1, scale = FALSE)^2)
    WCSS<-c(WCSS, wss)        
  }  
  
  #Total sum of within-class sum-of-squares
  Tot.WCSS<-sum(WCSS)
  
  #Between class sum of squares
  BCSS<-TSS-Tot.WCSS
  
  sumofsquares=list("TotalSumSquares"=TSS, "WithinClassSumSquares"=WCSS, "TotWithnClassSumSquares"=Tot.WCSS, "BetweenClassSumSquares"=BCSS)
  return(sumofsquares)
}
somR2<-function(som.obj){
  x=data.frame(som.obj$data)
  k=as.factor(som.obj$unit.classif)
  
  #Set range of functions to evaluate sum-of-squares 
  #Total Sum of Squares
  TSS <- sum(scale(x, scale = FALSE)^2)
  
  #Within-class sum of squares
  k.n<-length(levels(k))
  WCSS=NULL
  
  for (i in 1:k.n){
    x1=subset(x, k==i)
    wss<-sum(scale(x1, scale = FALSE)^2)
    WCSS<-c(WCSS, wss)        
  }  
  
  #Total sum of within-class sum-of-squares
  Tot.WCSS<-sum(WCSS)
  
  #Between class sum of squares
  BCSS<-TSS-Tot.WCSS
  
  #Within-Between Ratio 
  Pseudo.F<-BCSS / Tot.WCSS
  
  
  rsq=BCSS/TSS
  
  return(rsq)
}
somFstat<-function(som.obj){
  x=data.frame(som.obj$data)
  k=as.factor(som.obj$unit.classif)
  
  #Set range of functions to evaluate sum-of-squares 
  #Total Sum of Squares
  TSS <- sum(scale(x, scale = FALSE)^2)
  
  #Within-class sum of squares
  k.n<-length(levels(k))
  WCSS=NULL
  
  for (i in 1:k.n){
    x1=subset(x, k==i)
    wss<-sum(scale(x1, scale = FALSE)^2)
    WCSS<-c(WCSS, wss)        
  }  
  
  #Total sum of within-class sum-of-squares
  Tot.WCSS<-sum(WCSS)
  
  #Between class sum of squares
  BCSS<-TSS-Tot.WCSS
  
  #Within-Between Ratio 
  Pseudo.F<-BCSS / Tot.WCSS
  
  return(Pseudo.F)
}
somAIC = function(som.obj){
  
  m = ncol(data.frame(som.obj$codes))
  n = length(som.obj$unit.classif)
  p = nrow(data.frame(som.obj$codes))
  
  x=data.frame(som.obj$data)
  k=as.factor(som.obj$unit.classif)
  
  #Set range of functions to evaluate sum-of-squares 
  #Total Sum of Squares
  TSS <- sum(scale(x, scale = FALSE)^2)
  
  #Within-class sum of squares
  k.n<-length(levels(k))
  WCSS=NULL
  
  for (i in 1:k.n){
    x1=subset(x, k==i)
    wss<-sum(scale(x1, scale = FALSE)^2)
    WCSS<-c(WCSS, wss)        
  }  
  
  #Total sum of within-class sum-of-squares
  Tot.WCSS<-sum(WCSS)
  
  #Between class sum of squares
  BCSS<-TSS-Tot.WCSS
  
  D = Tot.WCSS
  
  AIC = D + 2*m*p
  
  return(AIC)
}  

############################################################################################
#Explore grouping structure (i.e., clustering) strategies applied in cluster anaylsis 
grp.eval<-function(data, kmn, kmx, iter.max=100){
  #Set data
  data=data.matrix(data)
  #set Data distance matrix
  dist.mat<-data.matrix(dist(x=data, method="euclidean"))
  
  AW=NULL
  AB=NULL
  CH=NULL
  ASW=NULL
  
  
 for (i in kmn:kmx){
    set.seed(1000)
    ex.som<-som(X=data, grid=somgrid(xdim=i, ydim=1, topo=c("rectangular")), 
                rlen=iter.max, mode="online", dist.fcts="euclidean")
    #SOM-based classifications 
    som.class=ex.som$unit.classif
    
    
    CLUSTER.STATS=cluster.stats(d=dist.mat, clustering=as.integer(som.class))
    
    #Cluster Cohesion
    aw=CLUSTER.STATS$average.within
    
    #Cluster Seperation
    ab=CLUSTER.STATS$average.between
    
    #Calinski-Harabasz Index
    ch=CLUSTER.STATS$ch
    
    #Silhoutte Width
    sw.summ=CLUSTER.STATS$clus.avg.silwidths
    sw.summ
    asw=CLUSTER.STATS$avg.silwidth
    
    AW=c(AW, aw)
    AB=c(AB,ab)
    CH=c(CH,ch)
    ASW=c(ASW,asw)
    
  }
  
  #Examine group structure using kmeans and common cluster statistics
  CLUST.EVAL<-data.frame(K=seq(kmn,kmx,1), Cohesion=AW, Seperation=AB, CalH=CH, SilW=ASW)
  
  par(mfrow=c(3,2), mar=c(5,4,2,1), cex=1, family='serif')
  
  #Examine group structure via distribution of pairwise distances between objects
  hist(dist.mat, main="a) Pairwise Object Distances", xlab="Euclidean Distance", 
       cex=2, col="lightgrey")
  box()
  box(which="outer")
  
  #Apply multidimensional scaling to visualize clustering structure
  loc <- cmdscale(dist.mat)
  x <- loc[, 1]
  y <- -loc[, 2] # reflect so North is at the top
  ## note asp = 1, to ensure Euclidean distances are represented correctly
  plot(x, y, type = "p", xlab = "MDS X", ylab = "MDS Y", asp = 1, axes = TRUE,
       main = "b) Multi-dimensional Scaling", pch=20, col="darkblue")
  #text(x, y, rownames(loc), cex = 0.6)
  box()
  box(which="outer")
  
  ##########################################################################
    #Create Plots
  plot(x=CLUST.EVAL[,"K"], y=CLUST.EVAL[,"Cohesion"], pch=20, col="darkblue", xlab="Number of profiles (k)",
       ylab="Distance", main="c) Average Within-Distance ",type="b")
  
  plot(x=CLUST.EVAL[,"K"], y=CLUST.EVAL[,"Seperation"], pch=20, col="darkblue", xlab="Number of profiles (k)",
       ylab="Distance", main="d) Average Between-Distance ",type="b")
  
  plot(x=CLUST.EVAL[,"K"], y=CLUST.EVAL[,"CalH"], pch=20, col="darkblue", xlab="Number of profiles (k)",
       ylab="Pseudo F-Statistic", main="e) Calinski Harabaz ",type="b")
  
  plot(x=CLUST.EVAL[,"K"], y=CLUST.EVAL[,"SilW"], pch=20, col="darkblue", xlab="Number of profiles (k)",
       ylab="Avg Width", main="f) Avg Silhouette Width ",type="b", ylim=c(-1,1))
  
  
  cat("*******************************************************************", 
      "\n")
  print(CLUST.EVAL)
  cat("*******************************************************************", 
      "\n")
}
############################################################################################

############################################################################################
############################################################################################
#Functions to evaluate application of SOM
############################################################################################
#Find optimal seed value 
som.seed<-function(data, nstart=10, iter.max=100, somx, somy, grid.topo="rectangular"){
  #Set random initializations
  set.seed(1)
  ran.inits=sample(1:100,nstart)
  ran.inits
  #Set eval object
  init.eval<-NULL
  
  for(i in 1:nstart){
    set.seed(ran.inits[i])
    ex.som<-som(X=data.matrix(data), grid=somgrid(somx,somy,grid.topo), 
                rlen=iter.max, mode="online", dist.fcts="euclidean")
    QE=mean(ex.som$distances)
    init.eval=c(init.eval,QE)
  }
  
  seed.value=ran.inits[which(init.eval == min(init.eval))]
  return(seed.value)
}
############################################################################################

#Summary Function to for summarizing and assessing SOM via kohonen package
som.fit.summ<-function(som.obj){
  #Print Summary
  print(summary(som.obj))
  cat("*******************************************************************", 
      "\n")
  cat("Overall SOM Model Fit Statistics",
      "\n")
  #somsss=somSSS(som.obj)
  somr2=somR2(som.obj)
  somaic=somAIC(som.obj)
  somfstat=somFstat(som.obj)
  QE<-mean(som.obj$distances)
  RMSE<-sqrt(mean(som.obj$distances^2))
  
  mod.stats=data.frame(QE, "R2"=somr2, RMSE)
  print(mod.stats)
  cat("*******************************************************************", 
      "\n")
  
  
  #Extract Training Data
  data=data.frame(som.obj$data)
  var.names<-names(data)
  
  #Extract SOM dimensions and create labels
  somx=som.obj$grid$xdim
  somy=som.obj$grid$ydim
  som.k=somx*somy
  XYs=paste(round(som.obj$grid$pts[,1],2),round(som.obj$grid$pts[,2],2), sep="")
  SOM_IDs=1:som.k
  #SOM Coordinates and reference table
  som.coords<-data.frame(X=som.obj$grid$pts[,1],Y=som.obj$grid$pts[,2])
  som.class.xy<-data.frame(SOM_ID=factor(SOM_IDs), XY=factor(XYs),X=as.numeric(som.obj$grid$pts[,1]),
                           Y=as.numeric(som.obj$grid$pts[,2]))
  
  #SOM Profiles
  profiles<-data.frame(som.obj$codes)
  row.names(profiles)<-SOM_IDs
  
  #print("Map Profiles and Coordinate Linkage Complete")
  
  #Create data table of SOM predictions 
  #Classifications
  som.class<-som.obj$unit.classif
  class.tab<-data.frame(OBS=1:length(som.class), SOM_ID=som.class)
  
  class.tab2<-merge(class.tab, som.class.xy, by="SOM_ID")
  class.tab2<-class.tab2[,c("OBS", "SOM_ID","XY","X","Y")]
  class.tab3<-class.tab2[do.call(order,class.tab2),]
  class.tab3$ERROR<-som.obj$distances
  rownames(class.tab3)<-class.tab3$OBS
  
  profiles2<-data.frame(SOM_ID=rownames(profiles), profiles)
  som.preds<-merge(class.tab3, profiles2, by="SOM_ID",all.x=TRUE)
  
  som.variable<-som.preds[order(som.preds$OBS),]
  
  #print('SOM Predictions Table Complete')
  
  #Frequency Summary
  som.n<-table(class.tab$SOM_ID)
  som.freq<-round(100*(table(class.tab$SOM_ID)/sum(table(class.tab$SOM_ID))),2)
  
  freq.tab<-data.frame(som.n, som.freq)
  freq.tab2<-merge(som.class.xy, freq.tab, by.x="SOM_ID", by.y="Var1", all.x=TRUE)
  freq.tab2$Var1.1<-NULL
  colnames(freq.tab2)<-c("SOM_ID","XY","X","Y","N","FREQ")
  freq.tab2<-freq.tab2[order(as.numeric(freq.tab2$SOM_ID)),]
  
  cat("*******************************************************************", 
      "\n")
  cat("SOM Profile Coordinates and Frequencies",
      "\n")
  print(freq.tab2)
  cat("*******************************************************************", 
      "\n")
  
  ####################################################################
  #Internal evaluation of SOM Fit using distance measures
  QE<-mean(som.obj$distances)
  error<-aggregate(x=som.obj$distances, by=list(som.obj$unit.classif), FUN=mean, na.omit=TRUE)
  error.tab<-data.frame(error) 
  colnames(error.tab)<-c("SOM_ID", "WITHIN")
  error.tab2<-merge(som.class.xy,error.tab, by.x="SOM_ID", by.y="SOM_ID", all.x=TRUE)
  profile.dist<-data.matrix(dist(data.frame(som.obj$codes)))
  profile.dist2<-apply(profile.dist, MARGIN=1, FUN=mean, na.rm=TRUE)
  error.tab2$BETWEEN<-profile.dist2
  error.tab2$WB.RATIO=error.tab2$WITHIN/error.tab2$BETWEEN
  error.tab3<-error.tab2[order(as.numeric(error.tab2$SOM_ID)),]
  
  #Extract cluster' evaluation statistics for each SOM profile 
  dist.summ<-error.tab3
  #print("Profile Distance Evaluation Complete")
  cat("*******************************************************************", 
      "\n")
  cat("Internal Class Distance Measures for SOM Profiles ",
      "\n")
  print(dist.summ)
  cat("*******************************************************************", 
     "\n")
  ####################################################################
  #Construct Evaluation Data Sets
  data.eval=data.frame(class.tab, data)
  data.eval.xy<-merge(data.eval, som.class.xy, by="SOM_ID")
  data.eval.xy.preds<-merge(data.eval.xy, som.preds, by="OBS",
                            suffixes=c("",".p"))
  ####################################################################
  ########################################################################
  #Variable statistics for SOM-based predictions
  ## Evaluate discriminatory power of SOM-based profiles scores
  R.p=NULL
  RMSE.p=NULL
  MAE.p=NULL
  
  #Asses correlation and precision measures
  for (j in 1:length(var.names)){
    var.n=var.names[j]
    D1<-data.eval.xy.preds
    Y=D1[,var.n]
    X=D1[,paste(var.n, ".p", sep="")]
    
    #Identify Correlation
    r.p<-cor(Y,X, method="pearson")
    
    #Extract RMSE
    err.p<-Y-X
    rmse.p<-round(sqrt(mean(err.p^2)),4)
    mae.p<-round(mean(abs(err.p)),4)
    RMSE.p<-c(RMSE.p, rmse.p)
    MAE.p<-c(MAE.p, mae.p)
    R.p<-c(R.p,round(r.p,2))
  }
  
  var.eval<-data.frame(VARS=var.names,R.P=as.numeric(R.p),
                       RMSE.P=RMSE.p, MAE.P=MAE.p)
  
  var.eval$R_P_Rank<-rank(-var.eval$R.P)
  var.eval$RMSE_P_Rank<-rank(var.eval$RMSE.P)
  var.eval$MAE_P_Rank<-rank(var.eval$MAE.P)
  
  colnames(var.eval)<-c("Variable", "R","RMSE", "MAE", 
                        "R_Rank","RMSE_Rank","MAE_Rank")
  
  #print("Profile-Variable Assessment Complete")
  cat("*******************************************************************", 
      "\n")
  cat("Variable-level Measures of SOM Model Fit",
      "\n")
  print(var.eval)
  cat("*******************************************************************", 
      "\n")
  #######################################################################
  ######################################################################
  #Evaluate 'spatial'organization of map via global tests of spatial autocorrelation
  if (somy>1){
    #Convert grid to spatial coordinates
    som.xy<-SpatialPoints(som.coords)
    
    #Global Network
    som.knn<-knn2nb(knearneigh(som.xy, k=round(som.k*.2, 0)))
    #Neighborhood Distances
    som.dist<-unlist(nbdists(som.knn, coords=som.xy))
    
    MI=NULL
    MI_P=NULL
    GI=NULL
    GI_P=NULL
    
    for (j in 1:length(var.names)){
      if (sd(profiles[,j]) > 0) {
      mt<-moran.test(x=profiles[,j], nb2listw(som.knn), na.action=na.exclude)
      mi<-as.numeric(mt$statistic)
      mi.p<-as.numeric(mt$p.value)
      
      gt<-geary.test(x=profiles[,j], nb2listw(som.knn))
      gi<-as.numeric(gt$statistic)
      gi.p<-as.numeric(gt$p.value)
      } else {
        mi=NA
        mi.p=NA
        gi=NA
        gi.p=NA
      }
      
      MI<-c(MI, mi)
      MI_P<-c(MI_P, mi.p)
      GI<-c(GI,gi)
      GI_P<-c(GI_P,gi.p)
    }
    
    sp.cor.tab<-data.frame(var.names,MI, MI_P, GI, GI_P)
    colnames(sp.cor.tab)<-c("Variable", "Morans I", "MI_Pval", "Gearys C", "GC_Pval")
    sp.cor.tab$MI_Rank<-rank(-sp.cor.tab$`Morans I`)
    sp.cor.tab$GI_Rank<-rank(-sp.cor.tab$`Gearys C`)
    
    #print("Map Structure/Autocorrelation Evaluation Complete")
  } else { sp.cor.tab<-NA}
    #print("1-d SOM NO Map Structure/Autocorrelation Evaluation")}
  cat("*******************************************************************", 
      "\n")
  cat("Variable-level Measures of SOM Organization (i.e., spatial structure)",
      "\n")
  print(sp.cor.tab)
  cat("*******************************************************************", 
      "\n")
  ###########################################################################
    #Pull Training Data
    D2<-data
    
    #Identify SOM Class Assignments
    GRP<-som.obj$unit.classif
    N<-table(GRP)
    FREQ<-round(100*(N/sum(N)),1)
    
    means<-aggregate(x=D2[,var.names], by=list(GRP), FUN=mean, na.rm=TRUE)
    means<-round(means[,-1],2)
    
    med<-aggregate(x=D2[,var.names], by=list(GRP), FUN=median, na.rm=TRUE)
    med<-round(med[,-1],2)
    
    stdev<-aggregate(x=D2[,var.names], by=list(GRP), FUN=sd, na.rm=TRUE)
    stdev<-round(stdev[,-1],2)
    
    stderr <- function(x, na.rm=FALSE) {
      if (na.rm) x <- na.omit(x)
      sqrt(var(x)/length(x))
    }
    
    std.err<-aggregate(x=D2[,var.names], by=list(GRP), FUN=stderr, na.rm=TRUE)
    std.err<-round(std.err[,-1],2)
    
    cv<-function (x) {sd(x)/mean(x)}
    
    CVar<-aggregate(x=D2[,var.names], by=list(GRP), FUN=cv)
    CVar<-round(CVar[,-1],2)
    
    mnv<-aggregate(x=D2[,var.names], by=list(GRP), FUN=min)
    mnv<-round(mnv[,-1],2)
    
    mxv<-aggregate(x=D2[,var.names], by=list(GRP), FUN=max)
    mxv<-round(mxv[,-1],2)
    
    wgv<-aggregate(x=D2[,var.names], by=list(GRP), FUN=var)
    wgv<-round(wgv[,-1],2)
    
    var.summ.stats<-list(Mean=means,Medians=med, St_Dev=stdev,St_Err=std.err,Coef_Var=CVar, Min=mnv, Max=mxv,
                     Var=wgv)
  #####################################################################################
  #####################################################################################
  summ.list<-list(MODEL_STATS=mod.stats, COORDINATES=som.class.xy, FREQUENCIES=freq.tab2, PROFILES=profiles, PREDICTIONS=som.variable,
                  PROFILE_DISTANCES=dist.summ, COMPONENT_EVALUATION=var.eval, MAP_ORGANIZATION=sp.cor.tab,
                  COMPONENT_STATISTICS=var.summ.stats)
}
############################################################################################

############################################################################################
############################################################################################
#Develop a set of functions to evaluate a range of common Self-Organizing Map sizes
map.size.eval<-function(data, kmx, nstart=10, iter.max=100, grid.topo="rectangular"){
  
  #Set common dimensions
  som.x=c(2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20)
  som.y=c(2,2,3,3,4,4,5,5,6,6,7,7,8,8,9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20)
  som.xy=som.x*som.y
  
  max.size<-which(abs(som.xy-kmx)==min(abs(som.xy-kmx)))
  
  #set Data
  data=data
  var.names=names(data.frame(data))
  n.vars=length(var.names)
  n.obs=dim(data)[1]
  
  MOD.TAB<-data.frame()
  EVAL.TAB<-data.frame()
  FREQ.TAB<-data.frame()
  MAP.TAB<-data.frame()
  
  
  for (i in 1:max.size){
    somx=som.x[i]
    somy=som.y[i]
    
    data=data.matrix(data)
    
    #Identify optimal seed value
    opt.seed<-som.seed(data=data, nstart=nstart, iter.max=iter.max, somx=somx, somy=somy, grid.topo=grid.topo)
    
    #Apply a SOM
    set.seed(opt.seed)
    som.mod<-som(X=data, grid=somgrid(somx,somy,grid.topo),
                  rlen=iter.max, alpha=c(0.05, 0.01))
    som.summ=som.fit.summ(som.mod)
    
    mod.tab=data.frame(SOMX=somx, SOMY=somy, k=somx*somy, 
                       R2=as.numeric(som.summ$MODEL_STATS$R2),
                       QE=as.numeric(som.summ$MODEL_STATS$QE),
                       RMSE=as.numeric(som.summ$MODEL_STATS$RMSE))
       
    eval.tab<-data.frame(VAR.NAM=var.names,SOMX=somx, SOMY=somy, k=somx*somy,
                       R.P=som.summ$COMPONENT_EVALUATION$R, 
                       RMSE.P=som.summ$COMPONENT_EVALUATION$RMSE,
                       MAE.P=som.summ$COMPONENT_EVALUATION$MAE)
    
    freq.tab<-data.frame(SOMX=somx, SOMY=somy,k=somx*somy, 
                         N=as.numeric(som.summ$FREQUENCIES$N), 
                         FREQ=as.numeric(som.summ$FREQUENCIES$FREQ))
    
    map.tab<-data.frame(VAR.NAM=var.names,SOMX=somx, SOMY=somy,k=somx*somy, 
                        MI=som.summ$MAP_ORGANIZATION$`Morans I`, 
                        MI_P=som.summ$MAP_ORGANIZATION$MI_Pval,
                        GC=som.summ$MAP_ORGANIZATION$`Gearys C`, 
                        GC_Pval=som.summ$MAP_ORGANIZATION$GC_Pval)
    
    MOD.TAB<-rbind(MOD.TAB,mod.tab)
    EVAL.TAB<-rbind(EVAL.TAB, eval.tab)
    FREQ.TAB<-rbind(FREQ.TAB,freq.tab) 
    MAP.TAB<-rbind(MAP.TAB, map.tab)
    
  }
  
  #Plots
  par(mfrow=c(3,2), mar=c(4,4,1.5,1), family='serif', ask=FALSE)
  
  plot(MOD.TAB$RMSE~MOD.TAB$k, ylab="Distance", xlab="Number of profiles (k)", 
       main="a) Overall Root-Mean-Square-Error (RMSE)", type="b", pch=19, cex=2, col="darkgrey", 
       ylim=c(0,max(MOD.TAB$RMSE)*1.2))
  
  plot(MOD.TAB$R2~MOD.TAB$k, ylab="BCSS/TSS", xlab="Number of profiles (k)", 
       main="b) Overall R-square", type="b", pch=19, cex=2, col="darkgrey",
       ylim=c(0,1))
  
  boxplot(EVAL.TAB$R.P~EVAL.TAB$k, ylab="R", xlab="Number of profiles (k)",
          main="c) Variablewise Correlations", ylim=c(-1,1))
  
  boxplot(EVAL.TAB$RMSE.P~EVAL.TAB$k, ylab="RMSE", xlab="Number of profiles (k)",
          main="d) Variablewise RMSE", ylim=c(0,max(EVAL.TAB$RMSE)*1.2))
  
  boxplot(FREQ.TAB$N~FREQ.TAB$k, ylab="N", xlab="Number of profiles (k)",
          main="e) Profile Frequency Distributions")
  
   boxplot(MAP.TAB$MI~MAP.TAB$k, ylab="Moran's I", xlab="Number of profiles (k)",
          main="f) Variablewise Map Spatial Autocorrelation")
  
  cat("--------------------------------------------------------------------------------",
      "\n")
  MAP.EVAL<-list(MOD.FIT=MOD.TAB, VAR.FIT=EVAL.TAB, MAP.FREQ=FREQ.TAB, MAP.STRUCTURE=MAP.TAB)
  print(MAP.EVAL)
  cat("--------------------------------------------------------------------------------",
      "\n")
}
############################################################################################

############################################################################################
#Functions for plotting SOM Results
############################################################################################
############################################################################################
#Generate Stars Plot
plot.somgrid <- function(x, xlim, ylim, ...){
  ## The following two lines leave equal amounts of space on both
  ## sides of the plot if no xlim or ylim are given
  if (missing(xlim)) xlim <- c(0, max(x$pts[,1]) + min(x$pts[,1]))
  if (missing(ylim)) ylim <-  c(max(x$pts[,2]) + min(x$pts[,2]), 0)
  eqscplot(xlim, ylim, axes = FALSE,
           type = "n", xlab = "", ylab = "", ...)
}
#Profile plots with frequency information
map.profile.plot<-function(som.obj, lab.cex=1, legend.loc="bottom", 
                           leg.sym.cex=1, leg.lab.cex=1, 
                           col.mod=diverge_hcl(p, h=c(130,43, c=100, I=c(70,90)))){
  #Extract info for plotting
  somx<-som.obj$grid$xdim
  somy<-som.obj$grid$ydim
  profiles<-as.data.frame(som.obj$codes)
  som.summ<-som.fit.summ(som.obj)
  som.topo=som.obj$grid$topo
  
  #Set labels
  IDs=1:dim(profiles)[1]
  p=dim(profiles)[2]
  var.labs<-dimnames(profiles)[[2]]
  
  #Set Color Scheme
  color.mod=col.mod
  
  #Create Profile plot
  par(mfrow=c(1,1),mar=c(.6,.6,.6,.6), family="serif", pty="m")
  plot.somgrid(som.obj$grid)
  stars(x=som.summ$PROFILES, locations=coordinates(som.summ$COORDINATES[,3:4]), 
        draw.segments=TRUE, axes=FALSE, scale=TRUE, 
        len = 0.4, col.segments=color.mod, 
        labels=NULL,ylim=c(0.5,somy), xlim=c(0.5,somx+.5), add=TRUE)
  
  symbols(coordinates(som.summ$COORDINATES[,3:4]), 
          squares = rep(1, nrow(som.summ$COORDINATES)),
          inches = FALSE, add = TRUE,
          fg = "black", bg = NA)
  
  axis(1, at=som.obj$grid$pts[,1], pos=0.5)
  text(x=mean(som.summ$COORDINATES$X), y=0.2, labels="SOM X") 
  axis(2, at=som.obj$grid$pts[,2], pos=0.5)
  text(x=0.25, y=mean(som.summ$COORDINATES$Y), labels="SOM Y", srt=90) 
  
  text(x=som.summ$COORDINATES$X, y=som.summ$COORDINATES$Y+.4,
       labels=paste("[",som.summ$COORDINATES$X,",", som.summ$COORDINATES$Y ,"]", sep=""),
       font=2, cex=lab.cex+.5)
  
  min.x<-min(som.summ$COORDINATES[,3])
  min.y<-min(som.summ$COORDINATES[,4])
  
  legend(legend.loc, legend=var.labs, pch=15, 
         col=col.mod,
         pt.cex=leg.sym.cex, cex=leg.lab.cex, angle=45, ncol=round(length(var.labs)/2),
         inset=0)
}

############################################################################################
############################################################################################
#Bar Plot with error bars of Profiles
#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.01,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, col="darkgrey", ...)
}

#Plot Function
map.bar.plot<-function(som.obj, lab.cex=1, label.loc="bottom", legend.cex=1, 
                       col.mod=diverge_hcl(p, h=c(130,43, c=100, I=c(70,90)))){
  #Get summary output
  som.summ<-som.fit.summ(som.obj)
  #Get SOM Dimensions
  somx<-som.obj$grid$xdim
  somy<-som.obj$grid$ydim
  #Get SOM Codes
  profiles<-as.data.frame(som.obj$codes)
  profiles
  
  #Get variable names
  var.names<-names(profiles)
  n.col=length(var.names)
  p=n.col
  
  #Set Bar Color Model
  color.mod<-col.mod
  
  #Extract data with Class Assignments
  FREQ<-som.summ$FREQUENCIES$FREQ
  N<-som.summ$FREQUENCIES$N
  #Generate a layout matrix
  par(mar=c(.6,3,2,.6), pty="m")
  M<-matrix(data=1:(somx*somy), 
            nrow=somy, ncol=somx, byrow=TRUE)
  M1=apply(t(M),1,rev)
  
  M2<-rbind(M1, rep(1+(somx*somy),somx))
  
  mf <- layout(M2, widths = rep(1, ncol(M2)),
               heights = c(rep(1, nrow(M2)-1), 0.5), respect = FALSE) 
  layout.show(mf)
  
  #Generate barplots
  plot.data<-som.summ$COMPONENT_STATISTICS$Mean
  plot.error<-som.summ$COMPONENT_STATISTICS$St_Dev
  plot.labs<-paste(som.summ$COORDINATES$X, ",",som.summ$COORDINATES$Y, sep="") 
  
  for(i in 1:(somx*somy)){
    plot.vals<-as.numeric(plot.data[i,])
    plot.err<-as.numeric(plot.error[i,])
    
    bp<-barplot(plot.vals, col=color.mod, axes=TRUE, axisnames=FALSE,
            ylim=c(1.2*min(plot.data), 1.2*max(plot.data)))
    error.bar(bp,plot.vals,plot.err)
    abline(h=0, lty=2, col="black")
    title(main=paste("[", plot.labs[i] ,"]", sep=""), font=2, cex.main=lab.cex+.5)
    legend(legend=c(paste(FREQ[i], "%",sep=""), paste("n=",N[i],sep="")), cex=lab.cex,text.font=2,
           x=label.loc, bty="n", horiz=FALSE)
    box()
    box(which="outer")
    
  }
  
  box(which="outer")
  par(mar=c(.1,.1,.1,.1))
  plot(x=1:n.col, y=rep(1.25,n.col), pch=15, cex=legend.cex+1, 
       col=color.mod, 
       xlab=var.names, ylim=c(1,1), axes=FALSE)
  box()
  box(which="outer")
  text(x=1:n.col, y=rep(1,n.col), var.names, pos=1, cex=legend.cex, srt=90)
  
}
############################################################################################
############################################################################################
#Generate Component Plot
map.comp.plot<-function(som.obj, col.mod=diverging_hcl, border.style="straight"){
  #Extract info for plotting
  somx<-som.obj$grid$xdim
  somy<-som.obj$grid$ydim
  profiles<-as.data.frame(som.obj$codes)
  som.summ<-som.fit.summ(som.obj)
  
  #Set labels
  IDs=1:dim(profiles)[1]
  p=dim(profiles)[2]
  var.labs<-dimnames(profiles)[[2]]
  color.mod=col.mod
  
  #Standardized component plots
  par(mfrow=c(sqrt(p)+1,sqrt(p)), mar=c(1,4,2,1), family='serif')
  for(i in 1:length(var.labs)){
    plot(som.obj, type="property", property=profiles[,i], 
         main=var.labs[i], palette.name=color.mod, shape=border.style)
    axis(1, at=som.obj$grid$pts[,1], pos=0.5)
    axis(2, at=som.obj$grid$pts[,2], pos=0.5)
  }
  
  plot(som.obj, type="property", property=rowSums(getCodes(som.obj,1)), 
       main="Cumulative Sum", palette.name=color.mod, shape=border.style)
  axis(1, at=som.obj$grid$pts[,1], pos=0.5)
  axis(2, at=som.obj$grid$pts[,2], pos=0.5)
 
}
############################################################################################
############################################################################################
#Map Fit Evaluation Functions
map.fit.plots<-function(som.obj, lab.cex=1){
  #Extract data
  data=data.frame(som.obj$data)
  str(data)
  var.names<-names(data)
  n.obs=dim(data)[1]
  n.vars=dim(data)[2]
  
  #set Data distance matrix
  dist.mat<-dist(x=data, method="euclidean")
  
  #Extract summary
  som.summ<-som.fit.summ(som.obj)
  
  #Get SOM Dimensions
  somx<-som.obj$grid$xdim
  somy<-som.obj$grid$ydim
  k=somx*somy
  
  #Get SOM Profiles
  som.profiles<-data.frame(som.obj$codes)
  som.obj$grid$pts
  row.names(som.profiles)<-paste(som.obj$grid$pts[,1], som.obj$grid$pts[,2], sep="")
  som.profiles
  som.class<-som.obj$unit.classif
  
  set.seed(100)
  k.mod<-kmeans(x=data, centers=k, nstart=10, algorithm="Hartigan-Wong", 
                iter.max=length(som.obj$changes))
  k.profiles=data.matrix(k.mod$centers)
  k.profiles
  k.class=k.mod$cluster
  
  #Compare class assigments for k-means and SOM
  class.comp=cluster.stats(d=dist.mat, clustering=som.class, alt.clustering = k.class, compareonly=TRUE)
  
  #Merge data for comparision
  data2=rbind(data, som.profiles,k.profiles) 
  #set Data distance matrix
  dist.mat2<-dist(x=data2, method="euclidean")
  
  
  par(mfrow=c(2,3), mar=c(5,4,2,1),ask=FALSE, family='serif', pty="m")
  
  #Apply Ward's Clustering to SOM map
  h.som<-hclust(dist(som.profiles, method="euclidean"), method="ward.D2")
  plot(h.som, xlab="", main="Ward's Clustering of SOM Profiles", sub="")
  box()
  
  #Generate a Sammon's Map to look at Map Distortion
  mapdist<-sammon(d=dist(data.matrix(som.summ$PROFILES)), y=data.matrix(som.summ$COORDINATES[,c("X","Y")]),
                  k=2, niter=1000)
  plot(mapdist$points, pch=20, col="lightgrey", cex=5, xlab="Sammons X", ylab="Sammons Y",
       main="b) Distance preserving projection of SOM Profiles")
  text(mapdist$points, labels = as.character(1:nrow(som.summ$PROFILES)))
  
  
  #Apply multidimensional scaling to visualize clustering structure
  loc <- cmdscale(dist.mat2)
  mds.x <- loc[, 1]
  mds.y <- -loc[, 2] # reflect so North is at the top
  ## note asp = 1, to ensure Euclidean distances are represented correctly
  plot(mds.x, mds.y, type = "p", xlab = "MDS X", ylab = "MDS Y", asp = 1, axes = TRUE,
       main = "c) SOM Profile Comparison", pch=c(rep(20, n.obs), rep(1,k), rep(3,k)), 
       col=c(rep("lightgrey",n.obs), rep("black", k), rep("black", k)))
  legend("bottomright", legend=c("OBS", "SOM", "K-Means"), pch=c(20,1,3), 
         col=c("lightgrey","black","black"), cex=1)
  
  box()
  box(which="outer")
  
  
  ## SOM-Variable Evaluation
  mapvareval<-som.summ$COMPONENT_EVALUATION
  barplot(height=mapvareval$R, names.arg=mapvareval$Variable, 
          las=2, cex.names=lab.cex, main="d) Variablewise Correlations", ylab="Correlation (R)", ylim=c(-1,1), col="lightgrey")
  abline(h=mean(mapvareval$R, na.rm=TRUE), col="darkgrey", lty=2)
  box()
  
  barplot(height=mapvareval$RMSE, names.arg=mapvareval$Variable, 
          las=2, cex.names=lab.cex, main="e) Variablewise RMSE", ylab="Distance", ylim=c(0,(max(mapvareval$RMSE)*1.2)), col="lightgrey")
  abline(h=mean(mapvareval$RMSE, na.rm=T), col="darkgrey", lty=2)
  box()
  
  ##Map Organization
  maporganization<-som.summ$MAP_ORGANIZATION
  barplot(height=maporganization$'Morans I', names.arg=maporganization$Variable, 
          las=2, cex.names=lab.cex, main="f) Variablewise Spatial Structure on SOM Grid", ylab="Moran's I",
          col="lightgrey")
  abline(h=mean(maporganization$'Morans I', na.rm=T), col="darkgrey", lty=2)
  box()
  
  print(class.comp)
  
}  
############################################################################################

############################################################################################
#Extract exposure variable from SOM analyses
############################################################################################
exp.metric<-function(som.obj){
  #Extract Data
  som.data<-data.frame(som.obj$data)
  OBS_IDs<-row.names(som.data)
  
  #Extract class assignments
  unit.classif<-som.obj$unit.classif
  class.tab<-data.frame(OBS_IDs=OBS_IDs,SOM_CLASS=unit.classif)
  
  #Coordinates
  XYs<-data.frame(SOM_CLASS=1:length(som.obj$grid$pts[,1]), 
                  SOMX=som.obj$grid$pts[,1],
                  SOMY=som.obj$grid$pts[,2])
  
  exp.tab<-merge(class.tab, XYs, by="SOM_CLASS", all.x=TRUE)
  exp.tab$SOM_ERROR<-som.obj$distances
  return(exp.tab)
}







