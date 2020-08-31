#SOMmix
# A Collection of Functions to support SOM analyses of complex environmental exposures
#Date: 8JULY2020
#Author: John Pearce

############################################################################################
#Install necessary packages
install.packages('kohonen')
install.packages('e1071')
install.packages('fpc')
install.packages('colorspace')
install.packages('spdep')
install.packages('corrplot')
install.packages('NBClust')

#Load necesary packages
library(kohonen)#self-organizing map algorithim
library(MASS)#Sammons Mapping
library(e1071) #k-nearest neighbor class assignments
library(fpc) #Cluster validation statistics
library(colorspace) #Color options
library(spdep) #Spatial Statistics
library(corrplot) #Correllogram
library(NbClust)
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
  cat("The var.eval() function is designed to enhance understanding of the distribution of each variable in your data, the relative variability, and the proportion of missing and zero values.",
      "\n")
  
  plot(summ.tab$Skewness, pch=20, col="darkblue", ylab="Skewness", 
       xlab="", main="b) Skewness", xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(0), col=c("lightgrey"), lwd=2, lty=2)
  cat("*******************************************************************", 
      "\n")
  cat("a) Skewness is asymmetry in a statistical distribution. Skewness for the normal distribution equals 0, with positive values indicating data are skewed right and negative values indicating data are skewed left.The rule of thumb seems to be: If the skewness is between -0.5 and 0.5, the data are fairly symmetrical; If the skewness is between -1 and - 0.5 or between 0.5 and 1, the data are moderately skewed; If the skewness is less than -1 or greater than 1, the data are highly skewed",
      "\n")
  
  plot(summ.tab$Kurtosis, pch=20, col="darkblue", ylab="Kurtosis",
       xlab="", main="b) Kurtosis", xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=0, col="lightgrey", lwd=2, lty=2) #Values below 0.05 not normal
  
  cat("*******************************************************************", 
      "\n")
  cat("b) The Kurtosis statistic reflects tail-heaviness of distribution. A distribution with a positive kurtosis value indicates that the distribution has heavier tails and a sharper peak than the normal distribution. A distribution with a negative kurtosis value indicates that the distribution has lighter tails and a flatter peak than the normal distribution.",
      "\n")
  
  plot(summ.tab$SD, pch=20, col="darkblue", ylab="Standard Deviation",
       xlab="", main="c) Standard Deviations", xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(0), col=c("lightgrey"), lwd=2, lty=2)
  
  cat("*******************************************************************", 
      "\n")
  cat("c) The standard deviation measures the amount of variation in a variable. Higher values indicate greater spread.",
      "\n")
  
  plot(summ.tab$CV, pch=20, col="darkblue", ylab="mean/sd", 
       xlab="", main="d) Coefficient of Variation (CV)", 
       xaxt="n", cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(1), col=c("lightgrey"), lwd=2, lty=2)
  
  cat("*******************************************************************", 
      "\n")
  cat("d) The coefficient of variation (CV) is the ratio of the standard deviation to the mean. it shows the variability, as defined by the standard deviation, relative to the mean. The coefficient of variation should typically only be used for data measured on a ratio scale. That is, the data should be continuous and have a meaningful zero. Measurement data in the physical sciences and engineering are often on a ratio scale. ",
      "\n")
  
  plot(summ.tab$NAs.per, pch=20, col="darkblue", ylab="Percentage", 
       xlab="", main="e) Missing Values (%)", ylim=c(0,100), xaxt="n",
       cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(90), col=c("lightgrey"), lwd=2, lty=2)
  
  cat("*******************************************************************", 
      "\n")
  cat("e) Missing values are problematic for many statistical techniques. Here we provide a relative measure using the percentage for each variable.",
      "\n")
  
  plot(summ.tab$ZEROs.per, pch=20, col="darkblue", ylab="Percentage", 
       xlab="", main="f) Values less than or equal to zero (%)", ylim=c(0,100), xaxt="n", 
       cex=sym.cex)
  axis(side=1, at=1:length(vars), labels=var.labs, las=2, cex.axis=lab.cex)
  abline(h=c(90), col=c("lightgrey"), lwd=2, lty=2)
  
  cat("*******************************************************************", 
      "\n")
  cat("f) Values less than or equal to Zero can impact data transformations and result in biased models. Here we provide a relative measure using the percentage for each variable.",
      "\n") 
 
  #return(summ.tab)
  cat("*******************************************************************", 
      "\n")
  cat("Coefficient of Variation (CV) Rankings provide insight into which variables in the data exhibit the most dispersion. This can help understand the potential for a variable to drive 'mixture' variability.",
      "\n")
  print.tab=summ.tab[,c("VARIABLE", "N", "AVG", "SD", "MN", "Q1", "Q3", "MX",  "CV", "CVrank", "MED", "IR")]
  print.tab2=print.tab[order(CVrank),]
  print(print.tab2)
}
############################################################################################

############################################################################################
#Explore pairwise correlations  
cor.eval<-function(data, cor.method="pearson", lab.cex=1){
  
  #Set up data
  data=data
  vars=names(data)
  
  #Set plot specifics

  var.labs=ifelse(nchar(vars)<8,vars,substr(vars, 0,8))
  lab.cex=lab.cex
  
  par(mfrow=c(1,1), mar=c(5,4,2,1), family='serif', cex=lab.cex, ask=TRUE, pty="m")
  
  #Examine group structure via distribution of pairwise distances between objects
  #Identy pairwise correlations
  cor.mat<-cor(data, method=cor.method, use="pairwise.complete.obs")
  
  write.table(x=cor.mat, file="cor.summ.tab.csv", sep=",")
  cat("Pairwise Correlation Summary Measures",
      "\n")
  cor.summ=data.frame(cor.mat)
  cor.summ$r.avg=rowMeans(abs(cor.mat))
  cor.summ$r.rank=rank(-cor.summ$r.avg)
  cor.summ2<-cor.summ[order(cor.summ$r.rank),]
  print(cor.summ2[,c("r.avg", "r.rank")])
  
  hist(cor.mat, main="Pairwise Correlations", xlab=cor.method, 
       cex=2, xlim=c(-1,1), col="lightgrey") 
  abline(v=c(-.8,0,.8), col=c("darkgrey","black", "darkgrey"), lwd=c(2), lty=c(2,1,2))
  box()
  box(which="outer")
  
  col.r <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(cor.mat, method="color", col=rev(col.r(200)),type="lower",
           order="hclust", tl.cex=lab.cex,
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
  print(summary(pca.mix))
  
  pca.loadings<-round(pca.mix$rotation,2)
  
  cat("*******************************************************************", 
      "\n")
  cat("Component Loadings",
      "\n")
  print(pca.loadings)
  write.table(x=pca.loadings, file="pca.loadings.summ.csv", sep=",")
  

  #Explore the variance explained by each component
  pr.var=pca.mix$sdev^2
  #Calculate the proportion of variance explained
  pve=round(pr.var/sum(pr.var),2)

  #Plot the component Standard deviations to see which are useful 
  plot(pca.mix$sdev^2, pch=20, ylab="Psuedo Eigenvalue", xlab="Principal Component",
     cex=lab.cex, type="b", main="a) PCA Scree Plot", col="darkblue")
  abline(h=1, col="darkgrey", lty=2)
  box(which="figure")

  plot(pve, xlab="Principal Component", ylab="Proportion of Variance", ylim=c(0,1), type='b', 
       pch=20, col="darkblue", main="b) Variance Explained", cex=lab.cex)
  points(cumsum(pve), type='b',pch=20, col="darkred")
  abline(h=0.85, col="darkgrey", lwd=2, lty=2)
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

############################################################################################
#Explore grouping structure (i.e., clustering) using kmeans and multiple cluster statistics 
grp.eval<-function(data, kmn, kmx){
  
  #set Data
  data=data.matrix(data)
  class(data)
  
  #Examine group structure using kmeans and common cluster statistics
  CLUST.EVAL<-NbClust(data, diss=NULL, distance = "euclidean", min.nc=kmn, max.nc=kmx, 
                      method = "kmeans", index=c("alllong")) 
  
  #CLUST.EVAL
  CH=CLUST.EVAL$All.index[,"CH"]
  SW=CLUST.EVAL$All.index[,"Silhouette"]
  GAP=CLUST.EVAL$All.index[,"Gap"]
  
  chb=CLUST.EVAL$Best.nc[1,"CH"]
  swb=CLUST.EVAL$Best.nc[1,"Silhouette"]
  gapb=CLUST.EVAL$Best.nc[1,"Gap"]
  
  par(mfrow=c(3,2), mar=c(5,4,2,1), cex=1, family='serif')
  
  #Examine group structure via distribution of pairwise distances between objects
  dist.mat<-dist(x=data, method="euclidean")
  
  hist(dist.mat, main="a) Pairwise Object Distances", xlab="Euclidean Distance", 
       cex=2, col="lightgrey")
  box()
  box(which="outer")
  
  #Apply Hierarchical clustering to visualize cluster dendrogram
  h.clust <- hclust(dist.mat, method = "ward.D")
  plot(h.clust, cex=0.1, main="b) Ward's Dendrogram", hang=-1, xlab="", col="darkgrey", sub="")
  box()
  box(which="outer")
  
  #Apply multidimensional scaling to visualize clustering structure
  loc <- cmdscale(dist.mat)
  x <- loc[, 1]
  y <- -loc[, 2] # reflect so North is at the top
  ## note asp = 1, to ensure Euclidean distances are represented correctly
  plot(x, y, type = "p", xlab = "", ylab = "", asp = 1, axes = FALSE,
       main = "c) Multi-dimensional Scaling", pch=20, col="darkblue")
  #text(x, y, rownames(loc), cex = 0.6)
  box()
  box(which="outer")
  
  ##########################################################################
  
  #Create Plots
  plot(x=seq(kmn,kmx,1), y=CH, pch=20, col="darkblue", xlab="Number of clusters (k)",
       ylab="Pseudo F-Statistic", main="d) Calinski Harabaz ",type="b")
  abline(v=chb, col="darkgrey", lty=2)
  
  plot(x=seq(kmn,kmx,1), y=SW, pch=20, col="darkblue", xlab="Number of clusters (k)",
       ylab="Width", main="e) Silhouette Width ",type="b", ylim=c(-1,1))
  abline(v=swb, col="darkgrey", lty=2)
  
  plot(x=seq(kmn,kmx,1), y=GAP, pch=20, col="darkblue", xlab="Number of clusters (k)",
       ylab="Gap", main="f) GAP Statistic",type="b")
  abline(v=gapb, col="darkgrey", lty=2)
  
  CLUST.STATS=data.frame(k=seq(kmn,kmx,1), Calinski=CH, Silhouette=SW, Gap=GAP)
  
  cat("*******************************************************************", 
      "\n")
  cat("*******************************************************************", 
      "\n")
  cat("Uncovering structure in uknown data with unsuperised learning tools is challenging and determining the number of groups in the data is often a difficult decision. The grp.str.eval() function is designed to provide as strategy for discovering potential clustering/natural grouping in your data if it exists. To achieve, the evaluation includes three commonly applied visualation techniques and three well respected cluster statistics.",
      "\n")
  cat("*******************************************************************", 
      "\n")
  
  cat("Graphical views of multivariate data can provide insights into the structure of the data.",
      "\n")
  cat("*******************************************************************", 
      "\n")
  cat("Panel a) The histogram can be used to search for modes in the pairwise object distances. Multiple modal peaks often suggest clustering.",
      "\n")
  cat("*******************************************************************", 
      "\n")
  cat("Panel b) provides a dendrogram that illustrates hierarchical clustering results as a tree.",
      "\n")
  cat("*******************************************************************", 
      "\n")
  cat("Panel c) Multidimensional scaling provides a visualization that supports examination of structure in the pairwise disimilarity data. Here distances between points reflect object similarity, where neighbors are more similar." ,
      "\n")
  
  cat("*******************************************************************", 
      "\n")
  
  cat("Clustering statistics provide internal clustering criteria that are often used to heuristicly determing cluster number.",
    "\n")
  cat("*******************************************************************", 
      "\n")
  cat("Panel d) The Calinski-Harabasz score, also known as the Variance Ratio Criterion, is defined as the ratio between within-cluster dispersion and between-cluster dispersion. Higher values indicate a 'better' solution.",
    "\n")
  cat("*******************************************************************", 
      "\n")
  cat("Panel e) The silhouette value is a measure of how similar an object is to its own cluster (cohesion) compared to other clusters (separation). The silhouette ranges from ???1 to +1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters.",
    "\n")
  cat("*******************************************************************", 
      "\n")
  cat("Panel f) The gap statistic compares the total within intra-cluster variation for different values of k with their expected values under null reference distribution of the data. The estimate of the optimal clusters will be value that maximize the gap statistic (i.e, that yields the largest gap statistic). This means that the clustering structure is far away from the random uniform distribution of points.",
    "\n")
  cat("*******************************************************************", 
      "\n")
  cat("*******************************************************************", 
      "\n")
  print(CLUST.STATS)
}
############################################################################################

############################################################################################
############################################################################################
#Functions to evaluate application of SOM
############################################################################################
#Find optimal seed value 
som.seed<-function(data, nstart=10, iter.max=500*(somx*somy), somx, somy, grid.topo="rectangular"){
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
  
  cat("Sum of Squares Summary:",
      "\n")
  
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
  cat("Adjusted R2:",
      "\n")
  
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
  
  cat("Pseudo.Fstat:",
      "\n")
  
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
  
  cat("AIC:",
      "\n")
  
  return(AIC)
}  
  

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
  
  mod.stats=data.frame("R2"=somr2, "AIC"=somaic, "PseudoF"=somfstat)
  print(mod.stats)
  
  
  #Extract Training Data
  data=data.frame(som.obj$data)
  var.names<-names(data)
  
  #Extract SOM dimensions and create labels
  somx=som.obj$grid$xdim
  somy=som.obj$grid$ydim
  som.k=somx*somy
  XYs=paste(round(som.obj$grid$pts[,1],2),round(som.obj$grid$pts[,2],2), sep="")
  IDs=1:som.k
  #SOM Coordinates and reference table
  som.coords<-data.frame(X=som.obj$grid$pts[,1],Y=som.obj$grid$pts[,2])
  som.class.xy<-data.frame(ID=factor(IDs), XY=factor(XYs),X=as.numeric(som.obj$grid$pts[,1]),
                           Y=as.numeric(som.obj$grid$pts[,2]))
  
  #SOM Profiles
  profiles<-data.frame(som.obj$codes)
  row.names(profiles)<-IDs
  
  #print("Map Profiles and Coordinate Linkage Complete")
  
  #Create data table of SOM predictions 
  #Classifications
  som.class<-som.obj$unit.classif
  class.tab<-data.frame(OBS=1:length(som.class), ID=som.class)
  
  class.tab2<-merge(class.tab, som.class.xy, by="ID")
  class.tab2<-class.tab2[,c("OBS", "ID","XY","X","Y")]
  class.tab3<-class.tab2[do.call(order,class.tab2),]
  class.tab3$ERROR<-som.obj$distances
  rownames(class.tab3)<-class.tab3$OBS
  
  profiles2<-data.frame(ID=rownames(profiles), profiles)
  som.preds<-merge(class.tab3, profiles2, by="ID",all.x=TRUE)
  
  som.variable<-som.preds[order(som.preds$OBS),]
  
  #print('SOM Predictions Table Complete')
  
  #Frequency Summary
  som.n<-table(class.tab$ID)
  som.freq<-round(100*(table(class.tab$ID)/sum(table(class.tab$ID))),2)
  
  freq.tab<-data.frame(som.n, som.freq)
  freq.tab2<-merge(som.class.xy, freq.tab, by.x="ID", by.y="Var1", all.x=TRUE)
  freq.tab2$Var1.1<-NULL
  colnames(freq.tab2)<-c("ID","XY","X","Y","N","FREQ")
  freq.tab2<-freq.tab2[order(as.numeric(freq.tab2$ID)),]
  
  cat("*******************************************************************", 
      "\n")
  cat("SOM Profile Coordinates and Frequencies",
      "\n")
  print(freq.tab2)
  
  ####################################################################
  #Internal evaluation of SOM Fit using distance measures
  QE<-mean(som.obj$distances)
  error<-aggregate(x=som.obj$distances, by=list(som.obj$unit.classif), FUN=mean, na.omit=TRUE)
  error.tab<-data.frame(error) 
  colnames(error.tab)<-c("ID", "WITHIN")
  error.tab2<-merge(som.class.xy,error.tab, by.x="ID", by.y="ID", all.x=TRUE)
  profile.dist<-data.matrix(dist(data.frame(som.obj$codes)))
  profile.dist2<-apply(profile.dist, MARGIN=1, FUN=mean, na.rm=TRUE)
  error.tab2$BETWEEN<-profile.dist2
  error.tab2$WB.RATIO=error.tab2$WITHIN/error.tab2$BETWEEN
  error.tab3<-error.tab2[order(as.numeric(error.tab2$ID)),]
  
  #Extract cluster' evaluation statistics for each SOM profile 
  dist.summ<-error.tab3
  #print("Profile Distance Evaluation Complete")
  cat("*******************************************************************", 
      "\n")
  cat("Internal Class Distance Measures for SOM Profiles ",
      "\n")
  print(dist.summ)
  
  ####################################################################
  #Construct Evaluation Data Sets
  data.eval=data.frame(class.tab, data)
  data.eval.xy<-merge(data.eval, som.class.xy, by="ID")
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
map.size.eval<-function(data, kmx, nstart=10, iter.max=10000, grid.topo="rectangular"){
  
  #Set common dimensions
  som.x=c(2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20)
  som.y=c(2,2,3,3,4,4,5,5,6,6,7,7,8,8,9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20)
  som.xy=som.x*som.y
  
  max.size<-which(abs(som.xy-kmx)==min(abs(som.xy-kmx)))
  
  #set Data
  data=data
  class(data)
  var.names=names(data.frame(data))
  n.vars=length(var.names)
  n.obs=dim(data)[1]
  
  MOD.TAB<-data.frame()
  EVAL.TAB<-data.frame()
  FREQ.TAB<-data.frame()
  MAP.TAB<-data.frame()
  CH.TAB<-data.frame()
  
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
                       QE=as.numeric(mean(som.mod$distances)),
                       AIC=as.numeric(som.summ$MODEL_STATS$AIC), 
                       PSEUDO.F=as.numeric(som.summ$MODEL_STATS$PseudoF))
       
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
    
    ch.tab=data.frame(SOMX=somx, SOMY=somy, k=somx*somy,
                      CH=round(calinhara(data,som.mod$unit.classif),digits=2))
      
    MOD.TAB<-rbind(MOD.TAB,mod.tab)
    EVAL.TAB<-rbind(EVAL.TAB, eval.tab)
    FREQ.TAB<-rbind(FREQ.TAB,freq.tab) 
    MAP.TAB<-rbind(MAP.TAB, map.tab)
    CH.TAB<-rbind(CH.TAB, ch.tab)
  }
  
  par(mfrow=c(4,2), mar=c(4,4,1.5,1), family='serif', ask=FALSE)
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("Determination of map size is an important aspect of SOM-based data analyses. Applications of the SOM algorithm varies broadly and thus we note that the measures provided here are generally designed for analyses of complex environmental exposure data. The map.size.eval() function presents six measures designed to assist with selection of map size (i.e., k). In brief, panels a and b provide common overall measures of model fit and for panels c-f boxplots are used to show distibutions evaluation metrics across all variables for each map size. We note that there is often no right answer here. That said, 'elbowing' or threshold based strategies are suggested for determinng optimal map size.",
      "\n")
  
  #Apply multidimensional scaling to visualize overall data structure
  dist.mat<-dist(data, method="euclidean")
  loc <- cmdscale(dist.mat)
  x <- loc[, 1]
  y <- -loc[, 2] # reflect so North is at the top
  ## note asp = 1, to ensure Euclidean distances are represented correctly
  plot(x, y, type = "p", xlab = "", ylab = "", asp = 1, axes = FALSE,
       main = "a) Multi-dimensional Scaling", pch=20, col="darkblue", cex=1.5)
  #text(x, y, rownames(loc), cex = 0.6)
  box()
  box(which="outer")
  
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("a)  A Multi-dimensional Scaling Plot",
      "\n")
  
  #PCA Plot
  pca.mix<-prcomp(data, center=FALSE, scale=FALSE)
  comp.SD<-pca.mix$sdev
  comp.n<-length(which(comp.SD > 0.99))
  comp.n
  
  #Explore the variance explained by each component
  pr.var=pca.mix$sdev^2
  #Calculate the proportion of variance explained
  pve=round(pr.var/sum(pr.var),2)
  
  plot(pve, xlab="Principal Component", ylab="Proportion of Variance", ylim=c(0,1), type='b', 
       pch=20, col="darkblue", main="b) Variance Explained", cex=2)
  points(cumsum(pve), type='b',pch=20, col="darkred")
  abline(v=comp.n, col="darkgrey", lwd=2, lty=2)
  legend("bottomright", pch=20, col=c("darkblue", "darkred"), legend=c("Individual", "Cumulative"))
  
  
  #Calinski Harabasz
  plot(CH.TAB$CH~CH.TAB$k, ylab="Calinski_Harabasz Index", xlab="Number of profiles (k)", 
       main="c) Grouping Structure (higher is better)", type="b", pch=19, cex=2, col="darkgrey")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("c)  The Calinski-Harabasz Index is an internal validation measure used as a cluster statistic to aid in determing number of clusters. Higher values indicate better group structure.",
      "\n")
  
  plot(MOD.TAB$R2~MOD.TAB$k, ylab="Adjusted R2", xlab="Number of profiles (k)", 
       main="d) SOM Model Adjusted R-squared", type="b", pch=19, cex=2, col="darkgrey")
  
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("d)  The R2 reflects the proportion of variation in the response variables explained by the SOM model at each map size. In this setting this measure reflects how well the mapping performed as a means of 'dimension reduction' aimed at retaining the variablity of the original data.",
      "\n")
  
  boxplot(EVAL.TAB$R.P~EVAL.TAB$k, ylab="R", xlab="Number of profiles (k)",
          main="e) SOM Profile-Variable Correlations")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("e) Profile-variable R values reflect association between SOM profiles at each map size for each variable in the dataset.",
      "\n")
  
  boxplot(EVAL.TAB$RMSE.P~EVAL.TAB$k, ylab="RMSE", xlab="Number of profiles (k)",
          main="f) SOM Profile-Variable RMSE")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("f) Profile-Varialbe RMSE reflect how accurately SOM profiles at each map size align with observations for each variable in the dataset.",
      "\n")
  
  boxplot(FREQ.TAB$N~FREQ.TAB$k, ylab="N", xlab="Number of profiles (k)",
          main="g) SOM-Profile Frequency Distributions")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("g) Profile frequency distributions reflect sample size of profile assignments for each map size. ",
      "\n")
  
  boxplot(MAP.TAB$MI~MAP.TAB$k, ylab="Moran's I", xlab="Number of profiles (k)",
          main="h) SOM Grid-Variable Spatial Autocorrelation")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("h) Variablewise assessments of spatial autocorrelation for each map size reflect organizational structure of mapping.If the values in the dataset tend to cluster spatially (high values cluster near other high values; low values cluster near other low values), the Moran's Index will be positive. When high values repel other high values, and tend to be near low values, the Index will be negative.",
      "\n")
  
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("--------------------------------------------------------------------------------",
      "\n")
  cat("--------------------------------------------------------------------------------",
      "\n")
 
  MAP.EVAL<-list(MOD.FIT=MOD.TAB, VAR.FIT=EVAL.TAB, MAP.FREQ=FREQ.TAB, MAP.STRUCTURE=MAP.TAB, GRP.STR=CH.TAB)
  #return(MAP.EVAL)
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
map.profile.plot<-function(som.obj, lab.cex=1, legend.loc="bottomright", 
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
          circles = rep(.5, nrow(som.summ$COORDINATES)),
          inches = FALSE, add = TRUE,
          fg = "black", bg = NA)
  
  text(x=som.summ$COORDINATES$X, y=som.summ$COORDINATES$Y+.4,
       labels=paste("[", IDs ,"]", sep=""), font=2, cex=lab.cex+.5)
  
  #text(x=som.summ$COORDINATES$X+.35, y=som.summ$COORDINATES$Y-.3,
  #     labels=paste(round(som.summ$FREQUENCIES$FREQ,1) ,"%", sep=""), font=2.5, cex=lab.cex)
  
  #text(x=som.summ$COORDINATES$X+.35, y=som.summ$COORDINATES$Y-.4,
  #     labels=paste("n=", som.summ$FREQUENCIES$N, sep=""), font=2.5, cex=lab.cex)
  
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
map.bar.plot<-function(som.obj, lab.cex=1, label.loc="bottomright", legend.cex=1, 
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
  for(i in 1:(somx*somy)){
    plot.vals<-as.numeric(plot.data[i,])
    plot.err<-as.numeric(plot.error[i,])
    
    bp<-barplot(plot.vals, col=color.mod, axes=TRUE, axisnames=FALSE,
            ylim=c(1.2*min(plot.data), 1.2*max(plot.data)))
    error.bar(bp,plot.vals,plot.err)
    abline(h=0, lty=2, col="black")
    title(main=paste("[", i ,"]", sep=""), font=2, cex.main=lab.cex+.5)
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
map.prof.comp.plot<-function(som.obj, col.mod=diverging_hcl){
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
         main=var.labs[i], palette.name=color.mod)
  }
  
  plot(som.obj, type="property", property=rowSums(getCodes(som.obj,1)), 
       main="Cumulative Sum", palette.name=color.mod)
 
}
############################################################################################
############################################################################################
#Map Fit Evaluation Functions
map.fit.plots<-function(som.obj, lab.cex=1){
  #Extract profiles
  som.summ<-som.fit.summ(som.obj)
  #Get SOM Dimensions
  somx<-som.obj$grid$xdim
  somy<-som.obj$grid$ydim
  
  
  par(mfrow=c(3,2), mar=c(5,4,2,1),ask=FALSE, family='serif', pty="m")
  ## Look at within cluster spread
  plot(som.obj, type="quality", main="a) Mean Assignment Distance", palette.name=diverge_hcl)
  #Look at between cluster spread
  plot(som.obj, type="dist.neighbours", main="b) Sum of Nearest Neighbor Distances (i.e., U-Matrix)", palette.name=diverge_hcl)
  
  ## SOM-Variable Evaluation
  mapvareval<-som.summ$COMPONENT_EVALUATION
  barplot(height=mapvareval$R, names.arg=mapvareval$Variable, 
          las=2, cex.names=lab.cex, main="c) Overall Profile-Component Correlations", ylab="Profile-Variable Correlation (R)", ylim=c(-1,1), col="lightgrey")
  abline(h=mean(mapvareval$R, na.rm=TRUE), col="darkgrey", lty=2)
  box()
  
  barplot(height=mapvareval$RMSE, names.arg=mapvareval$Variable, 
          las=2, cex.names=lab.cex, main="d) Overall Componentwise RMSE of SOM Profiles", ylab="Profile-Variable RMSE", ylim=c(0,(max(mapvareval$RMSE)*1.2)), col="lightgrey")
  abline(h=mean(mapvareval$RMSE, na.rm=T), col="darkgrey", lty=2)
  box()
  
  ##Map Organization
  maporganization<-som.summ$MAP_ORGANIZATION
  barplot(height=maporganization$'Morans I', names.arg=maporganization$Variable, 
          las=2, cex.names=lab.cex, main="e) Componentwise Spatial Structure on SOM Grid", ylab="Moran's I",
          ylim=c(min(maporganization$'Morans I'),max(maporganization$'Morans I')), col="lightgrey")
  abline(h=mean(maporganization$'Morans I', na.rm=T), col="darkgrey", lty=2)
  box()
  
  #Generate a Sammon's Map to look at Map Distortion
  mapdist<-sammon(d=dist(data.matrix(som.summ$PROFILES)), y=data.matrix(som.summ$COORDINATES[,c("X","Y")]),
                  k=2, niter=1000)
  plot(mapdist$points, pch=20, col="lightgrey", cex=5, xlab="Sammons X", ylab="Sammons Y", main="f) Distance Preserving SOM Grid Map")
  text(mapdist$points, labels = as.character(1:nrow(som.summ$PROFILES)))
  
}  
############################################################################################

############################################################################################
#Extract exposure variable from SOM analyses
############################################################################################
exp.metric<-function(som.obj){
  #Extract Data
  som.data<-data.frame(som.obj$data)
  IDs<-row.names(som.data)
  
  #Extract class assignments
  unit.classif<-som.obj$unit.classif
  class.tab<-data.frame(ROW.IDs=IDs,SOM.CLASS=unit.classif)
  
  #Coordinates
  XYs<-data.frame(SOM.CLASS=1:length(som.obj$grid$pts[,1]), 
                  SOMX=som.obj$grid$pts[,1],
                  SOMY=som.obj$grid$pts[,2])
  
  exp.tab<-merge(class.tab, XYs, by="SOM.CLASS", all.x=TRUE)
  exp.tab$SOM.ERROR<-som.obj$distances
  
  return(exp.tab)
}

############################################################################################
#Plot Function
freq.bar.plot<-function(som.obj, freq.tab, lab.cex=1, label.loc="bottomright", legend.cex=1)
{
  #Get summary output
  som.summ<-som.fit.summ(som.obj)
  #Get SOM Dimensions
  somx<-som.obj$grid$xdim
  somy<-som.obj$grid$ydim
  #Get SOM Codes
  profiles<-as.data.frame(som.obj$codes)
  profiles
  
  #Get Covariate Frequencies
  freq.table=freq.tab
  
  #Get variable names
  var.names<-rownames(freq.table)
  n.col=length(var.names)
  col.n=n.col
  
  
  #Set Bar Color Model
  color.mod<-grey.colors(col.n)
  
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
  
  plot.data<-freq.tab
  plot.n<-round(freq.tab/sum(freq.tab),2)
  for(i in 1:(somx*somy)){
    plot.vals<-as.numeric(plot.data[,i])
    
    bp<-barplot(plot.vals, col=color.mod, axes=TRUE, axisnames=FALSE,
                ylim=c(1.2*min(plot.data), 1.2*max(plot.data)))
    title(main=paste("[", i ,"]", sep=""), font=2, cex.main=lab.cex+.5)
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



#Create a function to plot Relative Risks from statistical model on SOM Grid
map.risk.plot<-function(som.obj, mod.fit, exp.var.nam, ref, title.lab=NULL){
  somx<-som.obj$grid$xdim
  somy<-som.obj$grid$ydim
  
  #Generate SOM Coordinate Reference
  xy.coords=data.frame(som.obj$grid$pts)
  xy.coords
  xy.coords$xy<-paste(xy.coords$x, xy.coords$y, sep=',')
  xy.coords$CLASS=1:(somx*somy)
  
  ref=ref
  #som.xy<-som.exp.variable(som.obj, ref)
  
  #Pull model summary
  mod.summ<-summary(mod.fit)
  mod.summ
  mod.coef<-mod.summ$coefficients
  
  
  #Identify exposure estimate labels
  ID<-rownames(mod.coef)
  ID
  ID.som<-grep(exp.var.nam, ID)
  ID.som
  ID[ID.som]
  
  #Extract Risk Estimates
  B<-round(mod.coef[,1],0)
  B
  
  SE<-round(mod.coef[,2],3)
  SE
  
  UL<-round(mod.coef[,1] + 1.96*mod.coef[,2],0)
  LL<-round(mod.coef[,1] - 1.96*mod.coef[,2],0)
  
  
  RR<-round(exp(mod.coef[,1]),3)
  RR
  
  RR_lb<-round(exp(mod.coef[,1] - 1.96*mod.coef[,2]),3)
  RR_lb
  
  RR_ub<-round(exp(mod.coef[,1] + 1.96*mod.coef[,2]),3)
  RR_ub
  
  p.val<-mod.coef[,4]
  p.val
  
  p.val.adj<-round(p.adjust(p.val, method="fdr"),5)
  p.val.adj
  
  risk.tab<-data.frame(ID, B, SE, LL, UL, RR, RR_lb, RR_ub, p.val, p.val.adj)
  print("Tab1")
  print(risk.tab)
  
  risk.tab$SOM_ID<-"NA"
  risk.tab$SOM_ID[ID.som]<- substring(risk.tab$ID[ID.som],nchar(exp.var.nam)+1, nchar(exp.var.nam)+2)
  print("Tab2")
  print(risk.tab)
  
  risk.tab.xy<-merge(xy.coords, risk.tab, by.x="CLASS", by.y="SOM_ID", 
                     all.x=TRUE, all.y=TRUE)
  
  risk.tab.xy<-na.omit(risk.tab.xy)
  print("Tab3")
  print(risk.tab.xy)
  #str(risk.tab.xy)
  
  #Create Plot
  print("step1")
  par(mfrow=c(1,1), pty="s")
  stars(som.obj$codes[[1]], draw.segments=TRUE, locations=som.obj$grid$pts, len=.375)
  print("step2")
  symbols(x=som.obj$grid$pts[,1], y=som.obj$grid$pts[,2], 
          circles=rep(0.5, length(som.obj$grid$pts[,1])),inches=FALSE, add=TRUE, bg="lightgrey", fg="lightgrey")
  print("step3")
  symbols(x=risk.tab.xy$x, y=risk.tab.xy$y, 
          circles=rep(0.5, length(risk.tab.xy$x)),inches=FALSE, add=TRUE, 
          bg=ifelse(risk.tab.xy$B < 0, "lightblue", "salmon"),
          fg=ifelse(risk.tab.xy$p.val.adj<0.059,"black","lightgrey"),
          lwd=2.5)
  print("step4")
  #Add labels
  text(risk.tab.xy$x, risk.tab.xy$y+0.35, cex=1,labels=paste("[",risk.tab.xy$CLASS,"]", sep=""), 
       font=ifelse(risk.tab.xy$p.val.adj<0.05,2,1))
  #Add Betas
  text(risk.tab.xy$x, risk.tab.xy$y+0.15, cex=1.25,labels=round(risk.tab.xy$B,2), 
       font=ifelse(risk.tab.xy$p.val.adj<0.05,2,1))
  print("step5")
  #Add Confidence intervals
  text(risk.tab.xy$x, risk.tab.xy$y, cex=1,
       labels=paste("(",round(risk.tab.xy$LL,2),",",round(risk.tab.xy$UL,2),")",sep=""), 
       font=1)
  print("step6")
  #text(risk.tab.xy$x, risk.tab.xy$y-0.15, cex=0.75,
  #     labels=paste("p=",round(risk.tab.xy$p.val.adj,5),sep=""), font=1)
  title(main=title.lab)
  rm(som.obj)
}












