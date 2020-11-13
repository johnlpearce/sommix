############################################################################################
#sommix: A function to apply Kohonen's Self-Organizing Map algorithm
############################################################################################

#'Standard Model Fit Statistics applied to Self-Organizing Maps (SOM) model
#'
#'somFitStats calculates common measures of model fit for a SOM model
#'@param som_obj an object with class 'sommix'
#'@return a dataframe of values
#'@return TSS     The total sum of squares
#'@return WCSS    The within-class sum of squares
#'@return TOTWCSS The total within-class sum of squares
#'@return BCSS    The between-class sum of squares, i.e. TSS-TOTWCSS
#'@return R2     The total variance explained by the SOM classes (BCSS/TSS)
#'@return ADJ_R2 The total variance explained by the SOM classes adjusted for class number
#'@return MAE     The mean absolute error (MAE) derived from the mean class assignment distance, i.e., quantization error
#'@return RMSE    The root-mean-square-error derived from class assignment distances, i.e., quantization error
#'@export

somFitStats<-function(som_obj){
  somx<-som_obj$grid$xdim
  somy<-som_obj$grid$ydim
  somk<-somx*somy

  data<-data.frame(som_obj$data)
  SOM_ID<-som_obj$unit.classif
  data2<-cbind(SOM_ID,data)
  n<-dim(data)[1]

  #Set range of functions to evaluate sum-of-squares a la kmeans
  #Total Sum of Squares
  TSS <- sum(scale(data[,-1], scale = FALSE)^2)

  #Within-class sum of squares
  WCSS=NULL

  for (i in 1:somk){
    datasub1=subset(data, SOM_ID == i)
    wss<-sum(scale(datasub1[,-1], scale = FALSE)^2)
    WCSS<-c(WCSS, wss)
  }

  #Total sum of within-class sum-of-squares
  Tot.WCSS<-sum(WCSS)

  #Between class sum of squares
  BCSS<-TSS-Tot.WCSS

  #calculate reduction in varaince explained by mapping (i.e., R-squared)
  R2<-BCSS/TSS
  ADJ_R2 <- 1-(Tot.WCSS*(n-1))/(TSS*(n-somk))

  ########################################################################
  #Calculate distance-based measures of error
  MAE<-mean(som_obj$distances)
  RMSE<-sqrt(mean(som_obj$distances^2))

  #######################################################################
  #Calculate AIC
  m <- dim(data.frame(som_obj$codes))[2]
  n <- length(som_obj$unit.classif)
  k <- dim(data.frame(som_obj$codes))[1]
  D <- Tot.WCSS
  somAIC<-(D + 2*m*k)


  ########################################################################
  #Internal evaluation of SOM nodes using distance summary measures
  #Within-class distance
  WCD=NULL

  for (i in 1:somk){
    datasub1=subset(data, SOM_ID == i)
    wcd<-mean(stats::dist(datasub1[,-1]))
    WCD<-c(WCD, wcd)
  }

  dist.tab<-data.frame(SOM_ID=1:length(som_obj$grid$pts[,1]),SOMX=som_obj$grid$pts[,1], SOMY=som_obj$grid$pts[,2], k=somk, WCD)

  BCD.dist<-data.matrix(stats::dist(data.frame(som_obj$codes)))
  BCD.dist2<-apply(BCD.dist, MARGIN=1, FUN=mean, na.rm=TRUE)
  dist.tab$BCD<-BCD.dist2
  dist.tab$WB.RATIO=dist.tab$WCD/dist.tab$BCD
  dist.summ<-dist.tab[order(as.numeric(dist.tab$SOM_ID)),]
  dist.summ

  modstats<-list("TOTSS"=TSS, "WCSS"=WCSS,
                "TOTWCSS"=Tot.WCSS, "BCSS"=BCSS,
                "R2"=R2, "ADJ_R2"=ADJ_R2, "DISTANCE_MAE"=MAE, "DISTANCE_RMSE"=RMSE,
                "AIC"=somAIC,"CLASS_DISTANCES"=dist.summ)
  return(modstats)
}










