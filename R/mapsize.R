############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################
###########################################################################################
#'Determining SOM size
#'
#'mapsize provides diagnostic plots and summaries of select criteria for determining how the size of the SOM influences characteristics of the SOM model.
#'@param x is a data object to train the map. Should be a numerical or factor based data matrix.
#'@param kmn is the minimum number of map units (aka nodes) to evaluate. Default value is 2.
#'@param kmx is maximum number of map units (aka nodes) to evaluate. Default is 5*sqrt(n).
#'@param itermax. maximum number of iterations passed to sommix. Default is 500*number of map nodes (k).
#'@param nstarts. number of random initializatons passed to sommix. Default is 5.
#'@param maptopo. map topology passed to sommix. Default is rectangular.
#'@param distmet. distance method passed to sommix. Default is Euclidean.
#'@param lmode. m initializatons passed to sommix
#'@param symsize sets symbol size on plots
#'@return Panels a-d on the diagnostic plot illustrate common model perfomance metrics as a function of map size. Panels e-f examine within-class-sum-of-squares and frequency distributions as a function of map size. A list of model and class-level perfomance statistics is also returned.
#'  \itemize{
#'  \item {SOMX} x dimension
#'  \item {SOMY} y dimension
#'  \item {K} Number of map nodes
#'  \item {R2} R2
#'  \item {ADJ_R2} Adjusted R2
#'  \item {MAE} mean absolute error based on class assignment distances
#'  \item {RMSE} root-mean-square-error based on class assignment distances
#'  \item {AIC} a form of Akaikes Information Criteria applied to clustering algorithms
#'  \item {TotWCSS}  Total Within-Cluster Sum-of-Squares
#'  \item {N} Number class assignments
#'  \item {FREQ} Proportion of class assignments
#'  \item {WCD} average within-class distance
#'  \item {BCD}  average between-class distacne
#'  \item {WB_Ratio}  WCD/BCD
#'
#' }
#'@details An important step in the application of SOM are the user provided inputs for the dimensions of the mapping (i.e., size). Here we provide common model performance metrics and class-level evaluations in effort to assist the user in determining an appropriate map size.
#'@export
#'@examples
#'#NIEHS Mixtures Workshop dataset1
#'data(dataset1)
#'mapsize(scale(dataset1[,3:9]), kmx=10, itermax.=10)


#Develop a set of functions to evaluate a range of common Self-Organizing Map sizes
mapsize<-function(x, kmn=NULL, kmx=NULL, itermax.=NULL, nstarts.=NULL, maptopo.=NULL, distmet.=NULL,
                  lmode.=NULL, symsize=1) {

  #set Data
  data<-x
  varnames<-names(data.frame(data))
  nvars<-length(varnames)
  nobs<-dim(data)[1]

  #Set training data
  data.trn<-data.matrix(data)
  #set Data distance matrix
  dist.mat<-stats::dist(x=data.trn, method="euclidean")

  #Set function defaults for sommix
  if (is.null(kmn)) kmn <- 2
  if (is.null(kmx)) kmx <- 5*sqrt(nobs)

  #Set common dimensions
  som_x<-c(2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20)
  som_y<-c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20)
  som_k<-som_x*som_y
  som_k

  #Set test range of sizes
  minsize<-which(abs(som_k-kmn)==min(abs(som_k-kmn)))
  maxsize<-which(abs(som_k-kmx)==min(abs(som_k-kmx)))

  #Set function defaults for sommix
  if (is.null(maptopo.)) maptopo. <- "rectangular"
  if (is.null(distmet.)) distmet.<-"euclidean"
  if (is.null(lmode.)) lmode. <- "online"
  if (is.null(nstarts.)) nstarts. <- 5

  MOD_EVAL<-data.frame()
  CLASS_EVAL<-data.frame()

  for (i in minsize:maxsize){
    somx<-som_x[i]
    somy<-som_y[i]
    somk<-somx*somy

    if (is.null(itermax.)) itermax. <- 500*(somx*somy)

    som_mix<-sommix(X=data.trn, somx=somx, somy=somy, itermax=itermax., nstarts=nstarts.,
                    maptopo=maptopo., distmet=distmet., lmode=lmode.)

    #Apply sommix summary function
    som_summ<-sommix_summ(som_mix)

    #Assess model fit
    som_fit<-somFitStats(som_mix)

    mod_eval<-data.frame(SOMX=somx, SOMY=somy, k=somk,
                       R2=as.numeric(som_fit$R2),
                       ADJ_R2=as.numeric(som_fit$ADJ_R2),
                       MAE=as.numeric(som_fit$DISTANCE_MAE),
                       RMSE=as.numeric(som_fit$DISTANCE_RMSE),
                       AIC=as.numeric(som_fit$AIC),
                       TotWCSS=as.numeric(som_fit$TOTWCSS))


    class_eval<-data.frame(SOMX=somx, SOMY=somy, k=somk,
                           N=as.numeric(som_summ$FREQ$N),
                           FREQ=as.numeric(som_summ$FREQ$FREQ),
                           WCD=som_fit$CLASS_DISTANCES$WCD,
                           BCD=som_fit$CLASS_DISTANCES$BCD,
                           WB_RATIO=som_fit$CLASS_DISTANCES$WB.RATIO,
                           WCSS=som_fit$WCSS)


    MOD_EVAL<-rbind(MOD_EVAL,mod_eval)
    CLASS_EVAL<-rbind(CLASS_EVAL, class_eval)
  }

  #Evaluation Plots
  opar<-graphics::par(mfrow=c(1,1), mar=c(5,4,3,2), pty="s", cex=1, xpd=FALSE)
  graphics::par(mfrow=c(3,2), mar=c(4,4,1.5,1), family='serif', ask=FALSE, xpd=FALSE)

  #Plots
  plot(MOD_EVAL$ADJ_R2~MOD_EVAL$k, ylab="Proportion of Variance", xlab="Number of nodes (k)",
       main="a) Adjusted R2", type="b", pch=19, cex=symsize, col="darkgrey",
       ylim=c(0,1))

  plot(MOD_EVAL$RMSE~MOD_EVAL$k, ylab="Distance", xlab="Number of nodes (k)",
       main="b) Root-Mean-Square-Error (RMSE)", type="b", pch=19, cex=symsize, col="darkgrey",
       ylim=c(0,max(MOD_EVAL$RMSE)*1.2))

  plot(MOD_EVAL$MAE~MOD_EVAL$k, ylab="Distance", xlab="Number of nodes (k)",
       main="c) Mean-Absolute-Error (MAE)", type="b", pch=19, cex=symsize, col="darkgrey",
       ylim=c(0,max(MOD_EVAL$MAE)*1.2))

  plot(MOD_EVAL$AIC~MOD_EVAL$k, ylab="AIC", xlab="Number of nodes (k)",
       main="d) Akaike Information Criterion (AIC)", type="b", pch=19, cex=symsize, col="darkgrey",
       ylim=c(0,max(MOD_EVAL$AIC)*1.2))

  graphics::stripchart(CLASS_EVAL$WCSS~CLASS_EVAL$k, ylab="", xlab="Number of nodes (k)",
          main="e) Within-Class Sum-of-Squares (WCSS)", vertical = TRUE, pch=19, cex=symsize, col="darkgrey")

  graphics::stripchart(CLASS_EVAL$N~CLASS_EVAL$k, ylab="N", xlab="Number of nodes (k)",
          main="f) Class Frequencies", vertical = TRUE, pch=19, cex=symsize, col="darkgrey")


  on.exit(graphics::par(opar))
  on.exit(graphics::layout(1))

  MAP_EVAL<-list(MOD_EVAL=MOD_EVAL, CLASS_EVAL=CLASS_EVAL)
  return(invisible(MAP_EVAL))

}
############################################################################################
