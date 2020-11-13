############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################
#A collection of plotting functions to support illustration of Self-Organizing Maps (SOM) results, modelled for compatability with sommix.

#Generate radial segment plots for SOM profiles
#'Plots SOM profiles using radial segment diagrams
#'@param som_obj a SOM object
#'@param varnames specifies variable names
#'@param colormod can be used to customize plot color schemes. The default is terrain.colors().
#'@param nodelab whether or not to display reference labels aon each node
#'@param labtype determines type of node referencing to display. "IDs" and "XYs" are available
#'@param labsize sets size of node labels
#'@param addFreq whether or not to display frequency labels
#'@param freqtype determines type of frequency labels. Relative frequencies can be specified with "frq" or counts with "cnts"
#'@param freqsize frequency label cex value
#'@param legsymsize legend symbol cex value
#'@param leglabsize legend text cex value
#'@param legtxtlas legend text position rotation passed to las. Numeric in {0,1,2,3}; the style of axis label.
#'@details Radial segment diagrams provide a compact way to illustrate profiles derived from SOM. Here diagrams are presented along the SOM grid in order provide a compact visualization of SOM results.s the map.
#'@return a kohenen map where radial segment diagrams are used to illustrate profiles for each map node
#'@export

sommix_star<-function(som_obj, varnames=NULL, colormod=NULL,
                      nodelab=TRUE, labtype="IDs", labsize=1,
                      addFreq=FALSE, freqtype="frq", freqsize=1,
                      legsymsize=2, leglabsize=1, legtxtlas=2)
                      {
  #Apply sommix summary function
  som_summ<-sommix_summ(som_obj)

  #Extract info for plotting
  somx<-som_summ$som_grid$xdim
  somy<-som_summ$som_grid$ydim

  profiles<-data.frame(som_summ$som_codes)
  p<-dim(profiles)[2]
  vars<-names(profiles)

  #Set Profile labels
  XYs<-paste("[", round(som_summ$som_coords$SOM_X,1), ",", round(som_summ$som_coords$SOM_Y,1),"]", sep="")
  IDs<-paste("[",som_summ$som_coords$SOM_ID,"]", sep="")


  #Set Frequency labels
  FREQ<-round(som_summ$som_freq$FREQ,1)
  N<-som_summ$som_freq$N

  #Set function defaults
  if (is.null(colormod)) colormod <- grDevices::terrain.colors(n=p)
  if (is.null(varnames)) varnames <- ifelse(nchar(vars)<7,vars,substr(vars, 0,7))

  #Specify label type
  if ( identical(labtype, "IDs")) {
    labtypes <- IDs
  }  else {
    labtypes<- XYs
  }

  #Determine which labels to use
  if ( identical(nodelab, TRUE)) {
    nodelabs<-labtypes
  }  else {
    nodelabs<- ""
  }

  #Specify Frequency labels
  if ( identical(addFreq, TRUE)) {
    freqlabs <- paste(FREQ,"%", sep="")
    nlabs<-paste("n=", N, sep="")
  }  else {
    freqlabs<-""
    nlabs<-""
  }

  #Specify which frequency labels to use
  if ( identical(freqtype, "frq")) {
    freqlab <- freqlabs
    }  else {
    freqlab<-nlabs
  }


  #Create Profile plot
  opar<-graphics::par()
  graphics::par(mfrow=c(1,1),mar=c(.6,.6,2.6,.6))
  graphics::layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(5,1), widths=c(1,1),respect = FALSE)

  graphics::stars(x=profiles, locations=som_summ$som_coords[,3:4],
        draw.segments=TRUE, axes=FALSE, scale=TRUE,
        len = 0.4, col.segments=colormod,
        labels=NULL,ylim=c(0.5,somy+.5), xlim=c(0.5,somx+.5))

  graphics::symbols(som_summ$som_coords[,3:4],
          squares = rep(1, nrow(som_summ$som_coords)),
          inches = FALSE, add = TRUE,
          fg = "black", bg = NA)

  graphics::text(x=som_summ$som_coords$SOM_X, y=som_summ$som_coords$SOM_Y+.425,
       labels=nodelabs, font=2, cex=labsize)

  graphics::text(x=som_summ$som_coords$SOM_X+.375, y=som_summ$som_coords$SOM_Y-.4,
       labels=freqlab,
       font=2, cex=freqsize)

  #Plot Legend on lower panel
  graphics::par(mar=c(rep(1,4)))
  plot(x=1:p, y=rep(.9,p), pch=22, cex=legsymsize, col="black", bg=colormod,
       ylim=c(0,1), xlim=c(1, (p +.25)), axes=FALSE)
  graphics::box(which="outer")
  graphics::axis(side=1, at=1:p, pos=0.85, labels=varnames, cex.axis=leglabsize, las=legtxtlas, tick=FALSE)

  #Reset plot window to normal
  suppressWarnings(graphics::par(opar))
  graphics::layout(1)

}

############################################################################################
############################################################################################
#Bar Plot with error bars of Profiles

#A function to add error bars to a barchart
error.bar <- function(x, y, upper, lower=upper, length=0.01,...){
  graphics::arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, col="darkgrey", ...)
}

############################################################################################
############################################################################################
#Generate bar plots for SOM profiles
#'Plots SOM profiles using bar charts
#'
#'Barplots are used to represent standardized summary values of components under each SOM profile using bars of various lengths.
#'@param som_obj a SOM object
#'@param varnames specifies variable names
#'@param colormod can be used to customize plot color schemes. The default is terrain.colors().
#'@param nodelab whether or not to display reference labels aon each node
#'@param labtype determines type of node referencing to display. "IDs" and "XYs" are available
#'@param labsize sets size of node labels
#'@param addFreq whether or not to display frequency labels
#'@param freqtype determines type of frequency labels. Relative frequencies can be specified with "frq" or counts with "cnts"
#'@param freqsize frequency label cex value
#'@param legsymsize legend symbol cex value
#'@param leglabsize legend text cex value
#'@param legtxtlas legend text position rotation passed to las. Numeric in {0,1,2,3}; the style of axis label.
#'@param barstat specifies summary statistic to plot as bars. "MED" and "AVG" are available
#'@details Barplots are a straighforward way to show patterns in component values under each SOM profile. The size of the bar represents the component summary under each profile. Error bars are also provided.  between obser
#'@return a kohonen plot that uses barplots to describe the profiles for each map node.
#'@export


sommix_bar<-function(som_obj, varnames=NULL, colormod=NULL,
                     nodelab=TRUE, labtype="IDS", labsize=1,
                     addFreq=FALSE, freqtype="frq", freqsize=1,
                     legsymsize=2, leglabsize=1, legtxtlas=2,
                     barstat="MED")
{
  #Apply sommix summary function
  som_summ<-sommix_summ(som_obj)

  #Extract info for plotting
  somx<-som_summ$som_grid$xdim
  somy<-som_summ$som_grid$ydim

  profiles<-data.frame(som_summ$som_codes)
  p<-dim(profiles)[2]
  vars<-names(profiles)

  #Set Profile labels
  XYs<-paste("[", round(som_summ$som_coords$SOM_X,1), ",", round(som_summ$som_coords$SOM_Y,1),"]", sep="")
  IDs<-paste("[",som_summ$som_coords$SOM_ID,"]", sep="")


  #Set Frequency labels
  FREQ<-round(som_summ$som_freq$FREQ,1)
  N<-som_summ$som_freq$N

  #Generate summary statistics for barplots
  summtab<-compstats(som_summ)

  #AVGS
  AVG<-data.frame(summtab$MEAN[,-c(1:4)])
  STDV<-data.frame(summtab$ST_DEV[,-c(1:4)])
  STDV[is.na(STDV)]<-0

  #MEDIANS
  MED<-data.frame(summtab$MED[,-c(1:4)])
  IQRS<-data.frame(summtab$IQRS[,-c(1:4)])
  IQRS[is.na(IQRS)]<-0

  #Set function defaults
  if (is.null(colormod)) colormod <- grDevices::terrain.colors(n=p)
  if (is.null(varnames)) varnames <- ifelse(nchar(vars)<7,vars,substr(vars, 0,7))

  #Specify label type
  if ( identical(labtype, "IDs")) {
    labtypes <- IDs
  }  else {
    labtypes<- XYs
  }

  #Determine which labels to use
  if ( identical(nodelab, TRUE)) {
    nodelabs<-labtypes
  }  else {
    nodelabs<- ""
  }

  #Specify Frequency labels
  if ( identical(addFreq, TRUE)) {
    freqlabs <- paste(FREQ,"%", sep="")
    nlabs<-paste("n=", N, sep="")
  }  else {
    freqlabs<-""
    nlabs<-""
  }

  #Specify which frequency labels to use
  if ( identical(freqtype, "frq")) {
    freqlab <- freqlabs
  }  else {
    freqlab<-nlabs
  }

  #Specify statistics for barplot
  if ( identical(barstat, "MED")) {
    plot.bars <- MED
    plot.error<- IQRS
  }  else {
    plot.bars <- AVG
    plot.error<- STDV
  }


  #Set plot window specifics
  opar<-graphics::par()
  graphics::par(mar=c(.5,1,1,.5), pty="m", cex=1)

  #Generate a layout matrix
  M<-matrix(data=1:(somx*somy),
            nrow=somy, ncol=somx, byrow=TRUE)

  #Flip matrix to match SOM grid
  M1=apply(t(M),1,rev)

  M2<-rbind(M1, rep(1+(somx*somy),somx))

  mf <- graphics::layout(M2, widths = rep(1, ncol(M2)),
               heights = c(rep(1, nrow(M2)-1), 1), respect = FALSE)


  #Run loop to create plots for each SOM node
  for(i in 1:(somx*somy)){

    plot.vals<-as.numeric(plot.bars[i,])
    plot.err<-as.numeric(plot.error[i,])

    bp<-graphics::barplot(height=plot.vals, col=colormod, axes=TRUE, axisnames=FALSE,
                ylim=c(min(plot.bars)+ -1.5, 1.2*max(plot.bars)))
    error.bar(bp,plot.vals,plot.err)
    graphics::abline(h=0, lty=2, col="black")
    graphics::title(main=nodelabs[i], font=2, cex.main=labsize)
    graphics::legend(legend=freqlab[i], cex=freqsize, text.font=2,
           x="bottomright", bty="n", horiz=TRUE)
    graphics::box()
    graphics::box(which="outer")

  }

  #Plot Legend on lower panel
  graphics::par(mar=c(rep(1,4)))
  plot(x=1:p, y=rep(.9,p), pch=22, cex=legsymsize, col="black", bg=colormod,
       ylim=c(0,1), xlim=c(1, (p +.25)), axes=FALSE)
  graphics::box()
  graphics::box(which="outer")
  graphics::axis(side=1, at=1:p, pos=0.85, labels=varnames, cex.axis=leglabsize, las=legtxtlas, tick=FALSE)

  #Reset plot window to normal
  suppressWarnings(graphics::par(opar))
  graphics::layout(1)
}
############################################################################################
############################################################################################

############################################################################################
############################################################################################
#

#'Generate multiple plots of SOM model results
#'
#'Plots SOM output using Parallel Coordinate Plots, Multi-dimensional Scaling (MDS), and Sammons mapping.  kmeans, wards, and PCA
#'@param som_obj a SOM object
#'@param legcex sets legend cex
#'@param labtype determines label type. "IDs" and "XYs" available.
#'@param labsize sets label size
#'@param colormod sets colorscheme
#'@return Panel a-c plot SOM profile results using a range of techniques designed to enhance understanding of results. A list is also returned that includes the following items.
#'  \itemize{
#'  \item {EVALTAB} a dataframe containing classification assignments for SOM, kmeans (KM), Ward (HC). Principal component scores for PC1 and PC2 are also included.
#'  \item {SOM_Profiles} SOM profiles (a.k.a., codebook vectors)
#'  \item {KM_Profiles} K-means profiles (a.k.a., cluster centers)
#'  \item {HC_Profiles} Ward's hierarchical clustering profiles (a.k.a., cluster centers)
#'  \item {PC1_Loadings} loadings from the first principal component
#'  \item {PC2_Loadings} loadings from the second principal component
#'  \item {CLASSIC_AGREEMENT} measures of class agreement ranging from 0:1 where 1 reflects perfect agreement
#' }
#'@details Additional visualization of SOM results can lead to stronger interpretion. Panel (a) presenst a parallel coordinates plot of SOM profile results in order visualize differences in the individual components. Each vertical bar represents a variable axis where minimum values are on the bottom and maximum values are on top. A line is then used to connect component values for each profile. Panel (b) presents a multidimensional scaling plot that allows visualization of how well the SOM profile results span the original data space and how they compare to profiles generated from kmeans and Ward's heirarchical clustering. Panel (c) presents a Sammon's mapping of SOM profiles that effectively visualizes the the SOM using a distance preserving map projection. This provides a clearer picture of similarity/dissimilarity among SOM profiles.
#'@export

sommix_eval<-function(som_obj, labtype="IDs", labsize=1, legcex=1, colormod=NULL){
  #Set up data
  data<-data.frame(som_obj$data)
  vars<-names(data)
  obs<-dim(data)[1]
  p<-dim(data)[2]

  #Get summary information
  som_summ<-sommix_summ(som_obj)

  #########################################################################
  #SOM grid
  somx<-som_summ$som_grid$xdim
  somy<-som_summ$som_grid$ydim
  somk<-som_summ$som_grid$somk
  somgrd<-som_summ$som_grid$topo

  #Set function defaults
  if (is.null(colormod)) colormod <- grDevices::terrain.colors(somk)

  #Set labels
  IDs<-som_summ$som_coords$SOM_ID
  XYs<-som_summ$som_coords$SOM_XY

  #Specify label type
  if ( identical(labtype, "IDs")) {
    labtypes <- IDs
  }  else {
    labtypes<- XYs
  }

  nodelabs<-labtypes

  #Get SOM Profiles
  som.profiles<-som_summ$som_codes
  som.coords<-data.frame(X=som_summ$som_coords$SOM_X, Y=som_summ$som_coords$SOM_Y)
  row.names(som.profiles)<-nodelabs
  som.classif<-som_summ$som_class$SOM_ID

  #Identify number of learning iterations
  niter<-length(som_obj$changes)

  #########################################################################
  #Apply kmeans for comparison
  k.mod<-stats::kmeans(x=data, centers=somk, nstart=10, algorithm="Hartigan-Wong",
                iter.max=niter)
  k.profiles=data.matrix(k.mod$centers)
  k.class=k.mod$cluster

  #Apply Ward's clustering for comparison
  h.mod<-stats::hclust(stats::dist(data), method="ward.D2")
  hc.class<-stats::cutree(h.mod, k=somk)
  hc.profiles <- NULL
  for(k in 1:somk){
    hc.profiles <- rbind(hc.profiles, colMeans(data[hc.class == k, , drop = FALSE]))
  }

  #########################################################################
  #Merge training data with discovered profiles for comparision/evaluation
  data2=rbind(data, som.profiles, k.profiles, hc.profiles)

  #set Data distance matrix for evaluation with MDS
  dist.mat2<-stats::dist(x=data2, method="euclidean")

  #Apply multidimensional scaling to visualize profiles over dataspace
  loc <- stats::cmdscale(dist.mat2)
  mds.x <- loc[, 1]
  mds.y <- -loc[, 2] # reflect so North is at the top

  #########################################################################
  #Compare class assignments
  g1 <- som.classif
  g2 <- k.class
  g3 <- hc.class
  tab1 <- table(g1, g2)
  tab2 <-table(g1,g3)

  somkm<-e1071::classAgreement(tab1)$rand
  somhc<-e1071::classAgreement(tab2)$rand

  classrand<-data.frame(SOM_KM=somkm, SOM_HC=somhc)

  ##########################################################################
  #Apply PCA and summarize scores by SOM Class
  pca.mix<-stats::prcomp(data, center=FALSE, scale=FALSE)
  pca.loadings<-round(pca.mix$rotation,2)

  IDdf<-data.frame(SOM_ID=IDs)


  pc1.sc<-pca.mix$x[,1]
  pc1.som<-stats::aggregate(pc1.sc, by=list(SOM_ID=som.classif), FUN=mean)
  pc1.som.summ<-merge(IDdf, pc1.som, by=c("SOM_ID"), all.x=TRUE)
  pc1.som.summ.t<-pc1.som.summ[,-1]


  pc2.sc<-pca.mix$x[,2]
  pc2.som<-stats::aggregate(pc2.sc, by=list(SOM_ID=som.classif), FUN=mean)
  pc2.som.summ<-pc2.som[,-1]
  pc2.som.summ<-merge(IDdf, pc2.som, by=c("SOM_ID"), all.x=TRUE)
  pc2.som.summ.t<-pc2.som.summ[,-1]

  ##########################################################################

  comptab<-data.frame(SOM_Classif=som.classif, KM_Classif=k.class, HC_Classif=hc.class,
                      PC1_Scores=pc1.sc, PC2_Scores=pc2.sc)
  evalout<-list(COMPTAB=comptab, SOM_Profiles=som.profiles, KM_Profiles=k.profiles, HC_Profiles=hc.profiles,
                PC1_Loadings=pca.loadings[,1], PC2_Loadings=pca.loadings[,2],
                CLASSIF_AGREEMENT=classrand)


  ##########################################################################
  #Generate Evaluation plots
  opar<-graphics::par()
  graphics::par(mfrow=c(1,1), mar=c(5,4,2,2), ask=FALSE, family='serif', pty="m", cex=0.75)
  graphics::layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), heights=c(1,1), widths=c(1,1),respect = FALSE)

  #Generate Parallel Coordinates plot
  MASS::parcoord(som.profiles, col=colormod,
                 lty=1,var.label=FALSE, lwd=2)
  graphics::title(main="a) Parallel Coordinate Plot for SOM Profiles")
  #graphics::box()
  graphics::box(which="figure")
  graphics::legend(x=0.05, y=1.5, legend=nodelabs, lty=1, pch=16,col=colormod, bty="n",
                   cex=legcex, title="SOM", xpd=NA, ncol=ifelse(somk>9,2,1))

  #Generate MDS Plot
  ## note asp = 1, to ensure Euclidean distances are represented correctly
  plot(mds.x, mds.y, type = "p", xlab = "MDS X", ylab = "MDS Y", asp = 1, axes = TRUE,
       main = "b) Multi-dimensional Scaling Plot for Profile Comparison", pch=c(rep(20, obs), rep(16,somk), rep(3,somk), rep(5,somk)),
       col=c(rep("lightgrey",obs), colormod, rep("black", somk), rep("black", somk)),
       cex=c(rep(1,obs), rep(1.75, somk), rep(1.25, somk), rep(1.25, somk)))

  graphics::legend("bottomleft", legend=c("OBS", "SOM", "k-Means", "Wards"), pch=c(20,16,3,5),
         col=c("lightgrey","black","black", "black"), cex=legcex)

  graphics::box(which="figure")

  #Generate a Sammon's Map to look at Map Distortion
  mapdist<-MASS::sammon(d=stats::dist(som.profiles), y=data.matrix(som.coords),
                  k=2, niter=1000, trace=FALSE)
  plot(mapdist$points, pch=20, col=colormod, cex=5, xlab="Sammons X", ylab="Sammons Y",
       main="c) Sammon's Map of SOM Profiles")
  graphics::text(mapdist$points, labels = nodelabs)
  graphics::box(which="figure")

  #Reset plot window to normal
  suppressWarnings(graphics::par(opar))
  graphics::layout(1)

  return(invisible(evalout))

}
##########################################################################

