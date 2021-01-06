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
  xdim<-som_summ$PARAMETERS$xdim
  ydim<-som_summ$PARAMETERS$ydim

  profiles<-data.frame(som_summ$PROFILES)
  p<-dim(profiles)[2]
  vars<-names(profiles)

  #Set Profile labels
  XYs<-paste("[", round(som_summ$GRID$SOM_X), ",", round(som_summ$GRID$SOM_Y),"]", sep="")
  IDs<-paste("[",som_summ$GRID$NODE,"]", sep="")

  #Set Frequency labels
  FREQ<-round(som_summ$FREQ$FREQ,1)
  N<-som_summ$FREQ$N

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
  opar<-graphics::par(mfrow=c(1,1), mar=c(5,4,3,2), pty="s", cex=1)
  graphics::par(mfrow=c(1,1),mar=c(.1,.1,.1,.1), pty="m", xpd=FALSE)
  graphics::layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(5,1), widths=c(1,1),respect = FALSE)

  graphics::stars(x=som_summ$PROFILES, locations=som_summ$GRID[,c("SOM_X", "SOM_Y")],
        draw.segments=TRUE, axes=FALSE, scale=TRUE,
        len = 0.4, col.segments=colormod,
        labels=NULL)

  #Specify which borders to use
  if(identical(som_obj$grid$topo, "rectangular")){
  graphics::symbols(som_summ$GRID[,3:4],
          squares = rep(1, nrow(som_summ$GRID)),
          inches = FALSE, add = TRUE,
          fg = "black", bg = NA)
  } else {graphics::symbols(som_summ$GRID[,c("SOM_X", "SOM_Y")],
                            circles = rep(.5, nrow(som_summ$GRID)),
                            inches = FALSE, add = TRUE,
                            fg = "black", bg = NA)}

  graphics::text(x=som_summ$GRID$SOM_X, y=som_summ$GRID$SOM_Y+.425,
       labels=nodelabs, font=2, cex=labsize)

  graphics::text(x=som_summ$GRID$SOM_X+.375, y=som_summ$GRID$SOM_Y-.4,
       labels=freqlab,
       font=2, cex=freqsize)

  #Plot Legend on lower panel
  graphics::par(mar=c(rep(1,4)))
  plot(x=1:p, y=rep(.9,p), pch=22, cex=legsymsize, col="black", bg=colormod,
       ylim=c(0,1), xlim=c(1, (p +.25)), axes=FALSE, xlab="", ylab="")
  graphics::box(which="outer")
  graphics::axis(side=1, at=1:p, pos=0.85, labels=varnames, cex.axis=leglabsize, las=legtxtlas, tick=FALSE)

  #Reset plot window to normal
  on.exit(graphics::par(opar))
  on.exit(graphics::layout(1))

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
                     nodelab=TRUE, labtype="IDs", labsize=1,
                     addFreq=FALSE, freqtype="frq", freqsize=1,
                     legsymsize=2, leglabsize=1, legtxtlas=2,
                     barstat="MED")
{
  #Apply sommix summary function
  som_summ<-sommix_summ(som_obj)

  #Extract info for plotting
  xdim<-som_summ$PARAMETERS$xdim
  ydim<-som_summ$PARAMETERS$ydim

  profiles<-data.frame(som_summ$PROFILES)
  p<-dim(profiles)[2]
  vars<-names(profiles)

  #Set Profile labels
  XYs<-paste("[", round(som_summ$GRID$SOM_X,1), ",", round(som_summ$GRID$SOM_Y,1),"]", sep="")
  IDs<-paste("[",som_summ$GRID$NODE,"]", sep="")


  #Set Frequency labels
  FREQ<-round(som_summ$FREQ$FREQ,1)
  N<-som_summ$FREQ$N

  #Generate summary statistics for barplots
  summtab<-compstats(som_obj)

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
  opar<-graphics::par(mfrow=c(1,1), mar=c(5,4,3,2), pty="s", cex=1, xpd=FALSE)
  graphics::par(mar=c(.5, 2, 1, 1), pty="m", cex=1)

  #Generate a layout matrix
  M<-matrix(data=1:(xdim*ydim),
            nrow=ydim, ncol=xdim, byrow=TRUE)

  #Flip matrix to match SOM grid
  M1=apply(t(M),1,rev)

  M2<-rbind(M1, rep(1+(xdim*ydim),xdim))

  mf <- graphics::layout(M2, widths = rep(1, ncol(M2)),
               heights = c(rep(1, nrow(M2)-1), 1), respect = FALSE)


  #Run loop to create plots for each SOM node
  for(i in 1:(xdim*ydim)){

    plot.vals<-as.numeric(plot.bars[i,])
    plot.err<-as.numeric(plot.error[i,])

    bp<-graphics::barplot(height=plot.vals, col=colormod, axes=TRUE, axisnames=FALSE,
                ylim=c(min(plot.bars)+ -1.5, 1.2*max(plot.bars)))
    suppressWarnings(error.bar(bp,plot.vals,plot.err))
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
  on.exit(graphics::par(opar))
  on.exit(graphics::layout(1))
}
############################################################################################
############################################################################################

############################################################################################
############################################################################################
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
  xdim<-som_summ$PARAMETERS$xdim
  ydim<-som_summ$PARAMETERS$ydim
  somk<-xdim*ydim
  somgrd<-som_summ$PARAMETERS$topo

  #Set function defaults
  if (is.null(colormod)) colormod <- grDevices::rainbow(somk)

  #Set labels
  IDs<-som_summ$GRID$NODE
  XYs<-paste(som_summ$GRID$SOM_X,som_summ$GRID$SOM_X, sep=",")

  #Specify label type
  if ( identical(labtype, "IDs")) {
    labtypes <- IDs
  }  else {
    labtypes<- XYs
  }

  nodelabs<-labtypes

  #Get SOM Profiles
  som.profiles<-som_summ$PROFILES
  som.coords<-data.frame(SOM_X=som_summ$GRID$SOM_X, SOM_Y=som_summ$GRID$SOM_Y)
  rownames(som.profiles)<-nodelabs
  rownames(som.coords)<-nodelabs
  som.classif<-som_summ$CLASSIF$NODE

  #Identify number of learning iterations and random initializations
  niter<-length(som_obj$changes)
  nstarts<-som_obj$grid$nstart

  #########################################################################
  #Apply kmeans for comparison
  k.mod<-stats::kmeans(x=data, centers=somk, nstart=nstarts, algorithm="Hartigan-Wong",
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

  IDdf<-data.frame(NODE=IDs)

  #PC1
  pc1.sc<-pca.mix$x[,1]
  pc1.som<-stats::aggregate(pc1.sc, by=list(NODE=som.classif), FUN=mean)
  pc1.som.summ<-merge(IDdf, pc1.som, by=c("NODE"), all.x=TRUE)
  pc1.som.summ.t<-pc1.som.summ[,-1]

  #PC2
  pc2.sc<-pca.mix$x[,2]
  pc2.som<-stats::aggregate(pc2.sc, by=list(NODE=som.classif), FUN=mean)
  pc2.som.summ<-pc2.som[,-1]
  pc2.som.summ<-merge(IDdf, pc2.som, by=c("NODE"), all.x=TRUE)
  pc2.som.summ.t<-pc2.som.summ[,-1]

  ##########################################################################

  comptab<-data.frame(SOM_Classif=som.classif, KM_Classif=k.class, HC_Classif=hc.class,
                      PC1_Scores=pc1.sc, PC2_Scores=pc2.sc)
  evalout<-list(COMPTAB=comptab, SOM_Profiles=som.profiles, KM_Profiles=k.profiles, HC_Profiles=hc.profiles,
                PC1_Loadings=pca.loadings[,1], PC2_Loadings=pca.loadings[,2],
                CLASSIF_AGREEMENT=classrand)


  ##########################################################################
  #Generate Evaluation plots
  opar<-graphics::par(mfrow=c(1,1), mar=c(5,4,3,2), pty="s", cex=1, xpd=FALSE)
  graphics::par(mfrow=c(1,1), mar=c(5,4,2,2), ask=FALSE, family='serif', pty="m", cex=0.75, xpd=TRUE)
  graphics::layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), heights=c(1,1), widths=c(1,1),respect = FALSE)

  #Generate Parallel Coordinates plot
  MASS::parcoord(som.profiles, col=colormod,
                 lty=1,var.label=FALSE, lwd=2)
  graphics::title(main="a) Parallel Coordinate Plot")
  #graphics::box()
  graphics::box(which="figure")
  graphics::legend("topleft", legend=nodelabs, lty=1, col=colormod, bty="n",
                   cex=legcex, xpd=NA, ncol=ifelse(somk>9,2,1),
                   inset=c(-.05,0))

  #Generate MDS Plot
  ## note asp = 1, to ensure Euclidean distances are represented correctly
  plot(mds.x, mds.y, type = "p", xlab = "MDS X", ylab = "MDS Y", asp = 1, axes = TRUE,
       main = "b) Multi-dimensional Scaling Plot", pch=c(rep(20, obs), rep(16,somk), rep(3,somk), rep(5,somk)),
       col=c(rep("lightgrey",obs), colormod, rep("black", somk), rep("black", somk)),
       cex=c(rep(1,obs), rep(2, somk), rep(1.5, somk), rep(1.5, somk)))

  graphics::legend("bottomleft", legend=c("SOM", "k-Means", "Wards"), pch=c(16,3,5),
         col=c("#FF0000","black", "black"), cex=legcex, bty="n")

  graphics::box(which="figure")

  #Generate a Sammon's Map to look at Map Distortion
  mapdist<-MASS::sammon(d=stats::dist(som.profiles), y=data.matrix(som.coords),
                  k=2, niter=1000, trace=FALSE)
  plot(mapdist$points, pch=20, col=colormod, cex=5, xlab="Sammons X", ylab="Sammons Y",
       main="c) Sammon's Map")
  graphics::text(mapdist$points, labels = nodelabs)
  graphics::box(which="figure")

  #Generate U matrix to show class distinction with cluster boundaries
  K<-round(somk/3)
  som.hc <- stats::cutree(stats::hclust(kohonen::object.distances(som_obj, "codes")), K)

  plot(som_obj, type="dist.neighbours", main = "d) Between-Class Distance (U-Matrix)", palette.name = grDevices::terrain.colors)
  graphics::text(x=som_obj$grid$pts[,1], y=som_obj$grid$pts[,2],
                 labels=nodelabs, font=2, cex=1)
  #kohonen::add.cluster.boundaries(som_obj, som.hc)
  #graphics::title(sub=paste("Cluster boundary k =", K))
  graphics::box(which="figure")

  #Generate plot illustrating class cohesion
  plot(som_obj, type="quality", palette.name = grDevices::terrain.colors, main="e) Within-Class Distance")
  graphics::text(x=som_obj$grid$pts[,1], y=som_obj$grid$pts[,2],
                 labels=nodelabs, font=2, cex=1)
  #kohonen::add.cluster.boundaries(som_obj, som.hc)
  #graphics::title(sub=paste("Cluster boundary k =", K))
  graphics::box(which="figure")

  #Reset plot window to normal
  on.exit(graphics::par(opar))
  on.exit(graphics::layout(1))

  return(invisible(evalout))

}
##########################################################################

