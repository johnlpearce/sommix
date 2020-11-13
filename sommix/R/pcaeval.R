############################################################################################
#pcaeval: Performs principal component analysis
############################################################################################
#'Principal Component Analysis with Diagnostic Plots
#'
#'Performs principal component analysis on the data and returns the results as an object of class prcomp. Scree, variance, and loading plots are provided.
#'
#'@param x a data
#'@param labsize sets label size
#'@param varlabs specifies x-axis variable labels
#'@param colormod can be used to customize plot color schemes. The default uses terrain palletes from grDevices
#'@param ... arguments passed to or from other methods
#'@details
#'Assessing patterns of intercorrelations between x-variables is an important aspect in multivariate analysis. This assessment tool can assist with understanding the complexitiy of the underlying data by identifiying the important primary modes of variance that exist. In most The evaluation tool is based upon application of prcomp for stats. The
#'@return a list of pca results with class "prcomp"
#'@export
#'@examples
#'#NIEHS Mixtures Workshop dataset1
#'data(dataset1)
#'pcaeval(scale(dataset1[,2:9]))


pcaeval<-function(x, labsize=1, colormod=NULL, varlabs=NULL, ...){
  #Set up data
  data<-data.frame(x)
  vars<-names(data)
  obs<-dim(data)[1]
  p<-dim(data)[2]

  #Set function defaults
  if (is.null(colormod)) colormod <- grDevices::terrain.colors(n=p)
  if (is.null(varlabs)) varlabs <- ifelse(nchar(vars)<7,vars,substr(vars, 0,7))

  #Conduct PCA
  pca.mix<-stats::prcomp(data, center=FALSE, scale.=FALSE)

  pca.loadings<-pca.mix$rotation

  #Explore the variance explained by each component
  pr.var=pca.mix$sdev^2
  #Calculate the proportion of variance explained
  pve=round(pr.var/sum(pr.var),2)

  #Evaluation Plots
  opar<-graphics::par()
  graphics::par(mar=c(4,4,2,1))
  graphics::layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), heights=c(2,2), widths=c(2,2),respect = TRUE)

  stats::biplot(pca.mix, xlabs=rep(".", dim(pca.mix$x)[1]), cex=c(2,1), col=c("darkgrey", "darkred"),
                main="a) PCA Biplot")


  #Plot the component Standard deviations to see which are useful
  plot(pca.mix$sdev^2, pch=19, ylab="Eigenvalue", xlab="Principal Component",
       cex=1, type="b", main="b) PCA Scree Plot", col="darkgrey")
  graphics::abline(h=1, col="black", lty=2)


  plot(pve, xlab="Principal Component", ylab="Proportion of Variance", ylim=c(0,1), type='b',
       pch=17, col="darkgrey", main="c) Variance Explained", cex=1)
  graphics::points(cumsum(pve), type='b',pch=19, col="darkgrey")
  graphics::legend("bottomleft", pch=c(17,19), col=c("darkgrey", "darkgrey"), legend=c("Individual", "Cumulative"))

  #Reset plot window to normal
  suppressWarnings(graphics::par(opar))
  graphics::layout(1)

  return(invisible(pca.mix))

}
