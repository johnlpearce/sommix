############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################

#'Kohonen's Self-Organizing Maps (SOM) algrithm
#'
#'Applies Kohonens self-organizing map via kohonen using robust tuning parameters and optimized initial values
#'@param X is a data object to train the map. Should be a numerical or factor based data matrix. data.frame objects are not supported.
#'@param somx is x-dimension of the SOM
#'@param somy is y-dimension of the SOM
#'@param maptopo specifies the topology of the SOM grid as either "rectangular" or "hexagonal"
#'@param distmet specifies dissimilarity metric. Current options include "sumofsquares", "euclidean", "manhattan", and "tanimoto". Default is to use "Euclidean" for continuous data, and "tanimoto" for factors.
#'@param lmode specifies the learning algorithm. The default is "online" but "batch" and "pbatch" are available via kohonen
#'@param itermax specifies the number if learning iterations
#'@param seedopt specifies if optimal initialization values are to be used. "Y" or "N" accepted.
#'@param nstarts specifies number of initialization schemes to test if inits are not provided
#'@return a 'som' object with class Kohonen
#'@description Kohonen's Self-Organizing Map algorithm is an artificial neural network (ANN) that applied competitive learning with an integraged neighborhood function in order to discover and visualize patterns in multivariate datasets. The 'map' is a low-dimensional representation that illustrates the discovered patterns (i.e., profiles) in a compact, spatially organized way. This unique feature allows for larger numbers of classes to be more easily understand, offer users the ability to construct high resolution classifications. This is a key distinction from traditional clustering/dimension reduction approaches which often seek to minimal representation.
#'@references Kohonen, T. (1995) Self-Organizing Maps. Springer-Verlag
#'@export
#'@examples
#'#NIEHS Mixtures Workshop dataset1
#'data(dataset1)
#'#Apply SOM
#'sommod<-sommix(scale(dataset1[,3:9]), somx=3, somy=2)
#'#Summarize Output
#'somsumm<-sommix_summ(sommod)
#'#Plot Results
#'sommix_bar(som_obj=sommod, varnames=NULL, colormod=NULL,nodelab=TRUE, labtype="IDS",
#'labsize=1, addFreq=TRUE, freqtype="frq", freqsize=1, legsymsize=2, leglabsize=1,
#'legtxtlas=2, barstat="MED")
#'#Evaluate SOM
#'sommix_eval(sommod, labtype = "IDs")


sommix <- function (X, somx=3, somy=2, maptopo=NULL, itermax=NULL, seedopt="Y",
                    distmet=NULL, nstarts=NULL, lmode=NULL)

{
  #Set k
  k<-somx*somy

  #Specify training data
  data.trn<-data.matrix(X)

  #Set function defaults
  if (is.null(maptopo)) maptopo <- "rectangular"
  if (is.null(itermax)) itermax <- 500*(somx*somy)
  if (is.null(distmet)) distmet<-"euclidean"
  if (is.null(lmode)) lmode <- "online"
  if (is.null(nstarts)) nstarts <- 5

  nstart<-nstarts

  #Identify intitial values using multiple starts
  set.seed(1)
  ran.vals=sample(1:100,nstarts)
  ran.vals

  #Set eval object
  QE<-NULL

  for(i in 1:nstarts){
    set.seed(ran.vals[i])
    init.som<-kohonen::som(X=data.trn, grid=kohonen::somgrid(somx,somy,maptopo),
                rlen=itermax, mode=lmode, alpha=c(0.05, 0.01),
                dist.fcts=distmet)
    qe<-mean(init.som$distances)
    QE<-c(QE,qe)
  }

  seedvals<-ran.vals[which(QE == min(QE))]
  seedval<-seedvals[1]

  #Fit SOM with optimal seed value
  if ( identical(seedopt, "Y")) {
    set.seed(seedval)
    }  else {
    set.seed(NULL)
   }

  #Apply SOM via Kohonens C code
  sommix_obj<-kohonen::som(data.trn,
                           grid=kohonen::somgrid(xdim=somx, ydim=somy, topo=maptopo),
                           rlen=itermax, mode=lmode, alpha=c(0.05, 0.01),
                           dist.fcts=distmet)

  #Add number of random initializations
  sommix_obj$grid$nstart<-nstart

  return(sommix_obj)

}
