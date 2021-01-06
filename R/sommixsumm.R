############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################
#Summarizes select output from Kohonen's Self-Organizing Maps (SOM) algorithm
#'SOM Summary Information
#'
#'sommix_summ returns summary information for Kohonens self-organizing map algorithm
#'@param som_obj is a som object with class 'kohonen'
#'@return 'sommix_summ'    a list of SOM results
#'\itemize{
  #'  \item {PARAMETERS} map specifics
  #'  \item {GRID}  map nodes and coordinates
  #'  \item {CODES} node profiles (a.k.a., cluster centers, classes, codebook vectors)
  #'  \item {FREQ} node assignment frequencies
  #'  \item {CLASSIF} a vector identifying the best matching unit (BMU) on the SOM grid
  #'  \item {DATA} input data merged with BMU assignments
  #'  }
#'@export

sommix_summ <- function (som_obj) {
  #Extract data info
  data<-data.frame(som_obj$data[[1]])
  vars<-names(data.frame(som_obj$codes[[1]]))
  nobs<-dim(data)[1]
  nvars<-length(vars)

  #########################################################################
  #SOM Parameters
  xdim<-som_obj$grid$xdim
  ydim<-som_obj$grid$ydim
  somk<-xdim*ydim
  som_parameters<-data.frame(xdim=xdim, ydim=ydim, niter=length(som_obj$changes),
                       topo=som_obj$grid$topo, neighborhood=som_obj$grid$neighbourhood.fct,
                       distance=som_obj$dist.fcts)

  #########################################################################
  #SOM GRID and coordinates
  XYs<-paste(round(som_obj$grid$pts[,1],1),round(som_obj$grid$pts[,2],1), sep="")
  NODE<-1:somk
  #SOM Coordinates and reference table
  som_grid<-data.frame(NODE=factor(NODE),
                         SOM_X=as.numeric(som_obj$grid$pts[,1]),
                         SOM_Y=as.numeric(som_obj$grid$pts[,2]))

  #########################################################################
  #Extract SOM Codes
  som_codes<-data.frame(som_obj$codes)
  row.names(som_codes)<-NODE

  #########################################################################
  #SOM Model
  fit<-somFitStats(som_obj)
  som_fit<-data.frame(N=nobs, VARS=nvars, k=somk, R2=fit$R2, ADJ_R2=fit$ADJ_R2,
                      MAE=fit$DISTANCE_MAE, RMSE=fit$DISTANCE_RMSE,
                      TOTSS=fit$TOTSS, TOTWCSS=fit$TOTWCSS, BCSS=fit$BCSS,
                      AIC=fit$AIC)

  #########################################################################
  #Summarize Frequency Assignments
  N<-table(som_obj$unit.classif)
  FREQ<-round(100*(table(som_obj$unit.classif)/sum(table(som_obj$unit.classif))),2)
  freq.tab<-data.frame(N, FREQ)
  freq.tab2<-merge(som_grid, freq.tab, by.x="NODE", by.y="Var1", all.x=TRUE)
  freq.tab2$Var1.1<-NULL
  colnames(freq.tab2)<-c("NODE","SOM_X","SOM_Y","N","FREQ")
  freq.tab2<-freq.tab2[order(as.numeric(freq.tab2$NODE)),]

  #########################################################################
  #Extract Classification Assignments
  classif<-data.frame(OBS=1:dim(data)[1],OBS_ID=rownames(data),
                        NODE=factor(som_obj$unit.classif), ERROR=som_obj$distances)

  classif2<-merge(classif, som_grid, by="NODE")
  classif2<-classif2[,c("OBS", "OBS_ID", "NODE","SOM_X","SOM_Y", "ERROR")]
  classif3<-classif2[do.call(order,classif2),]
  rownames(classif3)<-classif3$OBS

  ########################################################################
  #Set Training Data with SOM assignments
  dataclassif<-cbind(classif3, som_obj$data[[1]])

  ########################################################################
  #Specify output
  summ<-list(PARAMETERS=som_parameters, GRID=som_grid, FIT=som_fit, PROFILES=som_codes, FREQ=freq.tab2,
             CLASSIF=classif3, DATA=dataclassif)
  return(summ)

  }
