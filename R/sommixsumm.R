############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################
#Summarizes select output from Kohonen's Self-Organizing Maps (SOM) algorithm
#'SOM Summary Information
#'
#'sommix_summ returns summary information for Kohonens self-organizing map algorithm
#'@param som_obj is a som object with class 'kohonen'
#'@return 'sommix'    a list of SOM results
#'\itemize{
  #'  \item {som_grid} map specifics
  #'  \item {som_coords} map reference coordinates system
  #'  \item {som_codes} SOM profiles (a.k.a., cluster centers, classes, codebook vectors)
  #'  \item {som_freq} SOM class frequencies summaries
  #'  \item {som_class} a vector of SOM class assignments
  #'  \item {trn_data} Training data merged with class assignments
  #'  }
#'@export

sommix_summ <- function (som_obj) {
  #Set up data
  data<-data.frame(som_obj$data)
  vars<-names(data)
  obs<-dim(data)[1]
  p<-dim(data)[2]

  #########################################################################
  #SOM grid
  xdim<-som_obj$grid$xdim
  ydim<-som_obj$grid$ydim
  somk<-xdim*ydim
  som_grid<-data.frame(xdim=xdim, ydim=ydim, somk=somk,
                       topo=som_obj$grid$topo, neighborhood=som_obj$grid$neighbourhood.fct)

  #########################################################################
  #Extract SOM dimensions and coordinates
  XYs<-paste(round(som_obj$grid$pts[,1],2),round(som_obj$grid$pts[,2],2), sep="")
  SOM_IDs<-1:somk
  #SOM Coordinates and reference table
  som_coords<-data.frame(SOM_ID=factor(SOM_IDs), SOM_XY=factor(XYs),
                         SOM_X=as.numeric(som_obj$grid$pts[,1]),
                         SOM_Y=as.numeric(som_obj$grid$pts[,2]))

  #########################################################################
  #Extract SOM Codes
  som_codes<-data.frame(som_obj$codes)
  row.names(som_codes)<-SOM_IDs

  #########################################################################
  #Summarize Frequency Assignments
  som.n<-table(som_obj$unit.classif)
  som.freq<-round(100*(table(som_obj$unit.classif)/sum(table(som_obj$unit.classif))),2)
  freq.tab<-data.frame(som.n, som.freq)
  freq.tab2<-merge(som_coords, freq.tab, by.x="SOM_ID", by.y="Var1", all.x=TRUE)
  freq.tab2$Var1.1<-NULL
  colnames(freq.tab2)<-c("SOM_ID","XY","X","Y","N","FREQ")
  freq.tab2<-freq.tab2[order(as.numeric(freq.tab2$SOM_ID)),]

  #########################################################################
  #Extract Classification Assignments
  class.tab<-data.frame(OBS=1:length(som_obj$unit.classif), SOM_ID=factor(som_obj$unit.classif), ERROR=som_obj$distances)

  class.tab2<-merge(class.tab, som_coords, by="SOM_ID")
  class.tab2<-class.tab2[,c("OBS", "SOM_ID","SOM_XY","SOM_X","SOM_Y", "ERROR")]
  class.tab3<-class.tab2[do.call(order,class.tab2),]
  rownames(class.tab3)<-class.tab3$OBS

  ########################################################################
  #Set Training Data with SOM assignments
  trn.data<-data.frame(OBS=1:length(som_obj$unit.classif), SOM_ID=as.factor(som_obj$unit.classif),
                       SOM_XY=as.factor(class.tab3$SOM_XY), SOM_X=class.tab3$SOM_X, SOM_Y=class.tab3$SOM_Y,
                       som_obj$data)

  ########################################################################
  #Specify output
  summ<-list(som_grid=som_grid, som_coords=som_coords, som_codes=som_codes, som_freq=freq.tab2, som_class=class.tab3, trn_data=trn.data)
  return(summ)

  }
