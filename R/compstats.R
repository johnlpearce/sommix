############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################

#'Summary Statistics for observations assigned to each Self-Organizing Maps (SOM) node (i.e., map unit)
#'
#'Calculates summary statistics for data subsets assigned to each SOM node (i.e., map unit)
#'@param som_obj an object returned from sommix
#'@return a list containing summary statistics
#'@export

compstats<-function(som_obj){
  #########################################################################
  #Extract data from som object
  data<-data.frame(som_obj$data[[1]])
  vars<-names(data.frame(som_obj$codes[[1]]))
  node<-som_obj$unit.classif

  summ<-sommix_summ(som_obj)
  somm_xy<-summ$GRID

  #########################################################################
  #Summarize training data by each SOM node
  trn_data<-data

  means<-stats::aggregate(x=trn_data, by=list(NODE=node), FUN=mean, drop=TRUE)
  means_t<-merge(somm_xy, means, by="NODE", all.x=TRUE)
  means_t<-means_t[order(as.numeric(means_t$NODE)),]
  means_t

  stdev<-stats::aggregate(x=trn_data, by=list(NODE=node), FUN=stats::sd, drop=TRUE)
  stdev_t<-merge(somm_xy, stdev, by="NODE", all.x=TRUE)
  stdev_t<-stdev_t[order(as.numeric(stdev_t$NODE)),]
  stdev_t

  mn<-stats::aggregate(x=trn_data, by=list(NODE=node), FUN=min, drop=TRUE)
  mn_t<-merge(somm_xy, mn, by="NODE", all.x=TRUE)
  mn_t<-mn_t[order(as.numeric(mn_t$NODE)),]
  mn_t

  mx<-stats::aggregate(x=trn_data, by=list(NODE=node), FUN=max, drop=TRUE)
  mx_t<-merge(somm_xy, mx, by="NODE", all.x=TRUE)
  mx_t<-mx_t[order(as.numeric(mx_t$NODE)),]
  mx_t

  med<-stats::aggregate(x=trn_data, by=list(NODE=node), FUN=stats::median, drop=TRUE)
  med_t<-merge(somm_xy, med, by="NODE", all.x=TRUE)
  med_t<-med_t[order(as.numeric(med_t$NODE)),]
  med_t

  iqr<-stats::aggregate(x=trn_data, by=list(NODE=node), FUN=stats::IQR, drop=TRUE)
  iqr_t<-merge(somm_xy, iqr, by="NODE", all.x=TRUE)
  iqr_t<-iqr_t[order(as.numeric(iqr_t$NODE)),]
  iqr_t

  summ_stats<-list(MEAN=means_t, ST_DEV=stdev_t, MIN=mn_t, MAX=mx_t, MED=med_t, IQRS=iqr_t)
  return(summ_stats)
 }
