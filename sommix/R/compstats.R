############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################

#'Summary Statistics for observations assigned to each Self-Organizing Maps (SOM) node (i.e., map unit)
#'
#'Calculates summary statistics for data subsets assigned to each SOM node (i.e., map unit)
#'@param sommix_summ_obj an object returned from sommix_summ
#'@return a list containing summary statistics
#'@export

compstats<-function(sommix_summ_obj){
  #########################################################################
  #Extract training data
  somm_xy<-sommix_summ_obj$som_coords
  data<-data.frame(sommix_summ_obj$trn_data)
  varnames<-names(data[,-c(1:5)])
  node<-data$SOM_XY


  #########################################################################
  #Summarize training data by each SOM node
  trn_data<-data[,-c(1:5)]

  means<-stats::aggregate(x=trn_data[,varnames], by=list(SOM_XY=node), FUN=mean, drop=TRUE)
  means_t<-merge(somm_xy, means, by="SOM_XY", all.x=TRUE)
  means_t<-means_t[order(as.numeric(means_t$SOM_ID)),]
  means_t

  stdev<-stats::aggregate(x=trn_data[,varnames], by=list(SOM_XY=node), FUN=stats::sd, drop=TRUE)
  stdev_t<-merge(somm_xy, stdev, by="SOM_XY", all.x=TRUE)
  stdev_t<-stdev_t[order(as.numeric(stdev_t$SOM_ID)),]
  stdev_t

  mn<-stats::aggregate(x=trn_data[,varnames], by=list(SOM_XY=node), FUN=min, drop=TRUE)
  mn_t<-merge(somm_xy, mn, by="SOM_XY", all.x=TRUE)
  mn_t<-mn_t[order(as.numeric(mn_t$SOM_ID)),]
  mn_t

  mx<-stats::aggregate(x=trn_data[,varnames], by=list(SOM_XY=node), FUN=max, drop=TRUE)
  mx_t<-merge(somm_xy, mx, by="SOM_XY", all.x=TRUE)
  mx_t<-mx_t[order(as.numeric(mx_t$SOM_ID)),]
  mx_t

  med<-stats::aggregate(x=trn_data[,varnames], by=list(SOM_XY=node), FUN=stats::median, drop=TRUE)
  med_t<-merge(somm_xy, med, by="SOM_XY", all.x=TRUE)
  med_t<-med_t[order(as.numeric(med_t$SOM_ID)),]
  med_t

  iqr<-stats::aggregate(x=trn_data[,varnames], by=list(SOM_XY=node), FUN=stats::IQR, drop=TRUE)
  iqr_t<-merge(somm_xy, iqr, by="SOM_XY", all.x=TRUE)
  iqr_t<-iqr_t[order(as.numeric(iqr_t$SOM_ID)),]
  iqr_t

  summ_stats<-list(MEAN=means_t, ST_DEV=stdev_t, MIN=mn_t, MAX=mx_t, MED=med_t, IQRS=iqr_t)
  return(summ_stats)
 }
