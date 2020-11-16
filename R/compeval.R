############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################

#'Provides a range of statistics for evaluating how components are related to a SOM model.
#'
#'Calculates a range of statistics for components (i.e., variables) of the training data in order to better understand relationships with the SOM model.
#'@param sommix_summ_obj an object returned from sommix_summ
#'@details Understanding how well components relate to the SOM model can provide insights into variable importance. Here correlations, root means square error, and mean absolute error are calculated using observed vs predicted values from the SOM model. Measures of spatial autocorrelation are provided to enhance understanding of spatial patterning among components on the SOM.
#'@return a list containing component evaluation statistics
#' \itemize{
#'  \item {Component_Fit} a dataframe of model fit statistics for each component/variable
#'  \item {Component_Map} a dataframe of Moran's I statistics for each component/variable
#' }
#'@export

compeval<-function(sommix_summ_obj){

  #########################################################################
  #Extract training data
  data_trn<-data.frame(sommix_summ_obj$trn_data)
  varnames<-names(data_trn[,-c(1:5)])

  #Extract SOM Coordinate information
  som_xy<-data.frame(sommix_summ_obj$som_coords)
  somk<-sommix_summ_obj$som_grid$somk
  somx<-sommix_summ_obj$som_grid$xdim
  somy<-sommix_summ_obj$som_grid$ydim

  #Extract SOM Codes
  som_codes<-sommix_summ_obj$som_codes
  som_codes$SOM_ID<-rownames(som_codes)

  #Extract Classification Assignments
  som_classif_df<-sommix_summ_obj$som_class
  som_classif<-som_classif_df$SOM_ID

  #########################################################################
  #Generate SOM predictions
  som_pred<-merge(som_classif_df, som_codes, by="SOM_ID", all.x=TRUE)

  #########################################################################
  #Construct Evaluation Data Sets
  eval_data<-merge(data_trn, som_pred, by="OBS", suffixes=c("",".p"))

  #########################################################################

  #########################################################################
  #Variable statistics for SOM-based predictions
  ## Evaluate discriminatory power of SOM-based profiles assignments
  R.p=NULL
  RMSE.p=NULL
  MAE.p=NULL

  #Asses correlation and precision measures
  for (i in 1:length(varnames)){
      var.n=varnames[i]
      Y=eval_data[,var.n]
      X=eval_data[,paste(var.n, ".p", sep="")]

      #Identify Correlation
      r.p<-stats::cor(Y,X, method="pearson")

      #Extract RMSE
      err.p<-Y-X
      rmse.p<-round(sqrt(mean(err.p^2)),4)
      mae.p<-round(mean(abs(err.p)),4)
      RMSE.p<-c(RMSE.p, rmse.p)
      MAE.p<-c(MAE.p, mae.p)
      R.p<-c(R.p,round(r.p,2))
  }

  comp_eval<-data.frame(VARS=varnames,R=as.numeric(R.p),
                       RMSE=RMSE.p, MAE=MAE.p)


  comp_eval$R_P_Rank<-rank(-comp_eval$R)
  comp_eval$RMSE_P_Rank<-rank(comp_eval$RMSE)
  comp_eval$MAE_P_Rank<-rank(comp_eval$MAE)

  colnames(comp_eval)<-c("Variable", "R","RMSE", "MAE",
                          "R_Rank","RMSE_Rank","MAE_Rank")

  #############################################################################
  #Evaluate 'spatial'organization of map via global tests of spatial autocorrelation
  if (somy>1){
    #Convert grid to spatial coordinates
    som_coords<-sp::SpatialPoints(som_xy[,c("SOM_X", "SOM_Y")])

    #Global Network
    som_knn<-spdep::knn2nb(spdep::knearneigh(som_coords, k=round(somk*.2, 0)))

    #Neighborhood Distances
    som_dist<-unlist(spdep::nbdists(som_knn, coords=som_coords))

    #Set profiles
    profiles<-sommix_summ_obj$som_codes

    MI=NULL
    MI_P=NULL

    for (j in 1:length(varnames)){
      if (stats::sd(profiles[,j]) > 0) {
        mt<-spdep::moran.test(x=profiles[,j], spdep::nb2listw(som_knn), na.action=stats::na.exclude)
        mi<-as.numeric(mt$statistic)
        mi.p<-as.numeric(mt$p.value)

      } else {
        mi=NA
        mi.p=NA
      }

      MI<-c(MI, mi)
      MI_P<-c(MI_P, mi.p)

    }

    sp_cor_tab<-data.frame(varnames,MI, MI_P)
    colnames(sp_cor_tab)<-c("Variable", "Morans I", "Pval")
    sp_cor_tab$Rank<-rank(-sp_cor_tab$`Morans I`)

    #print("Map Structure/Autocorrelation Evaluation Complete")
  } else { sp_cor_tab<-NA}


  #############################################################################


  compeval_summ<-list(Component_Fit=comp_eval, Component_Map=sp_cor_tab)
  return(compeval_summ)

}

