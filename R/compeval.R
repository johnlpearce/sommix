############################################################################################
#sommix: A package to apply Kohonen's Self-Organizing Map algorithm
############################################################################################

#'Provides a range of statistics for evaluating how components are related to a SOM model.
#'
#'Calculates a range of statistics for components (i.e., variables) of the training data in order to better understand relationships with the SOM model.
#'@param som_obj an object returned from sommix_summ
#'@param labsize sets label size on plots
#'@param colormod sets colors for plots
#'@details Understanding how well components relate to the SOM model can provide insights into variable importance. Here correlations, root means square error, and mean absolute error are calculated using observed vs predicted values from the SOM model. Measures of spatial autocorrelation are provided to enhance understanding of spatial patterning among components on the SOM.
#'@return a list containing component evaluation statistics
#' \itemize{
#'  \item {Component_Fit} a dataframe of model fit statistics for each component/variable
#'  \item {Component_Map} a dataframe of Moran's I statistics for each component/variable
#' }
#'@export

compeval<-function(som_obj, labsize=1, colormod=NULL){

  #########################################################################
  #Extract training data
  data<-data.frame(som_obj$data[[1]])
  data$OBS_ID<-rownames(data)
  vars<-names(data.frame(som_obj$codes[[1]]))
  node<-som_obj$unit.classif
  profiles<-data.frame(som_obj$codes[[1]])
  p<-dim(profiles[2])

  #Extract SOM Coordinate information
  summ_obj<-sommix_summ(som_obj)
  som_xy<-data.frame(summ_obj$GRID)
  xdim<-summ_obj$PARAMETERS$xdim
  ydim<-summ_obj$PARAMETERS$ydim
  somk<-xdim*ydim

  #Set function defaults
  if (is.null(colormod)) colormod <- grDevices::terrain.colors(n=p)

  #Extract SOM Profiles
  som_codes<-summ_obj$PROFILES
  som_codes$NODE<-rownames(som_codes)

  #Extract Classification Assignments
  som_classif_df<-summ_obj$CLASSIF
  som_classif<-summ_obj$CLASSIF$NODE

  #########################################################################
  #Generate SOM predictions
  som_pred<-merge(som_classif_df, som_codes, by="NODE", all.x=TRUE)
  #########################################################################
  #Construct Evaluation Data Sets
  eval_data<-merge(data, som_pred, by="OBS_ID", suffixes=c("","_p"))
  #########################################################################

  #########################################################################
  #Variable statistics for SOM-based predictions
  ## Evaluate discriminatory power of SOM-based profiles assignments
  R.p=NULL
  RMSE.p=NULL
  MAE.p=NULL

  #Asses correlation and precision measures

  for (i in 1:length(vars)){
      varn=vars[i]
      Y=eval_data[,varn]
      X=eval_data[,paste(varn, "_p", sep="")]

      #Identify Correlation
      r.p<-stats::cor(x=X, y=Y, method="pearson", use="pairwise.complete.obs")

      #Extract RMSE
      err.p<-Y-X
      rmse.p<-round(sqrt(mean(err.p^2, na.rm=TRUE)),4)
      mae.p<-round(mean(abs(err.p), na.rm=TRUE),4)
      RMSE.p<-c(RMSE.p, rmse.p)
      MAE.p<-c(MAE.p, mae.p)
      R.p<-c(R.p,round(r.p,2))
  }

  comp_eval<-data.frame(VARS=vars,R=as.numeric(R.p),
                       RMSE=RMSE.p, MAE=MAE.p)


  comp_eval$R_P_Rank<-rank(-comp_eval$R)
  comp_eval$RMSE_P_Rank<-rank(comp_eval$RMSE)
  comp_eval$MAE_P_Rank<-rank(comp_eval$MAE)

  colnames(comp_eval)<-c("Variable", "R","RMSE", "MAE",
                          "R_Rank","RMSE_Rank","MAE_Rank")

  #############################################################################
  #Evaluate 'spatial'organization of map via global tests of spatial autocorrelation
  if (ydim>1){
    #Convert grid to spatial coordinates
    som_coords<-sp::SpatialPoints(som_xy[,c("SOM_X", "SOM_Y")])

    #Global Network
    som_knn<-spdep::knn2nb(spdep::knearneigh(som_coords, k=round(somk*.2, 0)))

    #Neighborhood Distances
    som_dist<-unlist(spdep::nbdists(som_knn, coords=som_coords))

    #Assess profiles

    MI=NULL
    MI_P=NULL

    for (j in 1:length(vars)){
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

    sp_cor_tab<-data.frame(vars,MI, MI_P)
    colnames(sp_cor_tab)<-c("Variable", "Morans I", "Pval")
    sp_cor_tab$Rank<-rank(-sp_cor_tab$`Morans I`)

    #print("Map Structure/Autocorrelation Evaluation Complete")
  } else { sp_cor_tab<-NA}

  #############################################################################
  #Evaluation Plots
  opar<-graphics::par(mfrow=c(1,1), mar=c(5,4,3,2), pty="s", cex=1, xpd=FALSE)
  graphics::par(mfrow=c(2,2), mar=c(4,4,1.5,1), family='serif', ask=FALSE, xpd=FALSE)

  #Plots
  #Pearson
  graphics::barplot(height=comp_eval$R, names.arg=comp_eval$Variable,
          las=2, cex.names=labsize, main="a) Predicted Correlations", ylab="R",
          ylim=c(-1,1), col=colormod)
  graphics::box()

  #RMSE
  graphics::barplot(height=comp_eval$RMSE, names.arg=comp_eval$Variable,
                    las=2, cex.names=labsize, main="b) Predicted RMSE", ylab="RMSE",
                    ylim=c(0,max(comp_eval$RMSE)*1.2), col=colormod)
  graphics::box()


  #MAE
  graphics::barplot(height=comp_eval$MAE, names.arg=comp_eval$Variable,
                    las=2, cex.names=labsize, main="c) Predicted MAE", ylab="RMSE",
                    ylim=c(0,max(comp_eval$MAE)*1.2), col=colormod)
  graphics::box()

  #Morans I
  graphics::barplot(height=sp_cor_tab$'Morans I', names.arg=sp_cor_tab$Variable,
                    las=2, cex.names=labsize, main="c) Profile Spatial Autocorrelation", ylab="Morans I",
                    col=colormod)
  graphics::box()

  on.exit(graphics::par(opar))
  on.exit(graphics::layout(1))

  #############################################################################
  #Set list with feature evaluation results
  compeval_summ<-list(Component_Fit=comp_eval, Component_Map=sp_cor_tab)
  return(invisible(compeval_summ))

}

