#Example Application of Self-organizing Map (SOM) for multidimensional exposure characterization
#Date: 31AUG2020
#Author: John Pearce

#Load Necessary Packages using SOMmix Source File

#Examle data is frm NIEHS Mixtures Workshop
#https://www.niehs.nih.gov/about/events/pastmtg/2015/statistical/index.cfm

#Read in example data
#setwd() #This needs to be assigned to folder with dataset2.
#dir()
ex.data<-read.csv("dataset2.csv", header=TRUE, sep=",")
colnames(ex.data)<-toupper(names(ex.data))
rownames(ex.data)<-ex.data$OBS
summary(ex.data)
head(ex.data)

#Set up exposure variables as training data
exp.vars<-names(ex.data[,3:16])
data.trn=na.omit(ex.data[,exp.vars])
names(data.trn)

#Standardize data for multivariate analyses
data.trn.sc<-scale(data.trn, center=TRUE, scale=TRUE)
head(data.trn.sc)

################################################################################################################################
#Step 1 Explore individual variables
var.eval(data.trn)

###################################################################################################################################
#Step 2 Assess correlational structure/covariation of data
cor.eval(data.trn)

###############################################################################################################################
#Step 3 Examine primary modes of variance using PCA
pca.eval(data.trn.sc)

################################################################################################################################
#Step 4 Examine grouping structure using common strategies applied in cluster analysis
grp.eval(data.trn.sc, kmn=2, kmx=20, iter.max=1500)

###############################################################################################################################
#Step 5 Explore properties of muliple Self Organizing Map (SOM) map dimensions
map.size.eval(data=data.trn.sc, kmx=20, nstart=10, iter.max=1500, grid.topo="rectangular")

###########################################################################

################################################################################################################################
#Develop SOM (e.g., data2)
#Set up map size
somx=4
somy=3

################################################################################################################################
#Set up for Final SOM
#Identify optimal seed for SOM
opt.seed<-som.seed(data.trn.sc, somx=4, somy=3, iter.max=5000)

#Fit example som
set.seed(opt.seed)
ex.som<-som(data.trn.sc, grid=somgrid(xdim=4,ydim=3, "rectangular"),
            rlen=10000, alpha=c(0.05,0.01), dist.fct="euclidean")

#Summarize Output
som.summ<-som.fit.summ(ex.som)
names(som.summ)
som.summ$COORDINATES
som.summ$FREQUENCIES

#Visualize SOM Profiles
par(mfrow=c(1,1))
#Mean Centered Barplot of SOM profiles with error bars
map.bar.plot(ex.som, lab.cex=1.5, label.loc="bottomright", legend.cex=2)
#Radial Secment Plot SOM profiles
map.profile.plot(ex.som)
#Plot component mappings
map.comp.plot(ex.som)
#Evaluate SOM fit
map.fit.plots(ex.som)

###############################################################################
#Create SOM metric for subsequent analyses
som.exp.metric<-exp.metric(ex.som)
head(som.exp.metric)
################################################################################################################################



