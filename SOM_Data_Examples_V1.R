#Basic Example Application of SOM to benchmark machine learning data

#Load packages
install.packages('mlbench')
library(mlbench)
data("Glass")
str(Glass)

ex.data<-Glass[,1:9]
ex.data.class<-Glass$Type
table(ex.data.class)

#Apply variable evaluation function
var.eval(ex.data)

#Apply correlation evaluation function
cor.eval(ex.data)

#Scale data for multivariate analyses
ex.data.sc<-scale(ex.data)

#Apply PCA evaluation function
pca.eval(ex.data.sc)

#Appy group structure evaluation function
grp.eval(ex.data.sc, kmn=2, kmx=10, iter.max=1000)

#Apply map size evaluation tools
map.size.eval(data=ex.data.sc, kmx=20, iter.max=1000)

#Identify optimal seed for SOM
opt.seed<-som.seed(ex.data.sc, somx=3, somy=2)

#Fit example som
set.seed(opt.seed)
ex.som<-som(ex.data.sc, grid=somgrid(xdim=3,ydim=2, "rectangular"), 
             rlen=10000, alpha=c(0.05,0.01), dist.fct="euclidean")

#Summarize Output
som.summ<-som.fit.summ(ex.som)

#Apply custom plotting functions
par(mfrow=c(1,1))
#Plot SOM profiles with error bars
map.bar.plot(ex.som)
#Plot SOM profiles
map.profile.plot(ex.som)
#Plot component mappings
map.comp.plot(ex.som)
#Evaluate SOM fit
map.fit.plots(ex.som)

#Create SOM metric for subsequent analyses
som.exp.metric<-exp.metric(ex.som)

