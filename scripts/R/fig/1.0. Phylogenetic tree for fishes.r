#############################################################################
#Script to create a phylogenetic tree from Open Tree of Life from time calibrated phylogenies.
#Author: Felix P Leiva
# Date: 20191125
# Modifications: YYYYMMDD (?)
#############################################################################
rm(list=ls()) #clear your work environment
today<-format(Sys.Date(),"%Y%m%d") #setting the date
#############################################################################
#From Windows
# setwd("C:/Users/Invunche/Dropbox/Aafke+Valerio") #Felix Lab
getwd()#to check
#############################################################################
#Libraries
library(car)
library(ape) # Paradis et al 2004
library(nlme) # regression modelling
library(phytools)# Revell et al 2012
library(geiger) #Harmon et al 2007
library(caper)
library(MASS)
library(stringr)
#############################################################################
#load data
data<-read.csv("example.csv")
fish<-subset(data,data$taxa=="fish")
fish<-fish[c(4:11)]
str(fish)

#############################################################################
# load phylogenetic tree
fish.tre<-read.tree("betancurt-R et al 2017.tre")#reading phylogenetic tree (Phylogenetic classification of bony fishes)
str(fish.tre)#Phylogenetic information for more than 11000 fish species!!!
tips<-fish.tre$tip.label#species labels contained in the tree
fish.tre$Nnode#Node numbers
is.binary.tree(fish.tre) # we want this to be TRUE, if FALSE, run line below
is.ultrametric(fish.tre) # if not run line below
fish.tre<-compute.brlen(fish.tre,method = "Grafen")
is.ultrametric(fish.tre) # check again
full.tree<-fish.tre

#############################################################################
#checking list of species in the tree
fish$species2 <- str_replace_all(fish$species," " , '_')
rownames(fish)<-fish$species2 

# Do all the species in the data are in the tree & vice versa?.
nombres.full<-name.check(full.tree,fish)
nombres.full

full.tree2<-drop.tip(full.tree,nombres.full$tree_not_data)
full.tree2$tip.label# 156 species in the data base with phylogenetic information

names.full2<-name.check(full.tree2,fish)
# But still there are 37 species (196-37=193) contained in the dataframe but NO in the tree
exclude.sp<-as.data.frame(names.full2$data_not_tree)

fish2<-fish[ !(fish$species2 %in% exclude.sp$`names.full2$data_not_tree`), ]
#Now are all species in the tree contained in the dataframe and viceversa

# Lets see how look the phylogeny
plotTree(full.tree2,type="fan",fsize=0.5,lwd=1,ftype="i")

#with the data of mass (in the branches) and CTmax
x3<-fish2[c(1,6,7,9)]  
x4<-setNames(x3[,3],rownames(x3))#MASS
y4<-setNames(x3[,2],rownames(x3))#CTMAX
obj2<-contMap(full.tree2,x4,plot=FALSE)
obj2<-setMap(obj2,invert=TRUE)

#pdf(file="for Affke and Valerio.pdf",width=10,height=10,useDingbats = FALSE)
plotTree.wBars(obj2$tre,y4,method="plotSimmap",
               tip.labels=FALSE,fsize=0.7,colors=obj2$cols,type="fan",scale=0.005)
add.color.bar(1.0,obj2$cols,title="Log (mass)",lims=obj2$lims,prompt=FALSE,
              x=0.9*par()$usr[1],y=0.9*par()$usr[3])
#dev.off()
#############################################################################
#Brownian mode of evolution
fit1<-gls(ctmax~body_mass*time_ctmax+abs_latitude,correlation=corPagel(value=1, phy=full.tree2,fixed = TRUE),
          method = "ML",data=fish2)
summary(fit1)
Anova(fit1)

fit1<-gls(ctmax~body_mass,correlation=corPagel(value=0, phy=full.tree2,fixed = TRUE),
          method = "ML",data=fish2)
summary(fit1)
summary(lm(ctmax~body_mass,data=fish2))
Anova(fit1)


fit0.5<-gls(ctmax~body_mass*time_ctmax+abs_latitude,correlation=corPagel(value=0.5, phy=full.tree2,fixed = TRUE),
          method = "ML",data=fish2)
summary(fit0.5)
Anova(fit0.5)

fit0<-gls(ctmax~body_mass*time_ctmax+abs_latitude,correlation=corPagel(value=0, phy=full.tree2,fixed = TRUE),
          method = "ML",data=fish2)
summary(fit0)
Anova(fit0)


anova(fit1,fit0.5,fit0)


# Fitting by usig linear model
fit.lm<-gls(ctmax~body_mass*time_ctmax+abs_latitude,
          method = "ML",data=fish2)
summary(fit.lm)
Anova(fit.lm)

anova(fit0,fit.lm)

####They are same!!!!!
