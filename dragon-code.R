###Github repository: https://github.com/RicardoQi/dragon-phylogeny.git

library(ape)
library(reshape2)
library(ggplot2)

#nex file imported
DragonNexus<-read.nexus.data("C:\\Users\\qi199\\Documents\\dragon-phylogeny\\input\\DragonMatrix.nex")

#trait coding achieved according to the table
toothless <- "100100111111000000000001001100011110000110000110000110000110011100000100011100"
mushu <- "100100001111000100100110101100110000000001101100001110000000111100010101100001"
smaug <- "100111110011000110000010010000110001001000001010001000100001111000001011111001"

#convert to data.frame and create matrix
DragonNexusDF<-data.frame(matrix(unlist(DragonNexus), ncol=78,byrow=T))
row.names(DragonNexusDF)<-names(DragonNexus)
DragonDist<-dist(DragonNexusDF,method='binary')
DragonDistMat<-as.matrix(DragonDist)

#introduce weight to rationalize the matrix
WeightsDat<-read.csv("C:\\Users\\qi199\\Documents\\dragon-phylogeny\\input\\Weights.csv")
Weights<-paste0(WeightsDat$Weight,collapse="")
Weights<-strsplit(Weights,split="")[[1]]
WeightsNum<-rep(NA,length(Weights))
for(i in 1:length(WeightsNum)){
  if(Weights[i] %in% LETTERS){
    WeightsNum[i]<-which(LETTERS==Weights[i])+9
  } else {
    WeightsNum[i]<-Weights[i]
  }
}
WeightsNum<-as.numeric(WeightsNum)

#couple weight with trait coding
WtDragonNexus<-DragonNexus
for (i in 1:length(DragonNexus)){
  RepWeight<-DragonNexus[[i]]==1
  WtDragonNexus[[i]][RepWeight]<-WeightsNum[RepWeight]
  RepWeight<-NA
}

#plot created based on weight and distance matrix
WtDragonNexusDF<-data.frame(matrix(unlist(WtDragonNexus),ncol=78,byrow=T))
row.names(WtDragonNexusDF)<-names(WtDragonNexus)
WtDragonDist<-dist(WtDragonNexusDF,method='euclidean')
WtDragonDistMat<-as.matrix(WtDragonDist)
WtPDat<-melt(WtDragonDistMat)
ggplot(data = WtPDat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+scale_fill_gradientn(colours=c("white","blue","green","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Phylogenic trees created 
library(ggtree)
WtDragonTree<-fastme.bal(WtDragonDist)
WtDragonTreeNJ<-nj(WtDragonDist)
ggtree(WtDragonTree,layout="circular")

#provide categeories to the tree and give distinct colours 
Country<-gsub("[0-9\\.]+([^X]+)X*","\\1",WtDragonTree$tip.label)
CountryGroups<-split(WtDragonTree$tip.label, Country)
WtDTcol<-groupOTU(WtDragonTree,CountryGroups)
ggtree(WtDTcol,layout="circular",aes(colour=group))+geom_tiplab(size=2,aes(angle=angle))
