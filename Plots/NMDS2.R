#NMDS using share_subsampled OTU table to 545 reads

library(ggplot2)
library(vegan)

shared_subsampled_545<-read.table("shared_subsampled_545.txt", header=TRUE, row.names=1)
shared_subsampled_545_percent<-shared_subsampled_545/(rowSums(shared_subsampled_545))
Env2<-read.table("Env2.txt", header=TRUE)

sol1_545<-metaMDS(shared_subsampled_545_percent, distance = "bray", k = 2, trymax = 50)
sol2_545<-metaMDS(shared_subsampled_545_percent, distance = "jaccard", binary=TRUE, k = 2, trymax = 50)

NMDS_bray_545=data.frame(x=sol1_545$point[,1],y=sol1_545$point[,2], Location=as.factor(Env2[,2]), Disinfection=as.factor(Env2[,3]), Country=as.factor(Env2[,4]), Platform=as.factor(Env2[,5]))
NMDS_jaccard_545=data.frame(x=sol2_545$point[,1],y=sol2_545$point[,2], Location=as.factor(Env2[,2]), Disinfection=as.factor(Env2[,3]), Country=as.factor(Env2[,4]), Platform=as.factor(Env2[,5]))

NMDS_bray_545_1<-ggplot(data=NMDS_bray_545,aes(x,y,colour=Location, shape=Disinfection)) + geom_point(size=4) +ggtitle("NMDS - Bray Curtis distances") 
ggsave("NMDS_bray_1.png")

NMDS_jaccard_545_1<-ggplot(data=NMDS_jaccard_545,aes(x,y,colour=Location, shape=Disinfection)) + geom_point(size=4) +ggtitle("PCoA - Jaccard distances") 
ggsave("NMDS_jaccard_1.png")
