#NMDS plots using share_subsampled OTU table subsampled to 11587 reads.
#Input data: shared_subsampled.txt, Env.txt

library(ggplot2)
library(vegan)

shared_subsampled<-read.table("shared_subsampled.txt", header=TRUE, row.names=1)
shared_subsampled_percent<-shared_subsampled/(rowSums(shared_subsampled))
Env<-read.table("Env.txt", header=TRUE)

sol1<-metaMDS(shared_subsampled_percent, distance = "bray", k = 2, trymax = 50)
sol2<-metaMDS(shared_subsampled_percent, distance = "jaccard", binary=TRUE, k = 2, trymax = 50)

NMDS_bray=data.frame(x=sol1$point[,1],y=sol1$point[,2], Location=as.factor(Env[,2]), Disinfection=as.factor(Env[,3]), Country=as.factor(Env[,4]), Platform=as.factor(Env[,5]))
NMDS_jaccard=data.frame(x=sol2$point[,1],y=sol2$point[,2], Location=as.factor(Env[,2]), Disinfection=as.factor(Env[,3]), Country=as.factor(Env[,4]), Platform=as.factor(Env[,5]))
                                                                                 
NMDS_bray_1<-ggplot(data=NMDS_bray,aes(x,y,colour=Location, shape=Disinfection)) + geom_point(size=4) +ggtitle("NMDS - Bray Curtis distances") 
ggsave("NMDS_bray_1.png")

NMDS_jaccard_1<-ggplot(data=NMDS_jaccard,aes(x,y,colour=Location, shape=Disinfection)) + geom_point(size=4) +ggtitle("PCoA - Jaccard distances") 
ggsave("NMDS_jaccard_1.png")
