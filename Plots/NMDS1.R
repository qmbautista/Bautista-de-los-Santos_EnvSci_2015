#NMDS of community compostion using phylotype table (output of mothur) subsampled to 545 reads.
#Using Bray Curtis and Jaccard distances
#From:http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

library(ggplot2)
library(vegan)

shared_subsampled<-read.table("shared_subs_545.txt", header=TRUE, row.names=1)
shared_subsampled_percent<-shared_subsampled/(rowSums(shared_subsampled))
Env<-read.table("Env.txt", header=TRUE)

sol1<-metaMDS(shared_subsampled_percent, distance = "bray", k = 2, trymax = 50)
sol2<-metaMDS(shared_subsampled_percent, distance = "jaccard", binary=TRUE, k = 2, trymax = 50)

NMDS_bray=data.frame(x=sol1$point[,1],y=sol1$point[,2], Location=as.factor(Env[,2]), Disinfection=as.factor(Env[,3]), Country=as.factor(Env[,4]), Platform=as.factor(Env[,5]))
NMDS_jaccard=data.frame(x=sol2$point[,1],y=sol2$point[,2], Location=as.factor(Env[,2]), Disinfection=as.factor(Env[,3]), Country=as.factor(Env[,4]), Platform=as.factor(Env[,5]))

######################Ellipses for Bray Curtis plot ###############################
plot.new()
ord1<-ordiellipse(sol1, as.factor(Env[,3]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell1 <- data.frame()
for(g in levels(NMDS_bray$Disinfection)){
  if(g!="" && (g %in% names(ord1))){
  df_ell1 <- rbind(df_ell1, cbind(as.data.frame(with(NMDS_bray[NMDS_bray$Disinfection==g,],
  veganCovEllipse(ord1[[g]]$cov,ord1[[g]]$center,ord1[[g]]$scale)))
  ,Disinfection=g))
  }
}

NMDS1.mean=aggregate(NMDS_bray[,1:2],list(group=NMDS_bray$Disinfection),mean)

shape_values<-seq(1,11)

#######################################################################################

######################Ellipses for Jaccard plot ######################################

plot.new()
ord2<-ordiellipse(sol2, as.factor(Env[,3]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell2 <- data.frame()
for(g in levels(NMDS_jaccard$Disinfection)){
  if(g!="" && (g %in% names(ord2))){
  df_ell2 <- rbind(df_ell2, cbind(as.data.frame(with(NMDS_jaccard[NMDS_jaccard$Disinfection==g,],
  veganCovEllipse(ord2[[g]]$cov,ord2[[g]]$center,ord2[[g]]$scale)))
 ,Disinfection=g))
  }
}

NMDS2.mean=aggregate(NMDS_jaccard[,1:2],list(group=NMDS_jaccard$Disinfection),mean)

shape_values<-seq(1,11)


#######################################################################################

NMDS_bray_1<-ggplot(data=NMDS_bray,aes(x,y, colour=Disinfection)) 
NMDS_bray_1<-NMDS_bray_1+geom_path(data=df_ell1, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
NMDS_bray_1<-NMDS_bray_1+geom_point(aes(shape=Country), size=4)+theme_bw()+ylab("NMDS2")+xlab("NMDS1")  
pdf("NMDS_bray_ell.pdf")
print(NMDS_bray_1)
dev.off()

NMDS_jaccard_1<-ggplot(data=NMDS_jaccard,aes(x,y, colour=Disinfection)) 
NMDS_jaccard_1<-NMDS_jaccard_1+geom_path(data=df_ell2, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
NMDS_jaccard_1<-NMDS_jaccard_1+geom_point(aes(shape=Country), size=4)+theme_bw()+ylab("NMDS2")+xlab("NMDS1")
pdf("NMDS_jaccard_ell.pdf")
print(NMDS_jaccard_1)
dev.off()
