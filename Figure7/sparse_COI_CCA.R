library("PMA")
library("kernlab")
library('variancePartition')
library('edgeR')
library('BiocParallel')
library("ade4")

setwd("~/bioinformatics_analysis/PGATES/MT/MT_species_corr/tpm_analysis/")


#######Lactobacillus crispatus
y_data <- read.csv("./input_files/Lactobacillus_crispatus_50cut_TPM.csv",sep=",")

y_data_trim <- y_data[,-1]
rownames(y_data_trim) <- y_data[,1]
y_data_trim <- y_data_trim[ rowSums(y_data_trim)!=0, ]
#y_data_trim <- y_data_trim[rowMeans(y_data_trim!=0)>0.1,]
y_data_trim.scale <- data.frame(scale(t(y_data_trim)))

keeps <- colnames(y_data[-1])
x_data <- read.csv("input_files/MG_taxa_104_cut_CLR.csv",sep=",")
x_data_trim <- x_data[ , (names(x_data) %in% keeps)]
rownames(x_data_trim) <- x_data[,1]
x_data_trim <- x_data_trim[ rowSums(x_data_trim)!=0,]
#x_data_trim <- x_data_trim[rowMeans(x_data_trim!=0)>0.75,]
#x_data_trim <- x_data_trim * 1000000
x_data_trim.scale <- data.frame(scale(t(x_data_trim)))

pca1 = dudi.pca(x_data_trim.scale, scal = FALSE,scann=FALSE)
pca2 = dudi.pca(y_data_trim.scale, scal = FALSE, scann = FALSE)

rv1 = RV.rtest(pca1$tab, pca2$tab, 999)
rv1

cca_output_perm <- CCA.permute(x=t(y_data_trim),z=t(x_data_trim),typex = "standard",typez = "standard",niter=100,trace=TRUE,standardize=TRUE,penaltyxs =c(0.05,0.05,0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.1,0.1,0.1,0.1,.1,0.125,0.125,0.125,0.125,0.125,.125,0.15,0.15,0.15,0.15,0.15,.15,.2,.2,.2,.2,.2,.2),penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))
cca_output_perm
plot(cca_output_perm)
cca_output <- CCA(x=t(y_data_trim),z=t(x_data_trim),typex="standard",typez="standard",penaltyx=0.1,penaltyz=0.2,niter=10000,trace=TRUE,standardize=TRUE,K=3)
cca_output

score.x <- scale(t(y_data_trim)) %*% cca_output$u
score.y <- scale(t(x_data_trim)) %*% cca_output$v
diag(cor(score.x, score.y))

x <- cca_output["u"] 
x <- data.frame(x[1])
row.names(x) <- row.names((y_data_trim))

write.csv(x,"./output_files/Lactobacillus_crispatus_CCA_VOG.csv")

y <- cca_output["v"] 
y <- data.frame(y[1])
row.names(y) <- row.names((x_data_trim))

write.csv(y,"./output_files/Lactobacillus_crispatus_CCA_VOGspecies.csv")

x_data_trim_pred = x_data_trim * y$v.1
x_pred <- colSums(x_data_trim_pred)
y_data_trim_pred = y_data_trim * x$u.1
y_pred <- colSums(y_data_trim_pred)
plot(y_pred,x_pred,pch=16)

predicted_scores <- data.frame(x_pred,y_pred)
write.csv(predicted_scores,"./output_files/Lactobacillus_crispatus_axesCov.csv")


#plotting?
rownames(x_data_trim)[which(cca_output$v[,1]>0)]
combined = cbind(t(y_data_trim[cca_output$u[,1] != 0, ]),
                 t(x_data_trim[cca_output$v[,1] != 0, ]))

pcaRes = dudi.pca(combined, scannf= TRUE,center=TRUE,scale=TRUE)

genotype     = substr(rownames(pcaRes$li), 1, 2)
sampleType  = substr(rownames(pcaRes$l1), 3, 4)
featureType = c(rep("Gene",41),rep("Taxa",2))
sampleInfo  = data.frame(pcaRes$li, genotype, diet=sampleType)
featureInfo = data.frame(pcaRes$c1,
                         feature = substr(colnames(combined), 1, 6))
library(ggplot2)
library(ggrepel)
ggplot() +  geom_point(data = sampleInfo,
                       aes(x = Axis1, y = Axis2), size = 3) +
  geom_label_repel(data = featureInfo,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = featureType),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = featureInfo,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = featureType),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed()+
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pcaRes$eig[1] / sum(pcaRes$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pcaRes$eig[2] / sum(pcaRes$eig), 2)),
       fill = "Feature Type", col = "Sample Type")

pca_component_coord <- pcaRes$co
pca_components <- pcaRes$li

write.csv(pca_components,"Lcrispatus_pcaAxes.csv")
write.csv(pca_component_coord,"Lcrispatus_pcaCoords.csv")


#########Lactobacillus iners
y_data <- read.csv("./input_files/Lactobacillus_iners_50cut_TPM.csv",sep=",")

y_data_trim <- y_data[,-1]
rownames(y_data_trim) <- y_data[,1]
y_data_trim <- y_data_trim[ rowSums(y_data_trim)!=0, ]
#y_data_trim <- y_data_trim[rowMeans(y_data_trim!=0)>0.1,]
y_data_trim.scale <- data.frame(scale(t(y_data_trim)))

keeps <- colnames(y_data[-1])
x_data <- read.csv("./input_files/MG_taxa_104_cut_CLR.csv",sep=",")
x_data_trim <- x_data[ , (names(x_data) %in% keeps)]
rownames(x_data_trim) <- x_data[,1]
x_data_trim <- x_data_trim[ rowSums(x_data_trim)!=0,]
#x_data_trim <- x_data_trim[rowMeans(x_data_trim!=0)>0.75,]
#x_data_trim <- x_data_trim * 1000000
x_data_trim.scale <- data.frame(scale(t(x_data_trim)))

pca1 = dudi.pca(x_data_trim.scale, scal = FALSE,scann=FALSE)
pca2 = dudi.pca(y_data_trim.scale, scal = FALSE, scann = FALSE)

rv1 = RV.rtest(pca1$tab, pca2$tab, 999)
rv1

cca_output_perm <- CCA.permute(x=t(y_data_trim),z=t(x_data_trim),typex = "standard",typez = "standard",niter=100,trace=TRUE,standardize=TRUE,penaltyxs =c(0.05,0.05,0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.1,0.1,0.1,0.1,.1,0.125,0.125,0.125,0.125,0.125,.125,0.15,0.15,0.15,0.15,0.15,.15,.2,.2,.2,.2,.2,.2),penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))
cca_output_perm
plot(cca_output_perm)
cca_output <- CCA(x=t(y_data_trim),z=t(x_data_trim),typex="standard",typez="standard",penaltyx=0.2,penaltyz=0.4,niter=10000,trace=TRUE,standardize=TRUE,K=3)
cca_output

score.x <- scale(t(y_data_trim)) %*% cca_output$u
score.y <- scale(t(x_data_trim)) %*% cca_output$v
diag(cor(score.x, score.y))

x <- cca_output["u"] 
x <- data.frame(x[1])
row.names(x) <- row.names((y_data_trim))

write.csv(x,"./output_files/Lactobacillus_iners_CCA_VOG.csv")

y <- cca_output["v"] 
y <- data.frame(y[1])
row.names(y) <- row.names((x_data_trim))

write.csv(y,"./output_files/Lactobacillus_iners_CCA_VOGspecies.csv")

x_data_trim_pred = x_data_trim * y$v.1
x_pred <- colSums(x_data_trim_pred)
y_data_trim_pred = y_data_trim * x$u.1
y_pred <- colSums(y_data_trim_pred)
plot(y_pred,x_pred,pch=16)

#plotting?
rownames(x_data_trim)[which(cca_output$v[,1]>0)]
combined = cbind(t(y_data_trim[cca_output$u[,1] != 0, ]),
                 t(x_data_trim[cca_output$v[,1] != 0, ]))

pcaRes = dudi.pca(combined, scannf= TRUE,nf=2)

genotype     = substr(rownames(pcaRes$li), 1, 2)
sampleType  = substr(rownames(pcaRes$l1), 3, 4)
featureType = c(rep("Gene",111),rep("Taxa",10))
sampleInfo  = data.frame(pcaRes$li, genotype, diet=sampleType)
featureInfo = data.frame(pcaRes$c1,
                         feature = substr(colnames(combined), 1, 6))
library(ggplot2)
library(ggrepel)
ggplot() +  geom_point(data = sampleInfo,
                       aes(x = Axis1, y = Axis2), size = 3) +
  geom_label_repel(data = featureInfo,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = featureType),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = featureInfo,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = featureType),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed()+
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pcaRes$eig[1] / sum(pcaRes$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pcaRes$eig[2] / sum(pcaRes$eig), 2)),
       fill = "Feature Type", col = "Sample Type")

x_data_trim_pred = x_data_trim * y$v.1
x_pred <- colSums(x_data_trim_pred)
y_data_trim_pred = y_data_trim * x$u.1
y_pred <- colSums(y_data_trim_pred)

pca_component_coord <- pcaRes$co
pca_components <- pcaRes$li

write.csv(pca_components,"./output_files/Liners_pcaAxes.csv")
write.csv(pca_component_coord,"./output_files/Liners_pcaCoords.csv")



######Gardnerella vaginalis
y_data <- read.csv("./input_files/Gardnerella_vaginalis_50cut_TPM.csv",sep=",")

y_data_trim <- y_data[,-1]
rownames(y_data_trim) <- y_data[,1]
y_data_trim <- y_data_trim[ rowSums(y_data_trim)!=0, ]
#y_data_trim <- y_data_trim[rowMeans(y_data_trim!=0)>0.1,]
y_data_trim.scale <- data.frame(scale(t(y_data_trim)))

keeps <- colnames(y_data[-1])
x_data <- read.csv("./input_files/MG_taxa_104_cut_CLR.csv",sep=",")
x_data_trim <- x_data[ , (names(x_data) %in% keeps)]
rownames(x_data_trim) <- x_data[,1]
x_data_trim <- x_data_trim[ rowSums(x_data_trim)!=0,]
#x_data_trim <- x_data_trim[rowMeans(x_data_trim!=0)>0.75,]
#x_data_trim <- x_data_trim * 1000000
x_data_trim.scale <- data.frame(scale(t(x_data_trim)))

pca1 = dudi.pca(x_data_trim.scale, scal = FALSE,scann=FALSE)
pca2 = dudi.pca(y_data_trim.scale, scal = FALSE, scann = FALSE)

rv1 = RV.rtest(pca1$tab, pca2$tab, 999)
rv1

cca_output_perm <- CCA.permute(x=t(y_data_trim),z=t(x_data_trim),typex = "standard",typez = "standard",niter=100,trace=TRUE,standardize=TRUE,penaltyxs =c(0.05,0.05,0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.1,0.1,0.1,0.1,.1,0.125,0.125,0.125,0.125,0.125,.125,0.15,0.15,0.15,0.15,0.15,.15,.2,.2,.2,.2,.2,.2),penaltyzs = c(0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.2,0.3,0.4,0.5,0.6,0.7))
cca_output_perm
plot(cca_output_perm)
cca_output <- CCA(x=t(y_data_trim),z=t(x_data_trim),typex="standard",typez="standard",penaltyx=0.2,penaltyz=0.35,niter=10000,trace=TRUE,standardize=TRUE,K=3)
cca_output

score.x <- scale(t(y_data_trim)) %*% cca_output$u
score.y <- scale(t(x_data_trim)) %*% cca_output$v
diag(cor(score.x, score.y))

x <- cca_output["u"] 
x <- data.frame(x[1])
row.names(x) <- row.names((y_data_trim))

write.csv(x,"./output_files/Gardnerella_vaginalis_CCA_VOG.csv")

y <- cca_output["v"] 
y <- data.frame(y[1])
row.names(y) <- row.names((x_data_trim))

write.csv(y,"./output_files/Gardnerella_vaginalis_CCA_VOGspecies.csv")

x_data_trim_pred = x_data_trim * y$v.1
x_pred <- colSums(x_data_trim_pred)
y_data_trim_pred = y_data_trim * x$u.1
y_pred <- colSums(y_data_trim_pred)
plot(y_pred,x_pred,pch=16)


#plotting?
rownames(x_data_trim)[which(cca_output$v[,1]>0)]
combined = cbind(t(y_data_trim[cca_output$u[,1] != 0, ]),
                 t(x_data_trim[cca_output$v[,1] != 0, ]))

pcaRes = dudi.pca(combined, scannf= TRUE)

genotype     = substr(rownames(pcaRes$li), 1, 2)
sampleType  = substr(rownames(pcaRes$l1), 3, 4)
featureType = c(rep("Gene",185),rep("Taxa",9))
sampleInfo  = data.frame(pcaRes$li, genotype, diet=sampleType)
featureInfo = data.frame(pcaRes$c1,
                         feature = substr(colnames(combined), 1, 6))
library(ggplot2)
library(ggrepel)
ggplot() +  geom_point(data = sampleInfo,
                       aes(x = Axis1, y = Axis2), size = 3) +
  geom_label_repel(data = featureInfo,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = featureType),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = featureInfo,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = featureType),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed()+
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pcaRes$eig[1] / sum(pcaRes$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pcaRes$eig[2] / sum(pcaRes$eig), 2)),
       fill = "Feature Type", col = "Sample Type")


x_data_trim_pred = x_data_trim * y$v.1
x_pred <- colSums(x_data_trim_pred)
y_data_trim_pred = y_data_trim * x$u.1
y_pred <- colSums(y_data_trim_pred)

pca_component_coord <- pcaRes$co
pca_components <- pcaRes$li

write.csv(pca_components,"./output_files/Gvag_pcaAxes.csv")
write.csv(pca_component_coord,"./output_files/Gvag_pcaCoords.csv")



