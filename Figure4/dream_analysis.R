library('variancePartition')
library('edgeR')
library('BiocParallel')

setwd("~/bioinformatics_analysis/PGATES/MT/MT_core_diff/")

data_test <- read.csv("MT_UMD_KEGGsums.txt",sep="\t",row.names = 1)
filter = rowSums(cpm(data_test)>0.1) >= 5

geneExpr = DGEList( data_test[filter,] )
geneExpr = calcNormFactors( geneExpr )

metadata <- read.csv("MT_core_key.csv",sep=",",row.names = 1)

param = SnowParam(4, "SOCK", progressbar=TRUE)
 
register(param)


# The variable to be tested must be a fixed effect
form <- ~ 0 + CST_MT_red + (1|Subject) 

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(geneExpr, form, metadata)

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata )
L1 = getContrast( vobjDream, form, metadata, c("CST_MT_redI", "CST_MT_redIV"))
L2 = getContrast( vobjDream, form, metadata, c("CST_MT_redI", "CST_MT_redIII"))
L3 = getContrast( vobjDream, form, metadata, c("CST_MT_redIII", "CST_MT_redIV"))
L = cbind(L1,L2,L3)

fit = dream( vobjDream, form, metadata,L,ddf="Kenward-Roger")

# Get results of hypothesis test on coefficients of interest
cst_I_IV <- topTable(fit, coef='L1',number = 6000,sort.by = "logFC")
cst_I_III <- topTable(fit, coef='L2',number = 6000,sort.by = "logFC")
cst_III_IV <- topTable(fit, coef='L3',number = 6000,sort.by = "logFC")


#pathway enrichment
drops <- c("AveExpr","t","P.Value","z.std")
CST_I_IV_red <- cst_I_IV[ , !(names(cst_I_IV) %in% drops)]
CST_I_IV_red$GeneID <- rownames(CST_I_IV_red)
CST_I_IV_red <- CST_I_IV_red[,c(3,1,2)]

form2 = ~ (1|Subject) + (1|CST_MT_red)
vp = fitExtractVarPartModel(vobjDream, form2, metadata)

plotVarPart(sortCols(vp))

write.csv(cst_I_IV,"MT_cstI_IV_kegg_de.csv",sep=",")
write.csv(cst_I_III,"MT_cstI_III_kegg_de.csv",sep=",")
write.csv(cst_III_IV,"MT_cstIII_IV_kegg_de.csv",sep=",")
write.csv(vp,"MT_kegg_variance_partition.csv")

####MG data analysis
data_test <- read.csv("MG_UMD_KEGGsums.txt",sep="\t",row.names = 1)
filter = rowSums(cpm(data_test)>0.1) >= 5

geneExpr = DGEList( data_test[filter,] )
geneExpr = calcNormFactors( geneExpr )

metadata <- read.csv("MT_core_key.csv",sep=",",row.names = 1)

param = SnowParam(4, "SOCK", progressbar=TRUE)

register(param)


# The variable to be tested must be a fixed effect
form <- ~ 0 + CST_MG_red + (1|Subject) 

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(geneExpr, form, metadata)

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata )
L1 = getContrast( vobjDream, form, metadata, c("CST_MG_redI", "CST_MG_redIV"))
L2 = getContrast( vobjDream, form, metadata, c("CST_MG_redI", "CST_MG_redIII"))
L3 = getContrast( vobjDream, form, metadata, c("CST_MG_redIII", "CST_MG_redIV"))
L = cbind(L1,L2,L3)

fit = dream( vobjDream, form, metadata,L)

# Get results of hypothesis test on coefficients of interest
cst_I_IV <- topTable(fit, coef='L1',number = 6000,sort.by = "logFC")
cst_I_III <- topTable(fit, coef='L2',number = 6000,sort.by = "logFC")
cst_III_IV <- topTable(fit, coef='L3',number = 6000,sort.by = "logFC")


#pathway enrichment
drops <- c("AveExpr","t","P.Value","z.std")
CST_I_IV_red <- cst_I_IV[ , !(names(cst_I_IV) %in% drops)]
CST_I_IV_red$GeneID <- rownames(CST_I_IV_red)
CST_I_IV_red <- CST_I_IV_red[,c(3,1,2)]

form2 = ~ (1|Subject) + (1|CST_MT_red)
vp = fitExtractVarPartModel(vobjDream, form2, metadata)

plotVarPart( sortCols(vp))

write.csv(cst_I_IV,"MG_cstI_IV_kegg_de.csv",sep=",")
write.csv(cst_I_III,"MG_cstI_III_kegg_de.csv",sep=",")
write.csv(cst_III_IV,"MG_cstIII_IV_kegg_de.csv",sep=",")
write.csv(vp,"MG_kegg_variance_partition.csv")