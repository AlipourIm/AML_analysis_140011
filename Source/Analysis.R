# Import libraries
library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(reshape2)
library(plyr)
library(ggplot2)
library(stringr)
library(gplots)

##############################################################
# set working directory to Downloads folder

setwd("~/Downloads/Bioinformatics_project")

##############################################################
# Now download the dataset and import the data into our code using GEOquery library.

series_accession_number <- "GSE48558"
platform_name <- "GPL6244"

gset <- getGEO(series_accession_number, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "Data/")

if (length(gset) > 1) {
  idx <- grep(platform_name, attr(gset, "names"))
} else {
  idx <- 1
}
gset <- gset[[idx]]

##################################################################
# Now label the normal patients and then label aml patients as test(as mentioned in instructions)
gset <-
  gset[, which(gset$source_name_ch1 == "AML Patient" |
                 gset$`phenotype:ch1` == "Normal")]

func <- function(x) {
  if (gset$source_name_ch1[x] == "AML Patient") {
    return("Test")
  } else {
    spll <- strsplit2(gset$source_name_ch1[x] , "\\+")[1, 1]
    return(paste0("Normal_" , spll))
  }
}

gr <- sapply(1:length(gset$`phenotype:ch1`) , func)

##################################################################
# Now get the expression matrix
expr <- exprs(gset)

##################################################################
# Check if the expression matrix is normalized or not
print(paste0("Maximum Expression: " , max(expr)))

##################################################################
# Since the max expression amount was 13.7..., the expression matrix was
# normalized, if if was a big number like 5000, it wasn't and we should've 
# normalized it, before performing a log operation on it we make sure
# there is no 0 in our expression matrix so that we won't have -infinity
# after performing the log operation. The operations are commented here since
# the data was already normalized, but if we wanted to perform the operations
# these were the operations.

# print(paste0("Minimum Expression: " , min(expr))) # we use this line to ckeck for potential 0 value in our expression matrix
# expr <- log2(1 + expr)  # this ensures there is no 0 value and no -infinity is formed
# exprs(gset)<- expr      # this finishes the process

##################################################################
##################################################################
## Quality control:
##################################################################
##################################################################
# First make a box plot and check if the data is normalized or not.
# Make a pdf file

pdf("Results/boxplot.pdf" , width = 90, height = 40)

boxplot(expr, ylab="Expression rate", xlab="Gene")

# Dump it into created file
dev.off()

# As you can see, the samples are normalized so we won't perform normalization on them

##################################################################
# If the data wasn't normalized however, we would normalize it with following commands

#expr <- normalizeQuantiles(expr) # normalizer the data
#exprs(gset) <- expr  # put it in exprs(gset)

#pdf("Results/boxplot_afterNormalize.pdf" , width = 32) # create a new pdf file
#boxplot(expr)  # make a boxplot for normalized data
#dev.off()      # dump it into the file

##################################################################
# Plot the correlation heat map and dump it into a file
pdf("Results/CorHeatmap.pdf", width=20, height = 20)
pheatmap(cor(expr), labels_row = gr, labels_col = gr, color = bluered(256), border_color = NA)
dev.off()

##################################################################
##################################################################
## Reducing the dimensions:
##################################################################
##################################################################
# Method 1: PCA(Principal Component Analysis)
pca <- prcomp(expr)
pdf("Results/Pca.pdf" , width = 10, height = 10)
plot(pca)
plot(pca$x[, 1:2])  # We can only use PC1 and PC2 according to previous plot
dev.off()

# Scale the data(by rotating it)
expr.scaled <- t(scale(t(expr) , scale = FALSE))

# Now perform PCA on rotated data as we did before and save it in a file
pca <- prcomp(expr.scaled)
pdf("Results/pca_scaled.pdf", width = 15, height = 15)
plot(pca)
plot(pca$x[, 1:2])
dev.off()

# Now perform what we did on samples from dataset
pcar <- data.frame(pca$rotation[, 1:3] , group = gr)
pdf("Results/PCA_Samples.pdf", width = 15 , height = 15)
ggplot(pcar , aes(x = PC1 ,y = PC2 , color = group), size = 4) + geom_point(size = 4) + theme_bw()
dev.off()

##################################################################
##################################################################
## Analyzing the differential expression:
##################################################################
##################################################################
# According to previous plots, we see that "Test" samples are close to "CD34" samples
# therefore we name them with AML_Near_CD34 for easier access and more meaning full name.
gr2 <- gr
gr2[which((pcar$PC2 > 0.11 & pcar$group == "Test"))] <-
  "AML_Near_CD34"

# Now we convert the group gr2 to a factor so that it will no longer be a normal string
gr2 <- factor(gr2)
gset$description <- gr2

# With these lines we make a matrix that shows each sample belongs to which group :))
onehot <- model.matrix(~ description + 0, gset)
colnames(onehot) <- levels(gr2)

colnames(onehot) <-
  c(sapply(colnames(onehot) , function(x)
    str_replace(x, " ", "_")))
# Fit a linear model to the data
fit <- lmFit(gset , onehot)

# Define what you are comparing with what! Here I want to compare AML_Near_CD34 with CD_34 samples
contrast <- makeContrasts(AML_Near_CD34 - Normal_CD34 , levels = onehot)

fit2 <- contrasts.fit(fit , contrast)
fit2 <- eBayes(fit2 , 0.01)
tT <- topTable(fit2, adjust = "fdr", sort.by = "B" , number = Inf)

# Remove excessive and unnecessary columns and data that was gathered from annotation
tT2 <- subset(tT , select = c("Gene.symbol", "Gene.ID" , "Gene.title", "adj.P.Val"  , "logFC"))
# Write the table into a file
write.table(tT2, file = "Results/aml_cd34.txt" , sep = "\t" , row.names = FALSE , quote = FALSE)

##################################################################
# Find genes with maximum expression difference
aml.up <- subset(tT2, logFC > 1 & adj.P.Val < 0.05)

# Find unique values from this set and remove Gene names with "///" and add all of them!
aml.up.genes <-unique(as.character(strsplit2((aml.up$Gene.symbol),"///")))

# Write results in a file(Genes with more difference in expression)
write.table(aml.up.genes, "Results/amd_cd34_up_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Now do the same for genes that have less expression(opposite of what we did before)
aml.down <- subset(tT2, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2((aml.down$Gene.symbol),"///")))

# Once again write the results into a file
write.table(aml.down.genes, "Results/amd_cd34_down_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
