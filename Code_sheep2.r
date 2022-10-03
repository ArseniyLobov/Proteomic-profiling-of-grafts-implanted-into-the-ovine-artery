###This is a code to reproduce statistical analysis of proteomics data from Proteomic profiling of tissue-engineered, small-diameter vascular grafts implanted into the ovine carotid artery

##Openinig_the_data
#We reccomend to set the working directory to make easy to reproduce the code
#setwd("your directory")

library(readxl)
dat <- data.frame(read_excel("dat_sheep_filt.xlsx"))

head(dat)
str(dat)
dat1 <- dat
rownames(dat) = make.names(dat[,1], unique=TRUE)
dat <- dat[,-1]
head(dat)
dat[dat == 0] <- NA

#Opening the sample info
fact <- data.frame(read_excel("sheep_info.xlsx"))

rownames(fact) <- fact[,1]

fact$Groups <- as.factor(fact$Groups)
fact$Groups
## A - Native carotid arteries; G - Biodegradable vascular grafts


## Qualititative analysis
library(VennDiagram)
library(RColorBrewer)

Art <- dat[which(rowMeans(!is.na(dat[,c(1:12)])) >= 0.75), ]
Gra <- dat[which(rowMeans(!is.na(dat[,c(13:24)])) >= 0.75), ]

venn.diagram(
  x = list(rownames(Art), rownames(Gra)),
  category.names = c("Art" , "Gra"),
  filename = '#venn_diagramm.png',
  resolution = 600,
  output=F
)

library(gplots)
v.table1 <- venn(list(rownames(Art), rownames(Gra)))
print(v.table1)


## Quantitative analysis
#Removing rows with a lot of missing values
dat2 <- dat[which(rowMeans(!is.na(dat)) >= 22/24), ]
mean(complete.cases(dat2))
colSums(is.na(dat2))

#Raw data
library(RColorBrewer)

#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$Groups]
boxplot(dat2, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$Groups), fill = pal, bty = "n", xpd = T)

#knn imputation of missng values
library(impute)
dat_knn <- impute.knn(t(dat2), k = 5)
dat_knn <- t(dat_knn$data)
mean(complete.cases(dat_knn))

#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
boxplot(dat_knn, outline = FALSE, col = cols, main = "Data after missed values impuration")
legend("topright", levels(fact$Groups), fill = pal, bty = "n", xpd = T)

#dev.off()
colSums(dat_knn)
head(dat_knn)

dat_log <- log2(dat_knn+1)
boxplot(dat_log, outline = FALSE, col = cols, main = "Data after missed values impuration and log-transf.")
legend("topright", levels(fact$Groups), fill = pal, bty = "n", xpd = T)

#Quantile normalization
library(limma)
data_norm <- normalizeQuantiles(dat_log)
boxplot(data_norm, outline = FALSE, col = cols, main = "Normalized data")
legend("topright", levels(fact$Groups), fill = pal, bty = "n", xpd = T)



#MA-plot
maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  scatter.smooth(x = X, y = Y,
                 main = main, pch = pch,
                 xlab = xlab, ylab = ylab,
                 lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

maplot(data_norm[, rownames(subset(fact,Groups=="A"))], data_norm[, rownames(subset(fact,Groups=="G"))], main = "Normalized data")


library(mixOmics)
dat_pca <- pca(t(data_norm), ncomp = 5, center = TRUE)

#tiff('PCA_group.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact$Groups, legend = TRUE, ellipse = T,
          title = 'PCA')

#dev.off()


#Limma
library(limma)

X <- model.matrix(~ fact$Groups)
X

fit <- lmFit(data_norm, design = X, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit <- eBayes(fit)

# Dif_expr_table
topTable(efit, coef = 2)
full_list <- topTable(efit, coef = 2, number = length(data_norm[,2]))
#write.csv(full_list,'Dif_expr_A_vs_G.csv')
head(full_list)


#Vulcano
library(EnhancedVolcano)

#tiff('Dif_expr.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(full_list,
                lab = rownames(full_list),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-5.5, 4.5),
                ylim = c(0, 28),
                title ="A vs G",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


p_aboveS <- full_list$adj.P.Val <= 0.05
sum(p_aboveS)

head(full_list)
#Verification of top-proteins
data_norm1 <- as.matrix(data_norm)
boxplot(as.matrix(dat_knn)[c("TAGLN"),] ~ Groups, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(as.matrix(dat_knn)[c("ITGB1"),] ~ Groups, data = fact,
  varwidth = TRUE, log = "y", las = 1)
boxplot(as.matrix(dat_knn)[c("FLNA"),] ~ Groups, data = fact,
  varwidth = TRUE, log = "y", las = 1)

boxplot(as.matrix(dat_knn)[c("CAPG"),] ~ Groups, data = fact,
  varwidth = TRUE, log = "y", las = 1)
boxplot(as.matrix(dat_knn)[c("CTSD"),] ~ Groups, data = fact,
  varwidth = TRUE, log = "y", las = 1)
boxplot(as.matrix(dat_knn)[c("VIM.1"),] ~ Groups, data = fact,
  varwidth = TRUE, log = "y", las = 1)
