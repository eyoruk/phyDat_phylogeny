
from clustalw .aln file download

http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi 		file convertion
output format a PHYLIP /sequential 

select and then copy your sequence to .phy file (you can use mediate file as .txt)

--r--

.txt can be read or convert from .aln to .txt 

                                                                  ## PCA_for_MUltiple_Alignment ##
> setwd('C:/Users/username/Desktop/yeni_PCA_Aygün/')
> dir()
> library(ape)
> library(phangorn)
> library(factoextra)
> ort = read.phyDat('31.phy')
> dm = dist.ml(ort)
> View(dm)
> xlsx::write.xlsx(as.matrix(dm), 'distance_for_31.xlsx')
> pca = prcomp(dm, scale = TRUE)
> xlsx::write.xlsx(pca$x, 'pcafor31.xlsx')
> str(pca)
> fviz_eig(pca)
> ct = c(rep("A",5),rep("B",16))
> fviz_pca_ind(pca,
+              col.ind = "cos2", # Color by the quality of representation
+              gradient.cols = c("blue", "red"),
+              repel = TRUE     # Avoid text overlapping
+ )

                                                                  ## Distance Analysis ##
> treeUPGMA  <- upgma(dm)
> plot(treeUPGMA, main="UPGMA")
> treeNJ  <- NJ(dm)
> library(phylocanvas)
> phylocanvas(treeNJ,treetype = "rectangular", alignlabels = T)
> library(ape)
> library(phangorn)
> parsimony(treeUPGMA, ort)
> treeRatchet  <- pratchet(ort, trace = 0, minit=100)
> parsimony(treeRatchet, ort)
> treeRatchet  <- acctran(treeRatchet, ort)
> if(inherits(treeRatchet, "multiPhylo")){
+     treeRatchet <- unique(treeRatchet)
+ }
> plotBS(midpoint(treeRatchet), type="phylogram")
> phylocanvas(treeRatchet,treetype = "rectangular", alignlabels = T)

> mt <- modelTest(ort, model=c("K80", "HKY", "GTR"),
+ control = pml.control(trace = 0))
> fit <- as.pml(mt, "BIC")
> bs <- bootstrap.pml(fit, bs=5, optNni=TRUE,
+ control = pml.control(trace = 0))