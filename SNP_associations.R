options(stringsAsFactors = F)
library(RSQLite) 
library(qtl2)

query_variants <- create_variant_query_func("cc_variants.sqlite")
query_genes <- create_gene_query_func("mouse_genes_mgi.sqlite")

###############################################################################
load("E:/DO/R-Analyses/Nov_2018/DO_Lipidomics.RData") 

tis = "plasma"
#Select normalized dataset
if(tis=="plasma") {
  pheno <- dataset.plasma.lipids$pheno
} else if (tis=="liver") {
  pheno <- dataset.liver.lipids$pheno
}

qtl.summary <- subset(qtl.summary,annot.id %in% colnames(pheno))

all=NULL

for(i in c(1:nrow(qtl.summary))) {
  
  # Get the QTL info for this peak.
  qtl.info = qtl.summary[i,]
  analyte = qtl.info$annot.id
  chr = qtl.info$qtl.chr
  pos = qtl.info$qtl.pos
  
  variants <- query_variants(chr, pos - 1, pos + 1)
  genes <- query_genes(chr, pos - 1, pos + 1)
  
  # Association mapping.
  start = max(1, pos - 2)
  end = pos + 2
  
  assoc = scan1snps(genoprobs = genoprobs[,chr], map = map, 
                    pheno = pheno[,analyte, drop = FALSE], 
                    kinship = K[[chr]], intcovar = dataset.liver.lipids$covar, query_func=query_variants,
                    chr = chr, start = start, end = end, keep_all_snps=TRUE)
  
  top <- top_snps(assoc[[1]], assoc[[2]], drop = 1.5, show_all_snps = T)
  
  if(nrow(top)>0){
    top$analyte <- analyte
    if(length(all)>0){all <- rbind(all,top)} 
    else {all <- top}
  }
}

write.csv(all, file=paste0(tis,"_snps.csv"))