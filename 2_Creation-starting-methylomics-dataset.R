#This script creates the starting methylation dataset for each glioma type by associating the exclusive pairs of genes (genes linked by relations that are exclusive for a glioma type) 

#load the desired JGL outcome:
load("~/Network-discovery-results/JGL-lam109-lam20001.RData")

# --- Matrix preparation: separate the three networks (one per each glioma type):
Theta_a=fgl.results$theta[[1]]
Theta_o=fgl.results$theta[[2]]
Theta_g=fgl.results$theta[[3]]

Theta_a=Theta_a-diag(dim(Theta_a)[1])*Theta_a
Theta_o=Theta_o-diag(dim(Theta_o)[1])*Theta_o
Theta_g=Theta_g-diag(dim(Theta_g)[1])*Theta_g

#put 1 if there is a link between two variables
Theta_a[Theta_a!=0] <-1
Theta_o[Theta_o!=0] <-1
Theta_g[Theta_g!=0] <-1

# --- Dimensionality reduction: detect the exclusive relations removing the ones in Astro and Oligo
Ex_GBM= Theta_g - Theta_a - Theta_o
Ex_Astro= Theta_a - Theta_g - Theta_o
Ex_Oligo= Theta_o - Theta_a - Theta_g

#This operation lead to matrices with only 3 possible values:
# 0 --> No edge in all glioma types OR shared edge with one other type (e.g., Ex_GBM = Theta_g - Theta_a - Theta_o = 0 - 0 - 0 OR 1 - 1 - 0)
# 1 --> Exclusive edge (the edge only exists in the first matrix: e.g., Ex_GBM = Theta_g - Theta_a - Theta_o = 1 - 0 - 0)
#-1 --> Exclusive edge in another glioma type (e.g., Ex_GBM = Theta_g - Theta_a - Theta_o = 0 - 1 - 0)
#-2 --> Shared edge in the two other glioma types (e.g., Ex_GBM = Theta_g - Theta_a - Theta_o = 0 - 1 - 1)

#Since we are interested in EXCLUSIVE edge, we only want to consider the pairs associated to 1:
ex_pairs_GBM=colnames(Ex_GBM)[which(Ex_GBM == 1, arr.ind=T)[, "col"]]
ex_pairs_GBM=ex_pairs_GBM[!duplicated(ex_pairs_GBM)]
ex_pairs_Astro=colnames(Ex_Astro)[which(Ex_Astro == 1, arr.ind=T)[, "col"]]
ex_pairs_Astro=ex_pairs_Astro[!duplicated(ex_pairs_Astro)]
ex_pairs_Oligo=colnames(Ex_Oligo)[which(Ex_Oligo == 1, arr.ind=T)[, "col"]]
ex_pairs_Oligo=ex_pairs_Oligo[!duplicated(ex_pairs_Oligo)]

# --- Creation of methylomics starting dataset: association of the genes with CpG sites 
load("~/map-genes-probes.RData")
starting_probes_GBM=CpG2Gene$CpG[which(CpG2Gene$Gene %in% ex_pairs_GBM)]
starting_probes_Astro=CpG2Gene$CpG[which(CpG2Gene$Gene %in% ex_pairs_Astro)]
starting_probes_Oligo=CpG2Gene$CpG[which(CpG2Gene$Gene %in% ex_pairs_Oligo)]

#load the methylomics dataset (samples classified according to 2021-WHO)
load("~/dataset/DNAMeth-completeDS.RData")
#load the RNASeq dataset (to match samples)
load("~/dataset/RNAseq_data.RData")

gbm_samples=rownames(GBM_RNA)
oligo_samples=rownames(oligo_RNA)
astro_samples=rownames(astro_RNA)

#Create the initial dataset for each glioma type by matching samples (accross omics) and variables (exclusive pairs)
astro_DNAMeth=astro_met[which(rownames(astro_met) %in% astro_samples),which(colnames(astro_met) %in% starting_probes_Astro)]
oligo_DNAMeth=oligo_met[which(rownames(oligo_met) %in% oligo_samples),which(colnames(oligo_met) %in% starting_probes_Oligo)]
GBM_DNAMeth=gbm_met[which(rownames(gbm_met) %in% gbm_samples),which(colnames(gbm_met) %in% starting_probes_GBM)]

#save(astro_DNAMeth,oligo_DNAMeth,GBM_DNAMeth, file= " DNAMeth_data-ExEdge_RNAseq.RData") #uncomment if you want to save your outcomes
