# MINGLE
This folder contains codes to reproduce the network-based variable selection process described in ... (Figure 1) 

Files are numbered based on the order of execution.

1 - Joint_graphical_lasso: first step of the defined framework. It loads RNA-sequencing data as input (RNAseq_data) and jointly estimates the networks for the three glioma types. The regularization parameters used here are lambda1=0.9 and lambda2=0.001. For personalized analysis, these parameters can be changed. 

2 - Creation_starting_methylomics_dataset: intermediate step needed to create the methylomics dataset to use for the second step of our study. It loads as input the results of the JGL (previous step), the complete DNA-methylation dataset (DNAMeth-completeDS) and the RNASeq_data, to match samples across the two omics, and the map between the two layers (map-genes-probes), to associate the CpG sites with the selected genes. The output is DNAMeth_data-ExEdge_RNAseq, containing the CpG sites associated with the genes involved in exclusive relations in the estimated transcriptomics networks.

3 - Graphical_lasso: second step of the defined framework. It loads DNAMeth_data-ExEdge_RNAseq data as input and estimates the network for a given glioma type. The desired type can be selected by changing the value of the variable GT, which can be 1 (astrocytoma), 2 (oligodendroglioma), or 3 (GBM). The regularization parameter used is rho=0.85. For personalized analysis, this parameter can be changed. 

4 - MINGLE: third and last step of the defined framework. It is an R function, which needs as input:
  - net: the methylomics network used to compute the edges (results_glasso-rho085-meth-gbm)
  - map: the functional map between genes and CpG sites (map-genes-probes)
  - edge.type: a label specifying which type of edge want to be used (strength or count, see Coletti et al. for details)
  - weights: if "count" is selected as edge.type, it must be a vector of the same dimension as the number of nodes of the MINGLE network. If "strength" is selected as edge.type, it might be 1 (or any other number).
Based on these input data, the MINGLE function computes the new network of genes based on methylomics relations.

For correct usage of these scripts, we suggest to download the dataset files available as "realises", and to store these in a folder named "dataset", saved in the same location as the R scripts.
