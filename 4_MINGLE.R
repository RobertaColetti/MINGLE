# title: ------ Create MINGLE Networks -------
# author: Roberta Coletti

# This script includes MINGLE function 
#(ref.: Coletti et al., A Novel Tool for Multi-Omics Network Integration and Visualization: A Study of Glioma Heterogeneity 2024)

# MINGLE function inputs:
# net: the underlying network used as reference to compute the edges (methylomics networks in Coletti et al.). Nodes of this network will be referred as "Nodes of underlying layer" (NodesUL)
# map: the functional link between the NodesUL and the new MINGLE Nodes (map between CpG and genes in Coletti et al.)
# edge.type: label indenfitying the desired function to compute the edge weights. There are two possibilities:
# 
# 1) strength: edges are computed by the median value of the weights of the links of all the probes related to a given gene;
#
# 2) count: edges are computed by counting the number of connection of all the probes related to a given gene. It can be weighted or not. 
#           If "count" is selected, MINGLE needs as input a vector of weights which dimension must match with the number of nodes of the MINGLE network.
# In Coletti et al. the MINGLE networks are constituted by the genes associated to the selected CpGs, and the weights the percentages of selected probes.
# 


MINGLE <- function(net, map, edge.type,weights){

  #Save the Node of the underlying network
  NodesUL=colnames(net)
  #Save the MINGLE nodes
  M.nodes=map[which(map[,1] %in% NodesUL),2]
  M.nodes=M.nodes[!duplicated(M.nodes)]
  
  #check the input are correct:
  if(edge.type != "strength" AND edge.type !="count")
  {  print("Error! Invalid edge.type.")
    return()
    #geterrmessage()
  }
  if (edge.type=="count")
  {
    if (length(weights)!=length(M.nodes))
    {  print("Error! Count measure has been selected; weights must be set.")
      return()
      #geterrmessage()
    }
  }
  
  #compute the MINGLE networks
  MINGLE.net=matrix(0, length(M.nodes),length(M.nodes))
  rownames(MINGLE.net)=colnames(MINGLE.net)=M.nodes
  last=length(M.nodes)-1
  for (i in c(1:last) )
  { first=i+1 #Only the upper triangular matrix is created
  for (j in c(first:length(M.nodes))) 
    #identification of the NodesUL linked to the i-th MINGLE node 
  { #rows of the NodesUL linked with i-th MINGLE node
    rows_NodesUL.i=which(map[,2] %in% M.nodes[i])
    #find the indexes of the adjacency matrix corresponding to the NodesUL associated to the i-th MINGLE node
    index_i=which(NodesUL %in% map[rows_NodesUL.i,1])
    #same with j-th MINGLE node
    rows_NodesUL.j=which(map[,2] %in% M.nodes[j])
    index_j=which(NodesUL %in% map[rows_NodesUL.j,1])
    
    vectorized_net=net[index_i,index_j]
    vectorized_net=vectorized_net[vectorized_net!=0]
    if (edge.type =="count")
    {
      w=(weights[i]+weights[j])/2
      MINGLE.net[i,j]=length(vectorized_net)*w
    }
    if(edge.type =="strength")
    { 
      MINGLE.net[i,j]=median(abs(vectorized_net))
      MINGLE.net[is.na(MINGLE.net)] <- 0
    }
  }
  }
  MINGLE.net=MINGLE.net+t(MINGLE.net)
 return(MINGLE.net) 
}


# Required libraries
library("qgraph")
library("base")

# ------------- EXAMPLE OF EXECUTION OF MINGLE ----------------- 
# ----- Script to reproduce the output of Coletti et al. -----


#1) Load data:
#methylomics glasso results (related to GBM):
load("~/Network-discovery-results/results_glasso-rho085-meth-gbm.RData") 

#map CpG --> Genes
load("~/map-genes-probes.RData")

#reduce the map to the starting GBM dataset:
red.CpG2Gene=CpG2Gene[which(CpG2Gene$CpG %in% colnames(gbm_met_i)),]
associated_genes=red.CpG2Gene$Gene[which(red.CpG2Gene$CpG %in% selected_meth)]
associated_genes=associated_genes[!duplicated(associated_genes)]

# ---- if Count measure is selected ---- 
# Weight computing:
# determining the percentage of selected probes for each gene 
percentage_sel_probes=matrix(0,nrow = length(associated_genes), ncol = 2) 
rownames(percentage_sel_probes)=associated_genes
colnames(percentage_sel_probes)=c("percentage","Tot_probes")
for(j in c(1:length(associated_genes))) #for each gene...
{  #count the total number of CpG sites 
  index=which(red.CpG2Gene$Gene %in% associated_genes[j]) 
  Total= length(index)
  #count the number of the selected probes 
  Selected= length(which(colnames(Sel_Theta)%in% red.CpG2Gene$CpG[index]))
  percentage_sel_probes[j,1]=Selected/Total
  percentage_sel_probes[j,2]=Total
}


MINGLE.net=MINGLE(Sel_Theta,red.CpG2Gene,edge.type = "strength",1)
#MINGLE.net=MINGLE(Sel_Theta,red.CpG2Gene,edge.type = "count", percentage_sel_probes[,1])

#Graphical part:
#Remove isolated nodes for graphical representation (i.e. genes that are selected for self-related probe relations)
MINGLE.net=MINGLE.net[apply(MINGLE.net, 2, var, na.rm=TRUE) != 0,apply(MINGLE.net, 2, var, na.rm=TRUE) != 0]
MINGLE.nodes=colnames(MINGLE.net)

qgraph(MINGLE.net, layout = "spring", 
       directed = FALSE, edge.labels = F, 
       esize = 8,
       labels=MINGLE.nodes,
       label.font=12,
       label.prop=1,
       legend = F,
       border.color = "black",border.width = 2,
       edge.width = 0.8,
       label.color = "black" )
title(sub = "Strength MINGLE network") #if edge.type = "strength"
#title(sub = "Count MINGLE network") #if edge.type = "count"



