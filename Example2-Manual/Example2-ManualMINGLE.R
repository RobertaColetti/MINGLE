#other examples of MINGLE application:

# ------ Example2: transcriptomics/proteomics -----

#  Example 2a: MINGLE protein network from transcriptomics
#----- library section
library("qgraph")


#1) Data loading
# transcriptomics network from Coletti and Lopes (2023): 
load("~/Example2-Manual/net_transcriptomics.RData") #loaded files: Sel_Theta

#map genes <--> proteins:
load("~/Example2-Manual/Gene-protein-map.RData") #loaded file: Gene2Prot

#2) Data preparation
#reduce the map to the genes of the transcriptomics network:
genes = colnames(Sel_Theta)
red.map = Gene2Prot[which(Gene2Prot$Gene %in% genes),]

# 3) Call MINGLE function
MINGLE.net = MINGLE(net = Sel_Theta, 
                    map = red.map, 
                    edge.type = "strength", 
                    weights = 1)

# 4) Graphical visualization
MINGLE.nodes=colnames(MINGLE.net)

#Graph:

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
title(sub = "Weight MINGLE network")
dev.off()


#  Example 2b: Comparison between MINGLE network and proteomics network
#----- library section
library("qgraph")

#1) Data loading
# proteomics network from Coletti and Lopes (2023): 
load("~/Example2-Manual/net-proteomics.RData") #loaded files: Sel_Theta
proteins = colnames(Sel_Theta)

#2) Graphical visualization:
# save common nodes: 
common = intersect(proteins, MINGLE.nodes)

#Graphs:
qgraph(MINGLE.net, layout = "spring", 
       directed = FALSE, edge.labels = F, 
       esize = 8,
       labels=MINGLE.nodes,
       label.font=12,
       label.prop=1,
       groups = list(Common_nodes =which(colnames(MINGLE.net) %in% common)),
       color = c("#ffe44d"),
       legend = F,
       border.color = "black",border.width = 2,
       edge.width = 0.8,
       label.color = "black" )
title(sub = "Weight MINGLE network")

qgraph(Sel_Theta, layout = "spring", 
       directed = FALSE, edge.labels = F, 
       esize = 8,
       labels=proteins,
       label.font=12,
       label.prop=1,
       groups = list(Common_nodes =which(colnames(Sel_Theta) %in% common)),
       color = c("#ffe44d"),#,"#ffb833","#39ac6b","#87CEFA"),
       legend = F,
       border.color = "black",border.width = 2,
       edge.width = 0.8,
       label.color = "black" )
title(sub = "Proteomics network")
dev.off()


