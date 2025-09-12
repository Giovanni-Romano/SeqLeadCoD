library(igraph)
library(ggraph)
library(tidyverse)

adj2graph = function(A, comm){
  g = graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  V(g)$ID = rownames(A) 
  V(g)$community = comm[match(rownames(A), names(comm))]
  V(g)$size = degree(g)
  info.R = ghe %>% select(ID, Region) %>% unique()
  V(g)$Region = (info.R$Region)[match(V(g)$ID, info.R$ID)]
  info.S = ghe %>% select(ID, Sex) %>% unique()
  V(g)$Sex = (info.S$Sex)[match(V(g)$ID, info.S$ID)]
  
  return(g)
}

layout_communities = function(g.ogpos, weight.comm){
  g_grouped = g.ogpos
  for(i in unique(V(g_grouped)$community)) {
    GroupV = which(V(g_grouped)$community == i)
    lg = length(GroupV)
    w = weight.comm
    if (lg > 2){
      g_grouped = add_edges(g_grouped, combn(GroupV, 2), attr=list(weight=w))
    }
  }
  out = layout_with_graphopt(g_grouped, charge = 0.01, mass = 1)#layout_with_fr(g_grouped)
  
  
  return(out)
}


plotgraph = function(HSM, thres.plot, thres.pos, comm, weight.comm, print.text = F,
                     positions = NULL){
  
  if (is.null(positions)){
    adj.pos = {
      thres = thres.pos
      tmp = HSM
      tmp[tmp <= thres] = 0
      tmp[tmp > thres] = 1
      tmp
    }
    g.pos = adj2graph(adj.pos, comm)
    LO = layout_communities(g.pos, weight.comm)
  } else {
    LO = positions
  }
  
  
  adj.plot = {
    thres = thres.plot
    tmp = HSM
    tmp[tmp <= thres] = 0
    tmp[tmp > thres] = (tmp[tmp > thres] - thres) / (19 - thres)
    tmp
  }
  
  g.plot = adj2graph(adj.plot, comm)
  
  V(g.plot)$x = LO[,1]
  V(g.plot)$y = LO[,2]
  
  # Create a data frame of node positions and communities
  node_df = data.frame(
    x = V(g.plot)$x,
    y = V(g.plot)$y,
    community = V(g.plot)$community,
    Region = V(g.plot)$Region,
    Sex = V(g.plot),
    ID = V(g.plot)$ID
  )
  
  
  plt = ggraph(g.plot, layout = "manual", x = V(g.plot)$x, y = V(g.plot)$y) +
    geom_edge_link(aes(edge_alpha = weight), edge_color = "grey20", edge_width = 0.5, show.legend = F) +
    geom_node_point(aes(fill = Region, shape = Sex), color = "black", size = 3) +
    scale_fill_manual(values = c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")) +
    guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
    scale_shape_manual(values = c(21, 22)) +
    ggforce::geom_mark_hull(data = node_df, aes(x = x, y = y, group = community),
                            fill = "lightblue", colour = "black",
                            alpha = 0.2, concavity = 5, expand = unit(2, "mm"), show.legend = FALSE)
  
  if (print.text) {
    plt = plt + 
      geom_text(aes(label = ID), repel = T, size = 3, vjust = 1.5, hjust = 0.5, color = "black") +
      theme(legend.position = "none")
  } else {
    plt = plt + theme(legend.position = "right")
  }
  
  return(plt)
}
