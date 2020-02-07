
################################################################################################

## Given a node, tree, and distance matrix,
## This function returns median pairwise patristic distance (MPPD) of all of its decendants
get.node.full.MPPD <- function(node,tree,distmat){
  nlist <- tips(tree,node)
  elist <- tree$edge[which.edge(tree,nlist),2]
  foo <- distmat[elist,elist]
  return(median(foo[upper.tri(foo,diag=FALSE)]))
}
################################################################################################

## Given a tree and (optionally) a distance matrix,
## This function returns a vector giving the median pairwise patristic distance 
## of the subtree under each internal node
pdist.clusttree <- function(tree,distmat){
  
  ntips<- Ntip(tree)
  nint <- tree$Nnode ## number of internal nodes
  return(sapply((ntips+1):(ntips+nint), function
                (x) {get.node.full.MPPD(x, tree,distmat)} ))
}
################################################################################################

## Given a tree and a threshold, and (optionally) a vector of the MPPD of the 
## subtree at each internal node, the reliability threshold and vector,
## and specify if you want the membership for all nodes or only the tips
## This function returns a vector indicating, for each internal node, which
##  cluster it belongs to
prosperi.cluster <- function(tree,thresh,distvec=NULL,
                             rthresh=NULL,reliabilityvec=NULL){
  
  ## set up clustering
  ntips<-Ntip(tree)
  cnum <- 0 ## cluster number
  assign <- rep(0,ntips+tree$Nnode) ## cluster assignment
  igraph.tree <- graph.edgelist(tree$edge) ## tree in igraph form
  dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out',
                   order=TRUE,dist=TRUE)
  ## travese the tree in depth first order
  for(i in 1:length(dfs$order)){
    node <- dfs$order[i]
    ## skip leaves
    if(node < ntips+1){ next }
    ## skip unreliable nodes (if reliability measure is available)
    if(! is.null(reliabilityvec) &&
         reliabilityvec[node-ntips] >= rthresh){ next }
    ## If the node's subtree is below the threshold, mark it and
    ## its subtree as members of a new cluster
    if(distvec[node-ntips]<=thresh && assign[node]<=0){
      cnum <- cnum+1
      subtree <- graph.dfs(igraph.tree,node,
                           neimode='out',unreachable=FALSE)$order
      assign[subtree] <- cnum
    }}
  ans <- list(membership=assign,allcsize=table(assign),
              leafclustsize=table(assign[1:ntips]),
              ntips=ntips,threshold=thresh)
  class(ans) <- c(class(ans),'p.cluster')
  return(ans)
}
################################################################################################

# Citation: Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete 
# character-taxon matrices: implementation, challenges, progress, and future directions. 
# Biological Journal of the Linnean Society, 118, 131-151.

FindAncestor <- function(descs, tree) {
  # Get tip numbers:
  tipnos <- match(descs, tree$tip.label)
  # Get ancestral nodes in order:
  anc.node <- sort(unique(tree$edge[, 1][match(tipnos, tree$edge[, 2])]))
  # Keep going until a single ancestral node is converged upon:
  while(length(anc.node) > 1) {
    
    # Get node with highest number (definitely not ancestor):
    highestnode <- anc.node[length(anc.node)]
    # Remove this node from the list:
    anc.node <- anc.node[-length(anc.node)]
    # Find its ancestor and add to unique list:
    anc.node <- sort(unique(c(anc.node, tree$edge[match(highestnode, tree$edge[, 2]), 1])))
  }
  # Return ancestral node:
  return(anc.node)
}
################################################################################################
