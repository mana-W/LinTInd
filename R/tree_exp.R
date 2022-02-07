library(DiagrammeR)



tree_nwk = ToNewick(treeinfo$tree)
write(tree_nwk,"tree.nwk")




node_dataframe<-ToDataFrameNetwork(treeinfo$tree)
clus<-unique(treeinfo$info$celltype)
node_tab<-data.frame(acast(treeinfo$info,tags~celltype))
node_tab<-data.frame(t(apply(node_tab,1,function(x){x/sum(x)})))
node_tab["N0",]<-rep(0,dim(node_tab)[2])
from_node<-node_dataframe$from[match(row.names(node_tab),node_dataframe$to)]
from_node_tab<-node_tab[c(from_node),]

node_tab_matrix<-as.matrix(node_tab)
from_node_tab_matrix<-as.matrix(t(from_node_tab))
netdf<-from_node_tab_matrix%*%node_tab_matrix

from_define<-function(x){
  max_index<-which(netdf[,x]==max(netdf[,x]))
  if(length(max_index)==1){
    p_node<-row.names(netdf)[max_index]
    if(p_node==x){
      return(paste0(x,"_pre"))
    }else{
      return(p_node)
    }
  }else{
    return("N0")
  }
}

to_list<-colnames(netdf)
from_list<-sapply(to_list,from_define)
df<-data.frame("from"=from_list,"to"=to_list)
df<-rbind(df,data.frame(from="N0",to=df$from[grep("_pre",df$from)]))
celltree<-FromDataFrameNetwork(df)

plot(celltree)
diagonalNetwork(ToListExplicit(celltree))




from_define<-function(x){
  max_index<-which(data_p[,x]<=0.05)
  if(length(max_index)>=1){
    return(row.names(data_p)[max_index])
  }else{
    return("N0")
  }
}

from_define<-function(x){
  max_index<-which(ds[,x]==min(ds[,x]))
  if(length(max_index)==1){
    return(row.names(ds)[max_index])
  }else{
    return("N0")
  }
}

to_list<-colnames(ds)
from_list<-sapply(to_list,from_define)
df<-data.frame("from"=from_list,"to"=to_list)
celltree<-FromDataFrameNetwork(df)
plot(celltree)

