library(stdbscanr)
library(tidyverse)
set.seed(0) # set random seed to get repeatable results

# make some trajectory data
dt <- data.table(X = rnorm(40),
                 Y = rnorm(40),
                 timestamp = 1:40)
#dt=df
# find clusters and density in it
out <- get_clusters_from_data(dt
                              ,x = "X"
                              ,y = "Y"
                              ,t = "timestamp"
                              ,eps = 1 # 1 metre distance threshold from other point
                              ,eps_t = 5 # 5 second time threshold from other point
                              ,minpts = 2)[  # minimum connected to 3 points to continue growing a cluster
                                ,cluster:=factor(cluster)] # make clusters a factor

#plot clusters output
ggplot(out,
       aes(x = X,
           y = Y,
           label = timestamp,
           color = cluster,
           alpha =))+
  geom_path(data = out[!is.na(cluster),],
            alpha = 0.5)+
  geom_path(data = out[is.na(cluster),],
            alpha = 0.2)+
  geom_text()+
  labs(title = "ST-DBSCAN Demo",
       subtitle = "minpts = 3, eps = 1, eps_t = 5")




get_clusters_from_data <- function(df
                                   ,x = "X"
                                   ,y = "Y"
                                   ,t = "timestamp"
                                   ,eps = 1 # metre from other point
                                   ,eps_t = 5 # seconds from other point
                                   ,minpts = 3){

  df=dt
  # set up input data - set it to be data.table
  data.table::setDT(df)

  #sort by time
  data.table::setkeyv(df,cols = t)

  #create id column
  df[,id:=.I]

  # get the right columns
  dt = df[,c("id", x, y, t),with = F]

  #name them how the other functions expect
  data.table::setnames(dt,
                       old = c("id",   x,   y,  t),
                       new = c("id", "x", "y", "time"))


  # find reachables
  reachables <- stdbscanr:::find_connected_points(dt = dt,
                                                  eps = eps,
                                                  eps_t = eps_t)


  point_density <- reachables[,.(density = .N),by = first]

  data.table::setnames(point_density,old = c("first"),new = c("id"))


  # remove noise (clusters with less than minpts)
  reachables_without_noise <- stdbscanr:::remove_noise(reachables = reachables,
                                                       min_number = minpts)


  #work out core and terminating points for debugging
  core_pts = unique(reachables_without_noise$first)
  terminating_pts = setdiff(unique(reachables_without_noise$second), core_pts)
  point_type = data.table(point = c(core_pts,
                                    terminating_pts),
                          point_type = c(rep("core", length(core_pts)),
                                         rep("terminating", length(terminating_pts))))

  # initially set cluster_id to first index
  clusters_core <- copy(reachables_without_noise[second %in% core_pts,])[,cluster:=first]
  clusters_terminating <- copy(reachables_without_noise[second %in% terminating_pts,])

  #iterate on core points until final clusters are here

  #get intermediate clusters
  iter_count <- 0
  clustersList<-list()
  while (T){
    #some debugging output - this can be removed later.
    iter_count <- iter_count + 1

    #store a copy of before the iteration to compare to afterwards
    clusters_core_old <- copy(clusters_core)
    clustersList[[iter_count]] <- clusters_core
    #iterate giving the cluster labels of firsts to their seconds
    #clusters_core <- stdbscanr:::cluster_iteration(clusters = clusters_core_old)
    clusters_core <- clusters_core[ , .(first,cluster=min(cluster)), by = second]
    clusters_core <- clusters_core[ , .(second,cluster=min(cluster)), by = first]
    #if nothing changed in an iteration, we are done!
    if (nrow(clusters_core_old)==nrow(clusters_core)){
      if (T==all.equal(clusters_core_old, clusters_core)){
        deduped_clusters_core <- clusters_core[,.(count = .N),by = list(second,cluster)]
        if (nrow(deduped_clusters_core)==length(unique(deduped_clusters_core$second))) {#also check that each second only has one cluster number
          break
          }else{
            stop("Error - multiple clusters found for each core point")
          }
      }
    }
  }

  deduped_clusters_core <- deduped_clusters_core[,.(second,cluster)]

  #deal with terminating points
  #get cluster number for every second from each core point it is connected to
  clusters_terminating2 <-
    data.table::merge.data.table(clusters_terminating,
                                 deduped_clusters_core,
                                 by.x = "first",
                                 by.y = "second",
                                 suffixes = c("1","2"),
                                 all.x = T)

  #keep only only min cluster number for each second
  clusters_terminating3 <-
    clusters_terminating2[
      clusters_terminating2[ , .I[which.min(cluster)],
                             by = second]$V1]


  clusters <- rbindlist(list(deduped_clusters_core,
                             clusters_terminating3[,.(second,cluster)]),
                        use.names = T)

  get_clus<-function(x){clustersList[[x]][,3]}
  clus <- bind_cols(clustersList[[1]][,2],lapply(1:length(clustersList),get_clus))
  names(clus)<-c("second",sapply(1:length(clustersList),function(x)paste0("cluster",x)))
  deduped_clus <- clus[,`:=`(second = paste0("Original pt. ",second)
                             ,cluster1 = paste0("First cluster ",cluster1)
                             ,cluster2 = paste0("Second cluster ",cluster2)
                             ,cluster3 = paste0("Final cluster ",cluster3))]
                            #,cluster4 = paste0("Fourth cluster ",cluster4))]


  connections<-
    rbindlist(list(clusters_terminating3[,.(second = paste0("Original pt. ",second),
                                            cluster = paste0("Final cluster ",cluster))],
                   deduped_clus[,.(second, cluster1)],
                   deduped_clus[,.(cluster1, cluster2)],
                   deduped_clus[,.(cluster2, cluster3)]#,
                   #deduped_clus[,.(cluster3, cluster4)]
    ),
    use.names = FALSE)


  nodevec <- unique(c(connections$second,connections$cluster))
  links = data.table(start = as.numeric(factor(connections$second, levels = nodevec))-1,
                     finish = as.numeric(factor(connections$cluster, levels = nodevec))-1,
                     val = 1)

  nodes = data.table(id = as.numeric(factor(nodevec, levels = nodevec))-1,
                     name = nodevec)

  networkD3::sankeyNetwork(Links = links,
                           Nodes = nodes,
                           Source = "start",
                           Target = "finish",
                           NodeID = "name",
                           NodeGroup = "name",
                           Value = "val",
                           fontSize = 12)

  #clusters <- find_equilibrium_clusters(clusters = clusters)

  #remove seconds with two different clusters - take the first
  #(this is becasue terminating points can belong to two disconnected clusters)
  deduped_clusters <- clusters[,.SD[1],by = second]

  #rename the clusters to be sequential starting from 1
  cluster_numbers_old <- na.omit(unique(deduped_clusters$cluster))
  cluster_numbers_new <- 1:length(cluster_numbers_old)
  names(cluster_numbers_new) <- cluster_numbers_old

  #extract just the point ids and the cluster they belong to
  point_to_clust <- unique(deduped_clusters[,.(point_id = second,
                                               cluster = cluster_numbers_new[as.character(cluster)])])

  #relabel the clusters to be sequential starting from 1
  temp <- data.table::merge.data.table(x = df,
                                       y = point_to_clust[,.(point_id,
                                                             cluster)],
                                       by.x = "id",
                                       by.y = "point_id",
                                       all.x = T)

  # give each point a type - core, terminating or noise.
  temp2 <- data.table::merge.data.table(x = temp,
                                        y = point_type,
                                        by.x = "id",
                                        by.y = "point",
                                        all.x = T)
  temp2$point_type[is.na(temp2$point_type)] <- "noise"
  temp2$point_type <- factor(temp2$point_type,
                             ordered = T,
                             levels = c("core","terminating","noise"))

  df_out <- data.table::merge.data.table(x = temp2,
                                         y = point_density,
                                         by.x = "id",
                                         by.y = "id",
                                         all.x = T)

  # test we get one row out per row in...
  if (nrow(df_out)!=nrow(df)){
    warning("Incorrect number of rows output - need to debug!")
  }

  return(df_out[,id:=NULL])
}

