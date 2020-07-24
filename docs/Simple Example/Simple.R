library(stdbscanr)
set.seed(0) # set random seed to get repeatable results

# make some trajectory data
dt <- data.table(X = rnorm(20),
                 Y = rnorm(20),
                 timestamp = 1:20)

# find clusters and density in it
out <- get_clusters_from_data(dt
                              ,x = "X"
                              ,y = "Y"
                              ,t = "timestamp"
                              ,eps = 1 # 1 metre distance threshold from other point
                              ,eps_t = 5 # 5 second time threshold from other point
                              ,minpts = 3)[  # minimum connected to 3 points to continue growing a cluster
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
