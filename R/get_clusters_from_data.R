#' Find temporally close points
#'
#' Find close times less than `eps_t` away from each other in an ordered vector of times.
#'
#' @param vec A sorted vector e.g. of timestamps.
#' @param eps_t The max distance between two points for them to be "connected".
#' @return A data.table with two columns, "first" and "second". These are the indicies of the connected point, and the points it is connected to.
#' @examples
#' find_temporally_connected_points(vec = 1:10, eps_t = 1.5)
find_temporally_connected_points<-function(vec, eps_t){

  sparse_sim <- vector("list", length(vec))

  for (i in 1:length(vec)){

    #count forwards
    fcounter <- 0
    while (T){
      if (length(vec) < i + fcounter + 1) break
      if (abs(vec[i + fcounter + 1] - vec[i]) <= eps_t) break
      fcounter <- fcounter + 1
    }


    #count backwards
    bcounter <- 0
    while (T){
      if (0 >= i - bcounter - 1 ) break
      if (abs(vec[i] - vec[i - bcounter - 1])  >= eps_t) break
      bcounter <- bcounter + 1
    }


    #write indicies of elements within range as a list element
    sparse_sim[[i]] <- (i-bcounter):(i+fcounter)



  }

  #collapse the sparse simmilarity list into a data.table
  names(sparse_sim) <- 1:length(vec)
  sparse_sim <- utils::stack(sparse_sim)
  data.table::setDT(sparse_sim)

  #return just the columns we need
  return(
    sparse_sim[,
               ind:=as.integer(as.character(ind))
               ][,.(first = ind, second = values)]
  )
}



#' Find spatially close points
#'
#' Find connected points less than `eps` away from each other in space,
#' given a list of those points which are sufficiently connected in time.
#'
#' @param dt data.table of points with columns `id`, `x`, `y`.
#' @param sls The sparse logical similarity matrix output from `find_temporally_connected_points()`.
#' @param eps The max distance between two points for them to be "connected".
#' @return A data.table with two columns, "first" and "second". These are the indicies of
#' the connected points, and the points it is connected to.
#' @examples
#' dt <- data.table(id = 1:100,
#'                 x=rnorm(100),
#'                 y = rnorm(100),
#'                 t = 1:100)
#' find_spatially_connected_points(dt = dt,
#'                       sls = find_temporally_connected_points(vec = dt$t,
#'                                              eps_t = 1.5),
#'                       eps = 10)
find_spatially_connected_points <- function(dt, sls, eps){

  #merge with coordinates of first point, second point
  int<-data.table::merge.data.table(sls, dt, by.x = "first", by.y = "id",suffixes = c("",".first"),all.x = T)
  fin<-data.table::merge.data.table(int, dt, by.x = "second", by.y = "id",suffixes = c("",".second"),all.x = T)

  #compute distance and check if under threshold
  fin[,dist:=sqrt((x-x.second)^2+(y-y.second)^2)][,close:=dist<=eps]
  #output only "reachable" points
  out<- fin[close==T,.(first,second)]
  return(out)
}

#' Find connected points less than `eps` away from each other in space,
#' and less than `eps_t` away from each other in time, given a data.table of those points
#' given a list of those points which are sufficiently connected in time.
#'
#' @param dt data.table of points with columns `id`, `x`, `y`, `time`.
#' @param eps The max distance between two points for them to be "connected".
#' @param eps_t The max time between two points for them to be "connected".
#' @return A data.table with two columns, "first" and "second". These are the indicies of
#' the connected points, and the points it is connected to.
#' @examples
#' dt <- data.table(id = 1:100,
#'                 x=rnorm(100),
#'                 y = rnorm(100),
#'                 t = 1:100)
#' find_connected_points(dt = dt,
#'                sls = find_temporally_connected_points(vec = dt$t,
#'                                       eps_t = 1.5),
#'                eps = 10)
find_connected_points <- function(dt, eps, eps_t){

  sls <- find_temporally_connected_points(vec = dt$time,
                                          eps_t = eps_t)

  reachables <- find_spatially_connected_points(dt = dt,
                                                sls = sls,
                                                eps = eps)

  return(reachables)
}


#' Remove noise
#'
#' Remove connections from a 'first' point if that point doesn't have enough points
#' around it (`min_points`) that it is connected to.
#'
#' @param reachables data.table of connections with columns `first`, `second` which
#' has been output from `find_connected_points()`
#' @param min_number The minimum number of connected `second` points for this `first` point to be kept.
#' @return A data.table with two columns, `first` and `second`. These are the indicies of
#' the connected points, and the points it is connected to.
#' @examples
#' dt <- data.table(id = 1:100,
#'                 x=rnorm(100),
#'                 y = rnorm(100),
#'                 t = 1:100)
#' reachables <- find_connected_points(dt = dt,
#'                sls = find_temporally_connected_points(vec = dt$t,
#'                                       eps_t = 1.5),
#'                eps = 10)
#' remove_noise(reachables = reachables, min_number = 4)
remove_noise <- function(reachables, min_number = 4){
  #remove all firsts which don't have min_number seconds
  noise_free_reachables <- reachables[,if(.N>=min_number) .SD,by = first]

  if(nrow(noise_free_reachables)==0){
    stop("Error: All points removed as noise. Try setting min points to be smaller or eps and eps_t to be biggger.")
  }

  return(noise_free_reachables)
}


#' Find final equilibrium clusters
#'
#' Iterates through expanding clusters until finished, then returns the final clusters.
#'
#' @param clusters data.table of connections with columns `first`, `second` and `cluster`
#' @return A data.table of connections with columns `first`, `second` and `cluster`,  in which the cluster numbers are final.
find_equilibrium_clusters <- function(clusters){

  iter_count <- 0
  while (T){
    #some debugging output - this can be removed later.
    iter_count <- iter_count + 1

    #store a copy of before the iteration to compare to afterwards
    clusters_old <- clusters

    #iterate by finding the min cluster label of a first or of a second
    clusters <- clusters[ , .(first,cluster=min(cluster)), by = second]
    clusters <- clusters[ , .(second,cluster=min(cluster)), by = first]

    #if nothing changed in an iteration, we are done!
      if (T==all.equal(clusters_old, clusters)){

        #also check that each second only has one cluster number
        deduped_clusters <- clusters[,.(count = .N),by = list(second,cluster)]
        if (nrow(deduped_clusters)==length(unique(deduped_clusters$second))) {
          break
        }else{
          stop("Error - multiple clusters found for each core point")
        }
      }
  }
  message(iter_count," iterations to equilibrium cluster configuration")

  return(deduped_clusters)
}


#' Get Clusters from Data
#'
#' Takes a dataframe/tibble/data.table with identified columns for x, y, t and finds
#' clusters of points that are connected through being no more than `eps` apart in space,
#' no more than `eps_t` apart in time, and connected through points with a minimum of `min_number` connections.
#'
#' @param df dataframe/data.table/tibble containing trajectories to cluster.
#' @param x string, name of the x spatial column
#' @param y string, name of the y spatial column
#' @param t string, name of the timestamp column
#' @param eps Largest distance apart points can be to be directly connected
#' @param eps_t Longest time apart points can be to be directly connected
#' @param minpts Smallest number of points a point must be connected to, to
#' not be excluded as a possible source point for connecting further points to a cluster.
#' @return A data.table that looks like `df` but contains two additional columns, `point_density` and `cluster`.
#' @export
#' @examples
#' dt <- data.table(id = 1:100,
#'                 X = rnorm(100),
#'                 Y = rnorm(100),
#'                 timestamp = 1:100)
#' get_clusters_from_data(dt
#'                        ,x = "X"
#'                        ,y = "Y"
#'                        ,t = "timestamp"
#'                        ,eps = 2 # 2 metre distance threshold from other point
#'                        ,eps_t = 10 # 10 second time threshold from other point
#'                        ,minpts = 3) # minimum connected to 3 points to continue growing a cluster
get_clusters_from_data <- function(df
                                   ,x = "X"
                                   ,y = "Y"
                                   ,t = "timestamp"
                                   ,eps = 1 # metre from other point
                                   ,eps_t = 5 # seconds from other point
                                   ,minpts = 7){

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
  reachables <- find_connected_points(dt = dt,
                                      eps = eps,
                                      eps_t = eps_t)


  point_density <- reachables[,.(density = .N),by = first]

  data.table::setnames(point_density,old = c("first"),new = c("id"))


  # remove noise (clusters with less than minpts)
  reachables_without_noise <- remove_noise(reachables = reachables,
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

  #iterate until final clusters are here
  clusters_core <- find_equilibrium_clusters(clusters = clusters_core)[,.(second, cluster)]

  #deal with terminating points
  #get cluster number for every second from each core point it is connected to
  clusters_terminating2 <-
    data.table::merge.data.table(clusters_terminating,
                                 clusters_core,
                                 by.x = "first",
                                 by.y = "second",
                                 suffixes = c("1","2"),
                                 all.x = T)

  #keep only only min cluster number for each second
  clusters_terminating3 <-
    clusters_terminating2[
      clusters_terminating2[ , .I[which.min(cluster)],
                             by = second]$V1]


  clusters <- rbindlist(list(clusters_core,
                             clusters_terminating3[,.(second,cluster)]),
                        use.names = T)


  #remove seconds with two different clusters - take the first
  #(this is becasue terminating points can belong to two disconnected clusters)
  #clusters <- clusters[,.SD[1],by = second]

  #rename the clusters to be sequential starting from 1
  cluster_numbers_old <- na.omit(unique(clusters$cluster))
  cluster_numbers_new <- 1:length(cluster_numbers_old)
  names(cluster_numbers_new) <- cluster_numbers_old

  #extract just the point ids and the cluster they belong to
  point_to_clust <- unique(clusters[,.(point_id = second,
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


#' Create demo data to test clustering algorithm
#'
#' @param n number of rows in demo data
#'
#' @return A new data.table with `n` rows to feed into the clustering algorithm
#' @export
#'
#' @examples
#' demo_data(n = 1000)
demo_data <- function(n = 1000){
  dt = data.table(timestamp = 1:n,
                  X = rnorm(n),
                  Y = rnorm(n))
}
