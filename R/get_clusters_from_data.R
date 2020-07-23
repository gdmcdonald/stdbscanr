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

#' One cluster iteration
#'
#' One iteration of finding connected groups of points and labelling them as the same cluster.
#' Takes a data.table of clusters and merges it with itself to find the next layer of connections,
#' and relables clusters with the minimum of their two previous cluster numbers.
#'
#' @param clusters data.table of connections with columns `first`, `second` and `cluster`
#' @return A data.table of connections with columns `first`, `second` and `cluster`,  in which the cluster numbers have been updated.
cluster_iteration <- function(clusters){

  # One row per t=(n-1) seed-cluster connection
  cluster2<-unique(clusters[,.(first,cluster)])

  # Merge with itself and take minimum cluster number
  # to get t=n seed-cluster connection
  clusters_second <- data.table::merge.data.table(clusters,
                                                  cluster2,
                                                  by.x = "second",
                                                  by.y = "first",
                                                  suffixes = c("1","2"),
                                                  all.x = T)

  clusters_third <- clusters_second[,cluster:=min(c(cluster1,
                                                    cluster2),
                                                  na.rm = T),by=first]

  return(clusters_third[,.(first, second, cluster)])
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

    #iterate giving the cluster labels of firsts to their seconds
    clusters <- cluster_iteration(clusters = clusters_old)

    #if nothing changed in an iteration, we are done!
    if (T==all.equal(clusters_old, clusters)) break
  }
  message(iter_count," iterations to equilibrium cluster configuration")

  return(clusters)
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

  # initially set cluster_id to first index
  clusters <- reachables_without_noise[,cluster:=first]

  #sort by first for faster merge
  setkey(clusters,first)

  #iterate until final clusters are here
  clusters <- find_equilibrium_clusters(clusters = clusters)

  #rename the clusters to be sequential starting from 1
  cluster_numbers_old <- na.omit(unique(clusters$cluster))
  cluster_numbers_new <- 1:length(cluster_numbers_old)
  names(cluster_numbers_new) <- cluster_numbers_old

  #extract just the point ids and the cluster they belong to
  point_to_clust <- unique(clusters[,.(point_id = second,
                                       cluster = cluster_numbers_new[as.character(cluster)])])

  temp <- data.table::merge.data.table(x = df,
                                       y = point_to_clust[,.(point_id,
                                                             cluster)],
                                       by.x = "id",
                                       by.y = "point_id",
                                       all.x = T)

  df_out <- data.table::merge.data.table(x = temp,
                                         y = point_density,
                                         by.x = "id",
                                         by.y = "id",
                                         all.x = T)

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
