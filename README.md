
# stdbscanr

<!-- badges: start -->
<!-- badges: end -->

The goal of stdbscanr is to find clusters in trajectory data in x, y, and t, such that
1. Points are directly connected to other points within `eps` distance in the x,y plane and within `eps_t` in the t axis
2. Points are not connnected through points which have less than minpts connections themselves (these are terminating points)
3. To do so in a way that is efficient in terms of memory and time.


## Installation

You can install the released version of stdbscanr from my github repo with:

``` r
install.packages("devtools")
devtools::install_github("gdmcdonald/stdbscanr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(stdbscanr)

# make some trajectory data
dt <- data.table(X = rnorm(1000),
                 Y = rnorm(1000),
                 timestamp = 1:1000)

# find clusters and density in it
out <- get_clusters_from_data(dt
                             ,x = "X"
                             ,y = "Y"
                             ,t = "timestamp"
                             ,eps = 2 # 2 metre distance threshold from other point
                             ,eps_t = 10 # 10 second time threshold from other point
                             ,minpts = 3) # minimum connected to 3 points to continue growing a cluster
```

