## code to prepare sample dataset

set.seed(1234)

# number of cells to simulate
n.cells <- 500

# number of genes to simulation
n.genes <- 5


# creating data from the Poisson, NB, ZIP, ZINB distributions

# Poisson Simulation
# initialize counts
poi_counts <- array(0, c(n.genes, n.cells))

# set mean count for expressed cells
nz.mean <- 2

# run simulation

for(j in 1:n.genes){

  # assign gene-specific mean expression
  lambda <- rgamma(1, nz.mean, 1)

  for(k in 1:n.cells){

    poi_counts[j,k] <- rpois(1,lambda)
  }
}

#################
# Poisson Simulation with zero inflation

# proportions of not expressed cells
p.zero <- 0.8

# initialize counts
zip_counts <- array(0, c(n.genes, n.cells))

# set mean count for expressed cells
nz.mean <- 2

# run simulation

for(j in 1:n.genes){

  # assign gene-specific mean expression for non-zero cells
  lambda <- rgamma(1, nz.mean, 1)

  for(k in 1:n.cells){

    # decide if cell is expressed (1 = not expressed)
    z <- rbinom(1, 1, p.zero)

    # determine counts for extra zero and non-zero cells
    if(z==1){
      zip_counts[j,k] <- 0
    } else{
      zip_counts[j,k] <- rpois(1,lambda)
    }
  }
}


#################
# Negative Binomial Simulation

# initialize counts
nb_counts <- array(0, c(n.genes, n.cells))

# set mean count for expressed cells
nz.mean <- 3

# set mean for theta in NB expressed cells
nz.theta <- 4

# run simulation

for(j in 1:n.genes){

  # assign gene-specific mean expression for non-zero cells
  mu <- rgamma(1, nz.mean, 1)

  # assign gene-specific theta (overdispersion) for non-zero cells
  theta <- rgamma(1, nz.theta, 1)

  for(k in 1:n.cells){

    nb_counts[j,k] <- MASS::rnegbin(1,mu, theta)
  }
}


############################
# Negative Binomial Simulation with zero inflation

# proportions of not expressed cells
p.zero <- 0.8

# initialize counts
zinb_counts <- array(0, c(n.genes, n.cells))

# set mean count for expressed cells
nz.mean <- 3

# set mean for theta in NB expressed cells
nz.theta <- 4

# run simulation

for(j in 1:n.genes){

  # assign gene-specific mean expression for non-zero cells
  mu <- rgamma(1, nz.mean, 1)

  # assign gene-specific theta (overdispersion) for non-zero cells
  theta <- rgamma(1, nz.theta, 1)

  for(k in 1:n.cells){

    # decide if cell is expressed (1 = not expressed)
    z <- rbinom(1, 1, p.zero)

    # determine counts for extra zero and non-zero cells
    if(z==1){
      zinb_counts[j,k] <- 0
    } else{
      zinb_counts[j,k] <- MASS::rnegbin(1,mu, theta)
    }
  }
}


sample_data <- rbind(poi_counts, nb_counts, zip_counts, zinb_counts)

# shuffle the order of rows in the matrix
rand <- sample(nrow(sample_data))
sample_data <- sample_data[rand,]

# assign gene names
gene_names <- paste0("Gene", sep = "_", seq(1:20))
rownames(sample_data) <- gene_names

# assign colnames
sample_names <- paste0("Cell", sep = "_", seq(1:500))
colnames(sample_data) <- sample_names

# add infromation on cell-types
cell_type <- sample(1:5, size = 500, replace = TRUE)
covar_info <- data.frame(cell_type, row.names = colnames(sample_data))

# add information on the library sizes
sample_lib.size <- sample(350:900, size = 500, replace = FALSE)

scData <- list(counts=sample_data, covariates=covar_info, lib_size=sample_lib.size)

usethis::use_data(scData, compress = "xz")
