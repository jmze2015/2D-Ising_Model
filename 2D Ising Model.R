###########################
## ----2D Ising Model----##
## --(Jose Muñoz-Lopez)--##
###########################

## Lattice sizes chosen N = 10,20,30,40 (or anything you want)
for (N in c(20)){
  ## Parameters
  kB <- 1
  
  ## Define a lattice with random spins (+1, -1) of size N^2
  lattice <- matrix(sample(c(-1,1), N^2, replace = TRUE),ncol = N)
  
  ## Metropolis Algorithm
  Metro <- function(A, Beta){
    ## Pick site
    i <- sample(seq(1, nrow(A)),1)
    j <- sample(seq(1, ncol(A)),1)
    
    ## Random number
    r <- runif(1,0,1)
    
    ## Propose flip
    Accept <- A
    Accept[i,j] <- -1*Accept[i,j]
    
    ## Periodic Condition
    nr <- nrow(A)
    nc <- ncol(A)
    ip1 <- if (i == nr) 1 else i + 1
    im1 <- if (i == 1) nr else i - 1
    jp1 <- if (j == nc) 1 else j + 1
    jm1 <- if (j == 1) nc else j - 1
    
    ## Local Energy Change
    H_pre <- -1* A[i,j] * (A[i, jp1] + A[i, jm1] + A[ip1, j] + A[im1, j])
    H_post <- -1 * Accept[i,j] * (Accept[i, jp1] + Accept[i, jm1] + 
                                    Accept[ip1, j] + Accept[im1, j])
    Delta_H <- H_post - H_pre
    
    prob_accept <- exp(-1*Beta*Delta_H)
    ## Decision
    if(Delta_H <= 0){
      return(Accept)
    }
    else if (r <= prob_accept){
      return(Accept)
    }
    else{
      return(A)
    }
  }
  
  
  
  ## Visualization of a config.
  Viz <- function(A){
    cols<- c("white", "lightblue")
    image(
      1:nrow(A), 1:nrow(A),
      t(apply(A, 2, rev)),
      col = cols,
      breaks = c(-1, 0, 1),
      xlab = "Column",
      ylab = "Row"
    )  
  }
  
  
  
  ## Energy of a given configuration (matrix A) which is periodic
  Total_E <- function(A){
    energy <- 0
    nr <- nrow(A)
    nc <- ncol(A)
    for (i in 1:nr){
      for (j in 1:nc){
        ip1 <- if (i == nr) 1 else i + 1
        im1 <- if (i == 1) nr else i - 1
        jp1 <- if (j == nc) 1 else j + 1
        jm1 <- if (j == 1) nc else j - 1
        energy <- energy + -1* A[i,j] * (A[i, jp1] + A[i, jm1] + A[ip1, j]
                                         + A[im1, j])
      }
    }
    return(energy / 2)
  }
  
  
  
  
  Temp_seq <- seq(0.2,4,0.2)
  ## vector length 20
  Beta_seq <- 1/(kB*Temp_seq)
  
  ## Add in temp as a parameter, # of burn sweeps, # of extra sweeps (recorded)
  Run <- function(M, burn, extra, Beta){
    Energies <- c()
    sweep <- N*N
    burn_sweeps <- burn*sweep 
    for (i in 1:burn_sweeps){
      M <- Metro(M, Beta)
    }
    for (j in seq(sweep, extra*sweep)){
      M <- Metro(M, Beta)
      if (j %% sweep  == 0){
        Energies <- c(Energies, Total_E(M))
      }
    }
    u_est <- mean(Energies) * (1/sweep)
    u_var_est <- var(Energies)
    C <- kB * Beta * Beta * u_var_est
    c <- C/(N^2)
    return(list(
      final_lattice = M,
      energies = Energies,
      avg_energy = u_est,
      beta = Beta,
      specific_heat = c
    ))
  }
  
  
  avg_E <- c()
  for (s in 1:length(Temp_seq)){
    avg_E[s] <- Run(lattice, 1000,2000, Beta_seq[s])$avg_energy
  }
  
  
  ## Temperature vs Internal Energy Plot
  plot(Temp_seq, avg_E, type = "b", main = paste("Lattice of Size N=",N),
       xlab = "Temperature (K)", ylab = "Internal Energy Estimate")
  #abline(v=2.269, col="red")
}



## Requires more sweeps than energies
# spec_H <- c()
# for (s in 1:length(Temp_seq)){
#   spec_H[s] <- Run(lattice, 100,1000, Beta_seq[s])$specific_heat
# }

## Temperature vs Specific Heat Plot
# plot(Temp_seq, spec_H, type = "b", main= paste("Specific Heat (N =",N,")"),
#      xlab = "Temperature", ylab = "Specific Heat")
# 
# abline(v=2.269, col="red")















# start <- Sys.time()
# ## in between code
# end <- Sys.time()
# time_elapsed <- end - start
# print(time_elapsed)













