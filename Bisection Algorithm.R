## Bisection Algorithm ##



G <- function(x){
  return(x^2 +x -1)
}

Bis_Alg <- function(H, a, b, max_steps){
  #H is the function and [a,b] is the interval.
  if(a < b){
    if(H(a)*H(b)<0){
      #init counter and interval
      i <- 0
      A <- a
      B <- b
      
      a_seq <- c()
      b_seq <- c()
      x_seq <- c()
      f_seq <- c()
      while (i <= max_steps){
        x <- (A + B)/2
        f <- H(x)
        
        # Saving the sequences
        a_seq <- c(a_seq, A)
        b_seq <- c(b_seq, B)
        x_seq <- c(x_seq, x)
        f_seq <- c(f_seq, f)
        
        if (H(x)*H(B) < 0){#right push case
          A <- x
          B <- B
        }
        else if (H(x)*H(B) == 0){#rare root finding case
          print(paste("The root is ",x))
        }
        else{#left push case
          A <- A
          B <- x
        }
        #counter
        i <- i+1
      }
      return(list(
        A = a_seq, B = b_seq, X = x_seq, f = f_seq, Estimate = x
      ))
    }
    else {print("Bisection Algorithm Fail: Try a different interval!")}
  }
  else{print("Error: Improper Interval")}
##end of function 
}

round(Bis_Alg(G,0,1,9)$f,4)










