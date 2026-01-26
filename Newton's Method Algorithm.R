## Newton's Method ##

## Must give input as expression in "x"
Newt <- function(expr, x0, tol = 1e-8, max_iter = 100){
  f_expr <- expr
  df_expr <- D(f_expr, "x")
  
  f <- function(x){
    eval(f_expr)
  }
  df <- function(x){
    eval(df_expr)
  }
  
  i <- 0
  x_seq <- c(x0)
  x_old <- x0
  
  for (i in 1:max_iter) {
    x_new <- x_old - f(x_old) / df(x_old)
    x_seq <- c(x_seq, x_new)
    
    if (abs(x_new - x_old) < tol) {
      break
    }
    
    x_old <- x_new
  }
  return(x_seq)
}

## Expression for Problem 2
f <- expression(x^2+x-1)
## Calling NM function
Newt(f, 1, 1e-6)
