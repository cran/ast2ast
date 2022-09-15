## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  f <- function(a) {
#    d_db = 1
#    ret <- a + 2 + d_db
#    return(ret)
#  }

## ---- eval = TRUE-------------------------------------------------------------
f <- function(a) {
  b <- a + 2
  return(b)
}  
library(ast2ast)
f_cpp <- translate(f, output = "XPtr", types_of_args = "sexp", return_type = "sexp")

## ---- eval = TRUE-------------------------------------------------------------
call_package(f_cpp)

## ---- eval = TRUE-------------------------------------------------------------
trash <- fct()

## ---- eval = TRUE-------------------------------------------------------------
trash <- fct()

## ---- eval = TRUE-------------------------------------------------------------
trash <- fct()

## ---- eval = TRUE-------------------------------------------------------------
f <- function(a) {
  a <- a + 2
}

library(ast2ast)
fa2a <- translate(f, reference = TRUE, output = "XPtr", types_of_args = "sexp", return_type = "void")
trash <- call_package(fa2a)

## ---- eval = TRUE-------------------------------------------------------------
library(RcppXPtrUtils)
library(Rcpp)
library(ast2ast)
library(r2sundials)
library(microbenchmark)

# R version
# ==============================================================================
ti <- seq(0, 5, length.out=101)
p <- list(a = 2)
p <- c(nu = 2, a = 1)
y0 <- 0
frhs <- function(t, y, p, psens) {
  -p["nu"]*(y-p["a"])
} 

res1 <- r2cvodes(y0, ti,
                    frhs, param = p)



# external pointer version
# ==============================================================================
ptr_exp=cppXPtr(code='
int rhs_exp(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
  NumericVector p(param);
  ydot[0] = -p["a"]*(y[0]-1);
  return(CV_SUCCESS);
}
', depends=c("RcppArmadillo","r2sundials","rmumps"),
includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
# For ease of use in C++, we convert param to a numeric vector instead of a list.

pv = c(a= 2)
# new call to r2cvodes() with XPtr pointer ptr_exp.
res3=r2sundials::r2cvodes(y0, ti, ptr_exp, param=pv)

## ---- eval = TRUE-------------------------------------------------------------
# ast2ast XPtr version
# ==============================================================================
ti <- seq(0, 5, length.out=101)
y0 <- 0

library(ast2ast)
ode <- function(y, ydot) {
  nu_db <- 2
  a_db <- 1
  ydot[1] <- -nu_db*(y[1] - a_db)
}
pointer_to_ode <- translate(ode,
                            reference = TRUE, output = "XPtr",
                            types_of_args = "sexp",
                            return_type = "void", verbose = TRUE)
res4 <- solve_ode(pointer_to_ode,
                      ti, y0)
res4 <- as.vector(res4)




# Rfunction-ast2ast version
# ==============================================================================
ode <- function(t, y, p, psens) {
  nu_db <- 2
  a_db <- 1
  return(-nu_db*(y - a_db))
}

odecpp <- translate(ode)

ti <- seq(0, 5, length.out=101)
y0 <- 0
res2 <- r2cvodes(y0, ti,
                    odecpp, param = p)


# check results and conduct benchmnark
# ==============================================================================
df <- data.frame(time = rep(ti, 4),
                 method = c(rep("R", 101), rep("R ast2ast", 101),
                            rep("XPtr", 101), rep("XPtr ast2ast", 101) ), 
                 y = c(res1, res2, res3, res4))

library(ggplot2)
ggplot() +
  geom_point(data = df, aes(x = time, y = y,
                            colour = method, shape = method,
                            group = method, size = method) ) +
  scale_colour_manual(values = c("#FFDB6D", "#D16103", "#52854C", "#293352") ) +
  scale_size_manual(values = c(6, 4, 2, 1.5)) +
  scale_shape_manual(values=c(19,17,4, 3)) +
  labs(colour = "")


r <- microbenchmark::microbenchmark(
                               r2cvodes(y0, ti,
                                        frhs, param = p),
                               r2cvodes(y0, ti,
                                        odecpp, param = p),
                               r2sundials::r2cvodes(y0, ti, ptr_exp, param=pv),
                               solve_ode(pointer_to_ode,
                                                     ti, y0))
boxplot(r, names = c("pure R", "ast2ast R", "C++", "ast2ast XPtr"))


