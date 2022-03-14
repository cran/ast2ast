## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## ---- eval = TRUE-------------------------------------------------------------
f <- function(a) {
  b <- a + 2
  return(b)
}  
library(ast2ast)
f_cpp <- translate(f)

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
fa2a <- translate(f, reference = TRUE)
trash <- call_package(fa2a)

## ---- eval = TRUE-------------------------------------------------------------
library(Rcpp)
library(ast2ast)
library(r2sundials)
library(RcppXPtrUtils)
library(microbenchmark)

# R version
ti <- seq(0, 5, length.out=101)
p <- list(a = 2)
p <- c(nu = 2, a = 1)
y0 <- 0
frhs <- function(t, y, p, psens) {
  -p["nu"]*(y-p["a"])
} 

res_exp <- r2cvodes(y0, ti,
                    frhs, param = p)
attributes(res_exp) <- NULL

# External pointer
ptr_exp <- cppXPtr(code = '
int rhs_exp(double t, const vec &y,
            vec &ydot,
            RObject &param,
            NumericVector &psens) {
            
  double a = 1;
  double nu = 2;
  ydot[0] = -nu*(y[0] - a);
  return(CV_SUCCESS);
}
',
depends=c("RcppArmadillo",
          "r2sundials","rmumps"),
includes="using namespace arma;\n
          #include <r2sundials.h>",
cacheDir="lib", verbose=FALSE)

pv <- c(a = 1)
res_exp2 <- r2cvodes(y0, ti,
                     ptr_exp, param = pv)
attributes(res_exp2) <- NULL

## ---- eval = TRUE-------------------------------------------------------------
# ast2ast version
ti <- seq(0, 5, length.out=101)
y0 <- 0

library(ast2ast)
ode <- function(y, ydot) {
  nu <- 2
  a <- 1
  ydot[1] <- -nu*(y[1] - a)
}
pointer_to_ode <- translate(ode,
                          reference = TRUE)
res_exp3 <- solve_ode(pointer_to_ode,
                      ti, y0)
attributes(res_exp3) <- NULL

stopifnot(identical(res_exp,
                res_exp2,
                res_exp3))
out <- microbenchmark(
  r2cvodes(y0, ti,
          frhs, param = p),
  r2cvodes(y0, ti,
          ptr_exp, param = pv),
  solve_ode(pointer_to_ode,
          ti, y0))

boxplot(out, names=c("R", "C++", "ast2ast"))

## ---- eval = TRUE-------------------------------------------------------------
#states
path <- system.file("examples",
                    package = "paropt")
states <- read.table(paste(
                      path,
                      "/states_LV.txt",
                      sep = ""),
                     header = T)

# parameter
lb <- data.frame(time = 0, 
                 a = 0.8,
                 b = 0.3,
                 c = 0.09,
                 d = 0.09)
ub <- data.frame(time = 0,
                 a = 1.3,
                 b = 0.7,
                 c = 0.4,
                 d = 0.7)

suppressMessages(library(paropt))
set.seed(1)
start_time <- Sys.time()
df_cpp <- optimizer_pointer(
  integration_times = states$time,
  ode_sys = test_optimization(),
  relative_tolerance = 1e-6,
  absolute_tolerances = c(1e-8, 1e-8),
  lower = lb, upper = ub, states = states,
  npop = 40, ngen = 1000, error = 0.0001,
  solvertype = "bdf")
end_time <- Sys.time()
cpp_time <- end_time - start_time

# ast2ast with at and _db
ode <- function(params, states) {
  a_db = at(params, 1)
  b_db = at(params, 2)
  c_db = at(params, 3)
  d_db = at(params, 4)
  n1_db = at(states, 1)
  n2_db = at(states, 2)
  at(states, 1) = n1_db*c_db*n2_db -
    n1_db*d_db;
  at(states, 2) = n2_db*a_db - 
    n2_db*b_db*n1_db;
}

pointer_to_ode <- ast2ast::translate(
                        ode,
                        reference = TRUE)
set.seed(1)
start_time <- Sys.time()
df_ast2ast <- optimize_paropt(pointer_to_ode,
                              states$time,
                              lb, ub, states)
end_time <- Sys.time()
a2a_time <- end_time - start_time

stopifnot(identical(df_cpp[[8]],
                    df_ast2ast[[8]]) )

times <- data.frame(cpp = cpp_time,
                    ast2ast = a2a_time)
kableExtra::kbl(times)

