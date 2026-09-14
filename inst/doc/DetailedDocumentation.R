## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# f <- function(a) {
#   if (a > 0) x <- 1.0   # x is only assigned here...
#   return(x)             # ...but read here regardless: returns 0.0 when a <= 0
# }

## ----eval = FALSE-------------------------------------------------------------
# f <- function(a, b, c) {
#   argtypes(
#     a |> type(vec(double)),
#     b |> type(mat(double)),
#     c |> type(double)
#   )
#   # ... body ...
# }
# f_cpp <- ast2ast::translate(f)

## ----eval = FALSE-------------------------------------------------------------
# f <- function(a, b, c) {
#   argtypes(
#     a |> type(borrow_vec(double)) |> ref(),            # mutable, passed by reference
#     b |> type(borrow_mat(double)) |> ref() |> const(), # read-only matrix reference
#     c |> type(double) |> ref()                         # scalar reference (XPtr only)
#   )
#   # ... body ...
# }

## ----eval = FALSE-------------------------------------------------------------
# a <- 1L
# a <- 2.5

## ----eval = FALSE-------------------------------------------------------------
# f <- function() {
#   a <- 1L      # looks like an integer assignment...
#   a <- 2.5     # ...but a is reassigned a double two lines later
#   return(a)
# }

## ----eval = FALSE-------------------------------------------------------------
# fcpp <- ast2ast::translate(f, derivative = "forward")
# fcpp <- ast2ast::translate(f, derivative = "reverse")

## ----eval = FALSE-------------------------------------------------------------
# f <- function(y, x) {
#   jac <- matrix(0.0, length(y), length(x))
#   for (i in 1L:length(x)) {
#     seed(x, i)
# 
#     y[[1L]] <- x[[1L]] * x[[2L]]
#     y[[2L]] <- x[[1L]] + x[[2L]] * x[[2L]]
# 
#     d <- get_dot(y)
#     jac[TRUE, i] <- d
# 
#     unseed(x, i)
#   }
#   return(jac)
# }
# 
# fcpp_forward <- ast2ast::translate(f, derivative = "forward")

## ----eval = FALSE-------------------------------------------------------------
# f <- function(y, x) {
#   y[[1L]] <- x[[1L]] * x[[2L]]
#   y[[2L]] <- x[[1L]] + x[[2L]] * x[[2L]]
#   jac <- deriv(y, x)
#   return(jac)
# }
# 
# fcpp_reverse <- ast2ast::translate(f, derivative = "reverse")

## ----eval = FALSE-------------------------------------------------------------
# f <- function(a) {
#   argtypes(a |> type(int))
#   factorial <- fn(
#     argtypes(a |> type(int) |> const()),
#     return(int),
#     {
#       if (a == 1L) return(a) else return(a * factorial(a - 1L))
#     }
#   )
#   return(factorial(a))
# }
# fcpp <- ast2ast::translate(f)

## ----eval = FALSE-------------------------------------------------------------
# sq <- fn(
#   argtypes(x |> type(double) |> const()),  # const -- accepts expressions
#   return(double),
#   return(x * x)
# )
# # sq(a + b) is fine; without const() on x, only sq(a) (a bare variable) would be.

## ----eval = FALSE-------------------------------------------------------------
# f <- function(interval) {
#   argtypes(interval |> type(vec(double)))
#   g <- fn(
#     argtypes(x |> type(double)),
#     return(double),
#     {
#       return(x^2 - 4)
#     }
#   )
#   res <- uniroot(g, interval, 1e-10, 1000)
#   return(res$root)
# }
# fcpp <- ast2ast::translate(f)

## ----eval = FALSE-------------------------------------------------------------
# f <- function(x) {
#   argtypes(
#     x |> type(vec(double))
#   )
#   sq <- fn(
#     argtypes(
#       a |> type(double) |> const()
#     ),
#     return(double),
#     return(a * a)
#   )
#   return(map(sq, x))
# }
# fcpp <- ast2ast::translate(f)
# fcpp(1:5)

## ----eval = FALSE-------------------------------------------------------------
# rosen <- function(p) {
#   argtypes(
#     p |> type(vec(double))
#   )
#   loss <- fn(
#     argtypes(
#       x |> type(vec(double)) |> const()
#     ),
#     return(double),
#     {
#       a <- 1.0 - x[[1L]]
#       b <- x[[2L]] - x[[1L]] * x[[1L]]
#       return(a * a + 100.0 * b * b)
#     }
#   )
#   lo <- c(-5.0, -5.0)
#   up <- c(5.0, 5.0)
#   res <- lbfgsb(loss, p, lo, up, 100L, 1e7, 1e-8, 5L)
#   return(res$par)
# }
# fcpp <- ast2ast::translate(rosen, derivative = "reverse")
# fcpp(c(-1.2, 1.0))

## ----eval = FALSE-------------------------------------------------------------
# types_f <- function() {
#   new_type(Point, slots(x |> type(double), y |> type(double)))
# }
# 
# f <- function(p) {
#   argtypes(p |> type(Point))
#   p$x <- p$x + 1
#   return(p)
# }
# 
# fcpp <- ast2ast::translate(f, types_f = types_f)
# 
# p <- structure(list(x = 1, y = 2), class = "Point")
# fcpp(p)

## ----eval = FALSE-------------------------------------------------------------
# f <- function() {
#   dep <- c(0, 1, 0.5, 2.5, 3.5, 4.5, 4)
#   indep <- 1:7
#   evalpoints <- c(
#     0.5, 1, 1.5, 2, 2.5,
#     3, 3.5, 4, 4.5, 5,
#     5.5, 6, 6.5
#   )
#   for (i in evalpoints) {
#     print(cmr(i, indep, dep))
#   }
# }

