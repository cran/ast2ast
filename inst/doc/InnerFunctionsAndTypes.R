## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## -----------------------------------------------------------------------------
# types_f <- function() {
#   new_type(Point, slots(x |> type(double), y |> type(double)))
# }

## -----------------------------------------------------------------------------
# f <- function(p, q) {
#   argtypes(
#     p |> type(Point),
#     q |> type(Point)
#   )
#   squared_dist <- fn(
#     argtypes(
#       a |> type(Point) |> const() |> ref(),
#       b |> type(Point) |> const() |> ref()
#     ),
#     return(double),
#     {
#       dx <- a$x - b$x
#       dy <- a$y - b$y
#       return(dx * dx + dy * dy)
#     }
#   )
#   return(squared_dist(p, q))
# }

## -----------------------------------------------------------------------------
# fcpp <- ast2ast::translate(f, types_f = types_f)
# 
# p <- structure(list(x = 0, y = 0), class = "Point")
# q <- structure(list(x = 3, y = 4), class = "Point")
# fcpp(p, q) # 25

