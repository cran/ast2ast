## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE, cache = FALSE)

## ----eval = FALSE, echo = TRUE, attr.source='.numberLines'--------------------
# bubbleSort <- function(a) {
#   size <- length(a)
#   for (i in 1:size) {
#     for (j in 1:(size - 1)) {
#       if (a[j] > a[j + 1]) {
#         temp <- a[j]
#         a[j] <- a[j + 1]
#         a[j + 1] <- temp
#       }
#     }
#   }
#   return(a)
# }

## ----eval = TRUE, echo = TRUE-------------------------------------------------
f <- function(a, b) {
  argtypes(
    a |> type(vec(double)) |> ref(),
    b |> type(vec(double)) |> ref()
  )
  c <- a + b
  return(c)
}
fcpp <- ast2ast::translate(f, output = "XPtr")

## ----eval = TRUE, echo = TRUE-------------------------------------------------
call_xptr(fcpp)

