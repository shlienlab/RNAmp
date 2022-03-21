#' Get confidence intervals
#' @export
get.ci <- function(x, interval = 0.95) {
    a <- mean(x)
    s <- sd(x)
    n <- length(x)
    error <- qnorm(1 - (1 - interval)/2) * s/sqrt(n)
    se <- s/sqrt(n)
    lower.ci <- a - error
    upper.ci <- a + error
    return(c(sdev = s, lower.ci = lower.ci, upper.ci = upper.ci, se = se))
}
