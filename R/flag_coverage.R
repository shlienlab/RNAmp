#' Apply threshold cutoff function
#' @param depth Value of interest
#' @param threshold Threshold to apply cutoff
#' @param test Type of test to perform.
#' @export
flag_coverage <- function(depth, threshold, test = c("greater eq", "greater", "less eq", "less")) {
    
    stopifnot(test[1] %in% c("greater eq", "greater", "less eq", "less"))
    
    if (test[1] == "greater eq") {
        return(factor(ifelse(depth >= threshold, 0, 1), levels = c(0, 1), labels = c("pass", "fail")))
    }
    if (test[1] == "greater ") {
        return(factor(ifelse(depth > threshold, 0, 1), levels = c(0, 1), labels = c("pass", "fail")))
    }
    if (test[1] == "less eq") {
        return(factor(ifelse(depth <= threshold, 0, 1), levels = c(0, 1), labels = c("pass", "fail")))
    }
    if (test[1] == "less") {
        return(factor(ifelse(depth < threshold, 0, 1), levels = c(0, 1), labels = c("pass", "fail")))
    }
}
