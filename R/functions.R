difference <- function(x){
  res <- x[['Predicted_VCF']]-x[['VCF']]
  return(res)
}

sq_diff <- function(x){
  res <- (x[['Predicted_VCF']]-x[['VCF']])^2
  return(res)
}