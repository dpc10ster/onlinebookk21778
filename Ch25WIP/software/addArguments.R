addArguments <- function(f,n){
  # adds n arguments to a function f; returns that new function 
  t <- paste("arg <- alist(", paste(names(formals(f)), collapse="=, "), "=, ",
            paste(sapply(1:n, function(i) paste("zetaFwd", i, "=", sep = "")), collapse = ", "),
            ")", sep= "")
  formals(f) <- eval(parse(text = t))
  return(f)
}