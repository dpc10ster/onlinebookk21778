AddZetaXFwd <- function(f,n){
  # adds n arguments to a function f; returns that new function 
  t <- paste("arg <- alist(", paste(names(formals(f)), collapse="=, "), "=, ",
            paste(sapply(1:n, function(i) paste("zetaXFwd", i, "=", sep = "")), collapse = ", "),
            ")", sep= "")
  formals(f) <- eval(parse(text = t))
  return(f)
}

AddZetaYFwd <- function(f,n){
  # adds n arguments to a function f; returns that new function 
  t <- paste("arg <- alist(", paste(names(formals(f)), collapse="=, "), "=, ",
             paste(sapply(1:n, function(i) paste("zetaYFwd", i, "=", sep = "")), collapse = ", "),
             ")", sep= "")
  formals(f) <- eval(parse(text = t))
  return(f)
}

AddZetaFwd <- function(f,n){
  # adds n arguments to a function f; returns that new function 
  t <- paste("arg <- alist(", paste(names(formals(f)), collapse="=, "), "=, ",
             paste(sapply(1:n, function(i) paste("zetaFwd", i, "=", sep = "")), collapse = ", "),
             ")", sep= "")
  formals(f) <- eval(parse(text = t))
  return(f)
}

AddZetaX <- function(f,n){
  # adds n arguments to a function f; returns that new function 
  t <- paste("arg <- alist(", paste(names(formals(f)), collapse="=, "), "=, ",
             paste(sapply(1:n, function(i) paste("zetaX", i, "=", sep = "")), collapse = ", "),
             ")", sep= "")
  formals(f) <- eval(parse(text = t))
  return(f)
}

AddZetaY <- function(f,n){
  # adds n arguments to a function f; returns that new function 
  t <- paste("arg <- alist(", paste(names(formals(f)), collapse="=, "), "=, ",
             paste(sapply(1:n, function(i) paste("zetaY", i, "=", sep = "")), collapse = ", "),
             ")", sep= "")
  formals(f) <- eval(parse(text = t))
  return(f)
}