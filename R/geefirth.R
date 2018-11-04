#' @title  Fitting Firth-type GEE for  correlated binary data with separation or near-to-separation
#'
#' @description geefirth fits GEE with Firth-type penalization to provide
#' bias-corrected finite  estimate of the regression coefficient  in case of separation or near-to-separation. In addition,
#' it provides bias-correctedsandwich estimate of the standard error of the penalized GEE estimate.
#'
#' @param formula Similar as geeglm
#'
#'
#' @param id cluster variable and requirs to be sorted
#'
#' @param corstr working correlation structure, "independence" (default), "ar1", "unstr", "exchangeable".
#'
#' @return NULL
#' @author Momenul Haque Mondol \email{mmondol@isrt.ac.bd}, M. Shafiqur Rahman \email{
#' shafiq@isrt.ac.bd}
#'
#'
#'
#' @examples
#'
#' # loading data
#' data(geefirth_data)
#'
#' # Fitting GEE for quasi-separated data
#'
#' geefirth(y ~ x + obstime, id=quasi_sep$id, data=quasi_sep, corstr = "exchangeable");
#'
#' # Fitting GEE for near to quasi-separated data
#'
#' geefirth(y ~ x + obstime, id=near_to_sep$id, data=near_to_sep, corstr = "exchangeable");
#'
#' @export


geefirth <- function(formula = formula(data), id = id, data = parent.frame(), corstr="independence"){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <-
    m$corstr <- m$Mv <- m$silent <- m$contrasts <-
    m$family <- m$scale.fix <- m$scale.value <- m$v4.4compat <- NULL
  if(is.null(m$id)) m$id <- as.name("id")
  if(!is.null(m$na.action) && m$na.action != "na.omit") {
    warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
    m$na.action <- as.name("na.omit")
  }
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m, contrasts)
  xx <- model.matrix(Terms, m, contrasts)
  xx <- as.data.frame(xx)
  id <- model.extract(m, id)
  if (dim(y)[2]==2) {
    N <- as.vector(y %*% c(1,1))
    y <- y[,1]
  }
  else {
    if (dim(y)[2]>2)
      stop("Only binomial response matrices (2 columns)")
  }

  if(is.null(id)) {
    stop("Id variable not found")
  }

  xnames <- dimnames(xx)[[2]]
  if(is.null(xnames)) {
    xnames <- paste("x", 1:p, sep = "")
    dimnames(xx) <- list(NULL, xnames)
  }

  k <- pp <- ncol(xx)
   ##########
  tx <- t(xx)
  tol<-2
  ####
  beta <- rep(0, k)
  while(tol > .0001){
    z <- as.matrix(xx)%*%beta
    pi <- as.vector(1/(1 + exp( - z )))
    W <- diag(pi*(1-pi))
    W05 <- W^.5
    XW2 <- (tx %*% W05)    #### X' (W ^ 1/2)
    XWXi <- solve(tx %*% W %*% as.matrix(xx))
    Hmat <- t(XW2) %*% XWXi %*% XW2
    UU <- y - pi + diag(diag(Hmat)) %*% (.5-pi)
    U.star <- tx %*% UU
    delta <- as.vector(XWXi %*% U.star)
    beta <- beta + delta
    tol <- max(abs(delta))
  }

  ############
  xx <- split(xx, id)
  y <- split(y, id)
  nc <- length(y)
  bet <- e <-  V <- DD <- part <- Z <- S1 <- S2 <- S3 <- FF <- sec.part <- list()
  p2 <- NULL
  del <- 10
  while(del > .0001){
    z <- lapply(xx, function(a) as.matrix(a) %*% beta)
    a <- lapply(z, function(a) exp(a))
    mu <- lapply(a, function(a) a/(1+a))
    var.mu <- lapply(mu, function(a) a * (1-a))
    W <- lapply(var.mu, function(x) diag(as.vector(x)))
    for(i in 1:nc){
      e[[i]] <- (y[[i]]-mu[[i]])/sqrt(var.mu[[i]])
    }

    if(corstr == "exchangeable"){
      cl.size <- sapply(e, function(x) length(x))
      frame <- lapply(e, function(x) data.frame(a=x, b=x))
      matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
      matr.sum <- lapply(matr, function(x){
        p1 <- as.matrix(x)
        ni <- nrow(p1)
        s <- sum(p1)-sum(diag(p1))
        mean.s <- s/(ni*(ni-1))
        return(mean.s)

      })
      alpha <- mean(unlist(matr.sum))


      R <- lapply(cl.size, function(x) xch(x, alpha))
    }
    if(corstr == "ar1"){
      e1 <- lapply(e, function(x) x[-length(x)])
      e2 <- lapply(e, function(x) x[-1])
      ni <- sapply(e, function(x) length(x))
      temp <- list()
      for(i in 1:nc){
        temp[[i]] <- (1/(ni[i]-1)) * sum(e1[[i]] * e2[[i]])
      }
      alpha <- Reduce("+", temp)/nc

      R <- lapply(ni, function(x) ar1(x, alpha))
    }
    if(corstr == "unstr"){
      ni <- sapply(e, function(x) length(x))
      frame <- lapply(e, function(x) data.frame(a=x, b=x))
      matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
      new.mat <- sapply(ni, function(x) max(ni)-x)
      matr.s <- list()
      for(i in 1:nc){
        matr.s[[i]] <- rbind(cbind(as.matrix(matr[[i]]), matrix(0, nrow =  ni[i], ncol = new.mat[i])),
                             matrix(0, nrow =  new.mat[i], ncol = max(ni)))
      }
      alpha <- Reduce("+", matr.s)/nc
      R <- lapply(ni, function(x) alpha[1:x, 1:x])
    }
    if(corstr=="independence"){
      ni <- sapply(e, function(x) length(x))
      R <- lapply(ni, function(x) diag(x))
    }

    phi <- mean(unlist(lapply(e, function(x) mean(x^2))))
    for(i in 1:nc){
      V[[i]] <- phi * (W[[i]])^(1/2) %*% R[[i]] %*% (W[[i]])^(1/2)
    }
    for(i in 1:nc){
      DD[[i]] <- t(as.matrix(xx[[i]])) %*% W[[i]]

    }
    for(i in 1:nc){
      xi <- as.matrix(xx[[i]]); txi <- t(xi); W12 <- W[[i]]^(1/2)
      part[[i]] <- txi %*% W12 %*% ginv(R[[i]]) %*% W12 %*% xi
    }
    I <- Reduce("+", part)/phi
    Q <- lapply(mu, function(x) diag(0.5-as.vector(x)))
    for(i in 1:nc){
      Z[[i]] <- list()
      xi <- as.matrix(xx[[i]])
      pi <- ncol(xi)
      for(j in 1:pi){
        temp <- as.vector(xi[,j])
        Z[[i]][[j]] <- diag(temp)
      }
    }

    for(i in 1:nc){
      W12 <- W[[i]]^(1/2)
      xi <- as.matrix(xx[[i]])
      txi <- t(xi)
      Qi <- Q[[i]]
      S1[[i]] <- list()
      for(j in 1:pp){
        S1[[i]][[j]] <- txi %*% W12 %*% ginv(R[[i]]) %*% W12 %*% Qi %*% Z[[i]][[j]] %*% xi
      }

    }
    for(j in 1:pp){
      FF[[j]] <- matrix(0, pp, pp)
      for(i in 1:nc){
        FF[[j]] <- FF[[j]] + S1[[i]][[j]]
      }
      FF[[j]] <- 2 * FF[[j]]/phi
    }
    for(i in 1:pp){
      p2[i] <-  sum(diag(ginv(I) %*% FF[[i]]))
    }
    for(i in 1:nc){
      xi <- as.matrix(xx[[i]]);  txi <- t(xi);  Wi <- W[[i]];  yi <- y[[i]]; mui <- mu[[i]]
      sec.part[[i]] <- txi %*% Wi %*% ginv(V[[i]]) %*% (yi - mui)

    }
    Ustar <- Reduce("+", sec.part) + (0.5 * p2)
    tt <- ginv(I) %*% Ustar
    del <- max(abs(tt))
    beta <- beta + tt
  }

  ##########################################
  temp <- e <- V <- DD <- sec.part <- tempS <- list()
  ni <- NULL
  z <- lapply(xx, function(a) as.matrix(a) %*% beta)
  mu <- lapply(z, function(a) exp(a)/(1+exp(a)))
  var.mu <- lapply(mu, function(a) a * (1-a))
  W <- lapply(var.mu, function(x) diag(as.vector(x)))
  for(i in 1:nc){
    e[[i]] <- (y[[i]]-mu[[i]])/sqrt(var.mu[[i]])
  }
  if(corstr == "exchangeable"){
    cl.size <- sapply(e, function(x) length(x))
    frame <- lapply(e, function(x) data.frame(a=x, b=x))
    ni <-  cl.size
    matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
    matr.sum <- lapply(matr, function(x){
      p1 <- as.matrix(x)
      nni <- nrow(p1)
      s <- sum(p1)-sum(diag(p1))
      mean.s <- s/(nni*(nni-1))
      return(mean.s)

    })
    alpha <- mean(unlist(matr.sum))
    R <- lapply(cl.size, function(x) xch(x, alpha))
  }
  if(corstr == "ar1"){
    e1 <- lapply(e, function(x) x[-length(x)])
    e2 <- lapply(e, function(x) x[-1])
    ni <- sapply(e, function(x) length(x))
    temp <- list()
    for(i in 1:nc){
      temp[[i]] <- (1/(ni[i]-1)) * sum(e1[[i]] * e2[[i]])
    }
    alpha <- Reduce("+", temp)/nc
    R <- lapply(ni, function(x) ar1(x, alpha))
  }
  if(corstr == "unstr"){
    ni <- sapply(e, function(x) length(x))
    frame <- lapply(e, function(x) data.frame(a=x, b=x))
    matr <- lapply(frame, function(x) with(x, sapply(a, function(x) x*b)))
    new.mat <- sapply(ni, function(x) max(ni)-x)
    matr.s <- list()
    for(i in 1:nc){
      matr.s[[i]] <- rbind(cbind(as.matrix(matr[[i]]), matrix(0, nrow =  ni[i], ncol = new.mat[i])),
                           matrix(0, nrow =  new.mat[i], ncol = max(ni)))
    }
    alpha <- Reduce("+", matr.s)/nc
    R <- lapply(ni, function(x) alpha[1:x, 1:x])
  }
  if(corstr=="independence"){
    ni <- sapply(e, function(x) length(x))
    R <- lapply(ni, function(x) diag(x))
  }

  phi <- mean(unlist(lapply(e, function(x) mean(x^2))))

  for(i in 1:nc){
    V[[i]] <- phi * (W[[i]])^(1/2) %*% R[[i]] %*% (W[[i]])^(1/2)
  }
  for(i in 1:nc){
    DD[[i]] <- t(as.matrix(xx[[i]])) %*% W[[i]]

  }
  for(i in 1:nc){
    temp[[i]] <- (DD[[i]]) %*% ginv(V[[i]]) %*% t(DD[[i]])
  }
  Vm <- ginv(Reduce("+",temp))
  for(i in 1:nc){
    xi <- as.matrix(xx[[i]]); txi <- t(xi); Wi <- W[[i]]; yi <- y[[i]]
    mui <- mu[[i]]
    sec.part[[i]] <- txi %*% Wi %*% ginv(V[[i]]) %*% (yi - mui)

  }
  d.bar <- Reduce("+", sec.part)/nc
  for(i in 1:nc){
    t1 <- (sec.part[[i]] - d.bar)
    tempS[[i]] <- t1 %*% t(t1)
  }
  N <- sum(ni)
  B <- Reduce("+", tempS) * (nc/(nc-1)) * ((N-1)/(N-pp))
  fi <- max(1, sum(diag(Vm %*% B))/nrow(B))
  deln <- min(0.5, (nrow(B)/(nc - nrow(B))))
  VsM <- (Vm %*% B %*% Vm) + (fi * deln * Vm)
  diag(VsM) <- abs(diag(VsM))
  VsMod=sqrt(diag(VsM))
  dimen <- sapply(R, function(x) dim(x)[2])
  ##########################################
  t.swM=round(beta/VsMod, 4)
  p.value=round(pt(t.swM, N-pp, lower.tail = F), 4)
  est.swM <- data.frame(coefficients=beta, std.err=VsMod, Wald=(t.swM)^2, p.val=p.value)
  row.names(est.swM) <- xnames

  fit <- list()
  attr(fit, "class") <- c("gee")
  fit$call <- call
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$contrasts <- attr(x, "contrasts")
  fit$coefficients <- est.swM
  fit$correlation <- R[dimen==max(dimen)][[1]]
  fit

}





ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}


# This xch function is from binarySimCLF package
xch <- function(n, rho){
  if (n <= 0){
    stop('n must be at least 1');
  }
  else if (n==1){
    r <- 1;
  }
  else if (n > 1){
    c1 = rep(rho,n-1);
    r = toeplitz(c(1, c1)) ;
  }
  return(r);
}


# This ar1 function is from binarySimCLF package
ar1 <- function(n, rho){
  if (n <= 0){
    stop('n must be at least 1');
  }
  if (n <= 2){
    return(xch(n,rho));
  }
  else{
    line1 <- c(1,rho^(1:(n-1)));
    return(toeplitz(line1));
  }
}

