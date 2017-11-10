## Class to describe a SPLMM model
SPLMM = setRefClass("SPLMM",
                    fields = list(x = "matrix", # d x n matrix
                                  y = "numeric", # of length n
                                  lambda_0 = "numeric",
                                  lambda_1 = "numeric",
                                  lambda_2 = "numeric",
                                  maxit = "integer",
                                  eps_abs = "numeric",
                                  eps_rel = "numeric",
                                  rho = "numeric",
                                  init_w = "numeric",
                                  init_w_prime = "numeric",
                                  init_z = "numeric",
                                  init_eta = "numeric",
                                  center_x = "numeric",
                                  scale_x = "numeric",
                                  useGrad = "logical",
                                  verbose = "logical",
                                  stopme = "logical")
)

##### Member functions of SPLMM #####

## Initialize fields including default values
SPLMM$methods(
  initialize = function(x, y, lambda_0, lambda_1,lambda_2, ...)
  {
    d = nrow(x)
    if(nrow(x) >= ncol(x))
      stop("ncol(x) must be greater than nrow(x)")
    if(ncol(x) != length(y))
      stop("ncol(x) should be equal to length(y)")

    for(i in 1:ncol(x)){
    if((y[i] != 1) & (y[i] != -1)) stop("y values should be equal to 1 or -1");
    }

    ## how about scaling and centering?
    xx = scale(t(x))
    .self$center_x = attr(xx, "scaled:center")
    .self$scale_x = attr(xx, "scaled:scale")

    attr(xx, "scaled:center")<-NULL
    attr(xx, "scaled:scale")<-NULL
    ## Transfrom x
    xx = t(t(xx)*y)


    .self$x = as.matrix(t(xx))
    .self$y = as.numeric(y)
    .self$lambda_0 = as.numeric(lambda_0)
    .self$lambda_1 = as.numeric(lambda_1)
    .self$lambda_2 = as.numeric(lambda_2)
    .self$maxit = 1000L
    .self$eps_abs = 1e-8
    .self$eps_rel = 1e-8
    .self$rho = 1
    .self$init_w = as.numeric(rep(1/(10*nrow(x)),d))
    .self$init_w_prime =as.numeric(rep(-1/(nrow(x)),d))
    .self$init_z = as.numeric(rep(0.01,d))
    .self$init_eta = as.numeric(rep(0.01,d))

}
)

## Print off SPLMM object
SPLMM$methods(
  show_common = function()
  {
    cat(sprintf("$x: <%d x %d> matrix\n", nrow(.self$x), ncol(.self$x)))
    cat(sprintf("$y: <%d x 1> vector\n", length(.self$y)))

    fields = setdiff(names(.refClassDef@fieldClasses), c("x", "y"))
    for(field in fields)
      cat("$", field, ": ", paste(.self$field(field), collapse = " "),
          "\n", sep = "")
  },
  show = function()
  {
    cat("SPLMM model\n\n")
    show_common()
  }
)


## Specify additional parameters
SPLMM$methods(
  opts = function(maxit = 10000, eps_abs = 1e-8, eps_rel = 1e-8,
                  rho = 1, center_x = 0, scale_x = 1, useGrad = TRUE, verbose = FALSE, stopme = FALSE,
                  #init_w = rep(1,nrow(x)), init_w_prime =  rep(1,nrow(x))
                    ...)
  {
    if(maxit <= 0)
      stop("maxit should be positive")
    if(eps_abs < 0 | eps_rel < 0)
      stop("eps_abs and eps_rel should be nonnegative")
    if(rho <= 0)
      stop("rho should be positive")
#  .self$init_w = as.numeric(init_w)
#    .self$init_w_prime = as.numeric(init_w_prime)
    .self$maxit = as.integer(maxit)
    .self$eps_abs = as.numeric(eps_abs)
    .self$eps_rel = as.numeric(eps_rel)
    .self$rho = as.numeric(rho)
    .self$useGrad = as.logical(useGrad)
    .self$verbose = as.logical(verbose)
    .self$stopme = as.logical(stopme)

    ## Transfrom x, keep y and use them to transform back for predictions
    xx= t(t(x)*y)
    ## how about scaling and centering?
    xx=scale(t(x))
    .self$center_x = as.numeric(attr(xx, "scaled:center"))
    .self$scale_x = as.numeric(attr(xx, "scaled:scale"))

    invisible(.self)
  }
)


## Fit model and conduct the computing
SPLMM$methods(
  fit = function(...)
  {
require(nloptr)
require(Matrix)
# Step1 Solve for w_prime
  i<-1;
  #lowerbounds = rep(-1000,nrow(x))
  #upperbounds = rep(1000,nrow(x))
  maxeval = 100
  L= function( x1, x2, x3, x4){
    out =-sum(log(pnorm(t(x)%*%(x1+x2)/sqrt(lambda_1))))+ (norm(as.matrix(x1),"2")^2)/(2*lambda_2)
    + lambda_0 * norm(as.matrix(x3),"1")
    + t(as.numeric(x4))%*%(x2-x3)
    + rho*(norm(as.matrix(x2-x3),"2")^2)/2;
    return(as.numeric(out))}
  L_z = function(x3){
    return(lambda_0 * norm(as.matrix(x3),"1")
    + t(as.numeric(init_eta))%*%(init_w-x3)
    + rho*(norm(as.matrix(init_w-x3),"2")^2)/2)


    }


  search_w_prime=function(){
  # w_prime, w, z, eta
  opts = list("algorithm"="NLOPT_LN_BOBYQA", "xtol_rel"=1.0e-8, "maxeval" = maxeval)
  ##gradient
  G_w_prime= function(xx){
    out = - colSums(as.numeric((dnorm(t(x)%*%(xx+init_w)/sqrt(lambda_1))/pnorm(t(x)%*%(xx+init_w)/sqrt(lambda_1))))*t(x)/sqrt(lambda_1)) + xx/lambda_2;
    return(as.numeric(out))}
  eval_f1 = function(xx){
    return (L(xx,init_w,init_z,init_eta))
  }

  if(isTRUE(.self$useGrad)){
    opts = list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-18, "ftol_abs" = 1.0e-18, "maxeval" = maxeval)
    temp1 = nloptr(x0 = as.numeric(init_w_prime),eval_f = eval_f1, eval_grad_f = G_w_prime, opts = opts) #lb = lowerbounds, ub =  upperbounds,
  }else{
    temp1 = nloptr(x0 = init_w_prime,eval_f = evalf1, lb=lowerbounds, ub =  upperbounds, opts = opts)
  }

  if(isTRUE(.self$verbose)){
  print("change in w_prime")
  print(t(as.matrix(init_w_prime-as.numeric(temp1$solution))))
  }
  if(norm(as.matrix(init_w_prime-as.numeric(temp1$solution)),"1")<1.0e-8) {
    if(isTRUE(verbose)){print("No more change init_w_prime")}
    .self$stopme = TRUE
  }
    init_w_prime <<- as.numeric(temp1$solution)
  }

search_w = function(){# Step 2 Solve for w

  # Gradient
  G_w= function(xx){
    #print("grad_W")
    out = - colSums(as.numeric((dnorm(t(x)%*%(init_w_prime+xx)/sqrt(lambda_1))/pnorm(t(x)%*%(init_w_prime+xx)/sqrt(lambda_1))))*t(x)/sqrt(lambda_1))+
      t(as.numeric(init_eta))+rho*(xx-init_z)
    return(out)}

  eval_f2 = function(xx){
    return (L(init_w_prime,xx,init_z,init_eta))}
  if(isTRUE(.self$useGrad)){
    opts = list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-18,"ftol_abs" = 1.0e-18, "maxeval" = maxeval)
    temp2 = nloptr(x0 = init_w, eval_f = eval_f2, eval_grad_f = G_w, opts = opts) # lb=lowerbounds, ub =  upperbounds,
  }else{
    temp2 = nloptr(x0 = init_w, eval_f = eval_f2,  opts = opts)
  }

  .self$init_w = temp2$solution}

repeat{
  search_w()
  search_w_prime()
  if(stopme){break}

# Step 3 solve for z  #NLOPT_LN_BOBYQA","NLOPT_GN_DIRECT_L",
    opts = list("algorithm"= "NLOPT_GN_DIRECT_L_RAND" ,
                "xtol_rel"=1.0e-18, "ftol_abs" = 1.0e-18, "maxeval" = 1000000L)
    #eval_f3 = function(xx){
    #  return (L(init_w_prime,init_w,xx,init_eta))}
    temp3 = nloptr(x0 = init_z, eval_f = L_z, opts = opts)
    z_diff = init_z - temp3$solution
#    if(norm(as.matrix(L(init_w_prime,init_w,init_z,init_eta)-L(init_w_prime,init_w,temp3$solution,init_eta)), "1")<1e-8){
#      print(i)
#      print("no imrpovement to loss value")
#      break}
    .self$init_z = temp3$solution

# Step 4 update eta
    .self$init_eta = init_eta + rho*(init_w-init_z)

# Step 5 Check eps_pri and eps_dual conditions
    eps_pri = sqrt(length(init_w)) * eps_abs + eps_rel * max(norm(as.matrix(init_w), "2"),norm(as.matrix(init_z),"2"))
    eps_dual = sqrt(length(init_w)) * eps_abs + eps_rel * norm(as.matrix(init_eta), "2")


# Step 6 Check top_level_eps condition
    if(isTRUE(.self$verbose)){
    print(i)
    print(c("eps_pri","eps_dual", "||w-z||_2 ", " ||z_diff||_2"))
    print(c( eps_pri,  eps_dual,  norm(as.matrix(init_w-init_z), "2"),    norm(as.matrix(z_diff),"2") ))
    print("w_prime")
    print(init_w_prime)
    print("w")
    print(init_w)
    print(L( init_w_prime, init_w, init_z, init_eta))
    }
    if ((norm(as.matrix(init_w-init_z), "2") <= eps_pri) & (norm(as.matrix(z_diff),"2") <= eps_dual)) {break}



    if(i<10){
    if(norm(as.matrix(init_w-init_z), "2")>10*norm(as.matrix(z_diff),"2")){rho<<-2*rho;}
    if(10*norm(as.matrix(init_w-init_z), "2")<norm(as.matrix(z_diff),"2")){rho<<-rho/2;}
    }


    i <- i+1;
  if(i >= maxit){    print("Max iteration reached with i="); print(i)  ;  break;}
}


    return("sufficient convergence")
}
)


## Plot SPLMM weights object
SPLMM$methods(
  plotCoefs = function(...)
  {
    require(ggplot2)
    coefs = as.numeric(.self$init_w_prime+.self$init_w)
    dat = data.frame(Index = seq_along(coefs),
                     Coefficients = coefs)
    g = ggplot(dat, aes(x = Index, y = Coefficients)) +
      geom_segment(aes(xend = Index, yend = 0))
    print(g)
    invisible(g)
  }
)


SPLMM$methods(
  weights= function(){return (as.numeric(.self$init_w))}
)

SPLMM$methods(
  predict= function(xx, epsil, eps_thres){
    n = ncol(xx)
    w =init_w  # init_w_prime +
    for(i in 1:length(w)){if (abs(w[i])<epsil){w[i]=0}}
    x_=(xx-center_x)/scale_x
    x_ = t(x_)%*%w
    y_<-NULL;
    for(i in 1:n){if (x_[i]>eps_thres){y_<-cbind(y_,1)}else{y_<-cbind(y_,-1)}}

      return(as.numeric(y_))},

  accuracy=function(xx, yy, EpsThre = 1/(10*nrow(x)), Eps= 1/(10*nrow(x))){
    if(ncol(xx) != length(yy))
      stop("ncol(x) should be equal to length(y)")
    for(i in 1:length(yy)){
      if((y[i] != 1) & (y[i] != -1)) stop("y values should be equal to 1 or -1");

    }
    y_hat = predict(xx, EpsThre, Eps)
    return(1 - 0.5*norm(as.matrix(y_hat-yy),"1")/norm(as.matrix(y_hat),"1"))
  }
)

splmm = function(x, y, ...)
{
  SPLMM(x, y, ...)
}

