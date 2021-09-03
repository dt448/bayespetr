############### Ising Sampling Functions #####
#' @export
pottsGibbsSampler = function(numIteration, n, m=n, J=0.2, stateSpace){
  numStates = length(stateSpace)
  probs = rep(1/numStates, numStates)

  # initialization
  config = sample(stateSpace, n*m, prob = probs, replace=TRUE)
  config = matrix(config, n, m)

  # Gibbs runs
  for (i in 1:numIteration){
    sampl1=sample(1:n)
    sampl2=sample(1:m)
    for (k in 1:n){
      for (l in 1:m){

        nN = rep(NA,numStates)
        for(neighType in 1:numStates){
          nN[neighType] = xneig4(config,sampl1[k],sampl2[l],stateSpace[neighType])
        }

        config[sampl1[k],sampl2[l]] = sample(stateSpace, 1, prob=exp(J*nN))

      }}}
  config
}

#' @export
xneig4=function(x,a,b,col){
  #dimensions of the matrix
  n=dim(x)[1];m=dim(x)[2]


  nei=c(x[a-1,b]==col,x[a,b-1]==col)

  if (a!=n)
    nei=c(nei,x[a+1,b]==col)
  if (b!=m)
    nei=c(nei,x[a,b+1]==col)

  sum(nei)
}

isingibbs3=function(niter,n,m=n,J=0.2){
  # initialization
  x=sample(c(1,2,3),n*m,prob=c(0.33,0.33,0.34),replace=TRUE)
  x=matrix(x,n,m)

  # Gibbs runs
  for (i in 1:niter){
    sampl1=sample(1:n)
    sampl2=sample(1:m)
    for (k in 1:n){
      for (l in 1:m){
        n1=xneig4(x,sampl1[k],sampl2[l],1)
        n2=xneig4(x,sampl1[k],sampl2[l],2)
        n3=xneig4(x,sampl1[k],sampl2[l],3)
        x[sampl1[k],sampl2[l]] = sample(c(1,2,3),1,prob=exp(J*c(n1,n2,n3)))
      }}}
  x
}

isingibbs2Toy=function(niter,n,m=n,J=0.2){
  # initialization
  x=sample(c(-5,5),n*m,prob=c(0.5,0.5),replace=TRUE)
  x=matrix(x,n,m)

  # Gibbs runs
  for (i in 1:niter){
    sampl1=sample(1:n)
    sampl2=sample(1:m)
    for (k in 1:n){
      for (l in 1:m){
        n1=xneig4(x,sampl1[k],sampl2[l],-5)
        n2=xneig4(x,sampl1[k],sampl2[l],5)
        x[sampl1[k],sampl2[l]] = sample(c(-5,5),1,prob=exp(J*c(n1,n2)))
      }}}
  x
}

isingibbs2Pet=function(niter,n,m=n,J=0.2){
  # initialization
  x=sample(c(1,2),n*m,prob=c(0.5,0.5),replace=TRUE)
  x=matrix(x,n,m)

  # Gibbs runs
  for (i in 1:niter){
    sampl1=sample(1:n)
    sampl2=sample(1:m)
    for (k in 1:n){
      for (l in 1:m){
        n1=xneig4(x,sampl1[k],sampl2[l],1)
        n2=xneig4(x,sampl1[k],sampl2[l],2)
        x[sampl1[k],sampl2[l]] = sample(c(1,2),1,prob=exp(J*c(n1,n2)))
      }}}
  x
}

isingPrior = function(i,j,value,matrix,J){
  #Unnormalized
  countNeigh = xneig4(matrix,i,j,value)
  exp(J*countNeigh)

}

pottsPrior = function(config,J,n,m){
  #Unnormalized
  countNeigh= 0
  for(i in 1:n){
    for(j in 1:m){
      countNeigh = countNeigh + xneig4(config,i,j,config[i,j])
    }
  }
  exp(J*countNeigh)
}

pottsPriorSimplified = function(i,j,value,config,J){
  #Unnormalized
  countNeigh = xneig4(config,i,j,value)
  exp(2*J*countNeigh)

}
