test_that("potts gibbs sampler works", {

  isingibbs3Test=function(niter,n,m=n,J=0.2, ranSeed = 6){
    set.seed(ranSeed)
    # initialization
    x=sample(c(1,2,3),n*m,prob=rep(1/3,3),replace=TRUE)
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

  isingibbs2ToyTest=function(niter,n,m=n,J=0.2, ranSeed = 6){
    set.seed(ranSeed)
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

  pottsGibbsSamplerTest = function(numIteration, n, m=n, J=0.2, stateSpace, ranSeed =6){
    set.seed(ranSeed)
    pottsGibbsSampler(numIteration,n,m,J,stateSpace)
  }

  expect_equal(pottsGibbsSamplerTest(10,20,20,0.6,c(1,2,3)), isingibbs3Test(10,20,20,J=0.6))
  expect_equal(pottsGibbsSamplerTest(10,20,20,0.6,c(-5,5)), isingibbs2ToyTest(10,20,20,J=0.6))
})
