# ###############################################################################
# ##############          Test Optim / nlm                          #############
# ###############################################################################
#
# # Specify "your" variables for the model (and lower & upper bounds for optim)
# N <- 30
# Time <- 5
#
# beta <- c(0.5, 0.5); K <- length(beta);
# delta=c(0.5, 0.5); R <- length(delta)
# sigma_u = 10; sigma_v = 10
# mu <- 0.5
# lower <- c(0.05,0.05, rep(-Inf,K), rep(-Inf,R)); upper <-  c(1000,1000,1000,1000,1000,1000)
#
# # Generates data for the optimization & does it "Iteration" times
# Iteration <- 1
# estimate <- matrix(c(rep(NA, (2+K+R)*Iteration)), nrow= 2 + K + R, ncol=Iteration)
# for(i in 1:Iteration){
#   # Generate dataset for optim based on inputs above
#   set.seed(sample(99999,1))
#   data <- SFM.generate(N = N, Time = Time, beta = beta, delta = delta, mu_u = mu, sigma_u = sigma_u, sigma_v = sigma_v)
#   x    <- as.matrix(data[, 1:length(beta)]); y <- as.matrix(data[, (length(beta)+1)])
#   z    <- as.matrix(data[, (length(beta)+2):(length(beta)+2+length(delta)-1)])
#   test <- nlminb(objective = SFM.within, lower=c(0.001,0.001,rep(-Inf, R+K)), start = c(sigma_u,sigma_v,beta,delta), mu = mu, Time = Time, N = N, xv = x, y = y, z = z, optim = T)
#
#   estimate[, i] <- test$par
# }
#   typeof(data)
# x <- as.matrix(data[, 1:2], ncol=2)
# z <- as.matrix(data[, 4:5], ncol=2)
# y <- as.matrix(data[,3], ncol=1)
# solve(t(x)%*%x) %*% t(x) %*% y
# solve(t(z)%*%z) %*% t(z) %*% y
#
# apply(data,2,var)
# apply(estimate,1,mean)
#
#
# # ?optim
# # # some other optimization functions
# # optim(par = c(1, 1, rep(2, K), rep(4,R)), method = c("L-BFGS-B"), hessian = T, fn = SFM.within, Time = Time, N = N, x = x, y = y, z = z, lower= lower, upper= upper)
# # SFM.within(x = x, y = y, z = z, N = N,  Time = Time, par = c(1, 2, 1, 2, 3, 4), mu=0, optim = T)
# # cl <- makeCluster(spec=detectCores(), outfile="")
# # setDefaultCluster(cl=cl)
# # optimParallel(par = c(1, 2, rep(2, K), rep(4,R)), method = c("L-BFGS-B"), hessian = T, fn = SFM.within, Time = Time, N = N, x = x, y = y, z = z, lower= lower, upper= upper)
# # nlminb(objective = SFM.within, hessian=T,lower=c(0.0001,0.0001,-Inf,-Inf,-Inf,-Inf), start = c(1, 1, rep(1, K), rep(1,R)),Time = Time, N = N, x = x, y = y, z = z)
#
#
#
# ###############################################################################
# ##############          Test sfmfep                               #############
# ###############################################################################
#
# N <- 2
# Time <- 30
#
# beta <- c(0.5,3); K <- length(beta);
# delta <- c(0.5,2); R <- length(delta)
# sigma_u <- 2; sigma_v <- 1
#
# set.seed(665) # der seed läuft aktuell
# data.test <- SFM.generate(N = N, Time = Time, beta = beta, delta = delta, mu_u = 0, sigma_u = sigma_u, sigma_v = sigma_v)
# data.test <- cbind(rep(c("farmer","company"),N*Time/2), data.test)
# colnames(data.test) <- c("gr","x1", "x2","y","z1", "z2")
#
# form.test <- formula(y  ~ x1 + x2 + (z1 + z2))
# formula <- form.test
# data <- data.test
# data.test[1,2] <- Inf
# Time = NULL; N = NULL
# # do it with N & Time
# (sfmfep(formula = form.test, N = N, Time = Time, data = data.test, mu=0, myPar = NULL))
# # do it with N & Time
# (sfmfep(formula = form.test, group = "gr", data = data.test, mu=0, myPar = NULL))
#
# profvis(sfmfep(formula = form.test, R = 2, N = N, Time = Time, data = data.test, mu=0, myPar = NULL, sigma = c(0.1,0.05)))
#
#
# ###############################################################################
# ##############          Test SFM alpha index                      #############
# ###############################################################################
#
# N.input <- 20
# Time.input <- 10
#
# beta <- c(1,2); K <- length(beta);
# delta=c(3,4); R <- length(delta)
# sigma_u = 1; sigma_v = 2
# lower <- c(0,0, rep(-Inf,K), rep(-Inf,R)); upper <-  c(1000,1000,1000,1000,1000,1000)
#
# # Generate dataset for optim based on inputs above
# set.seed(666)
# data <- SFM.generate(N = N.input, Time = Time.input, beta = beta, delta = delta, mu_u = 0, sigma_u = sigma_u, sigma_v = sigma_v)
# x.dat <- as.matrix(data[,1:length(beta)]); y.dat <- as.matrix(data[,(length(beta)+1)])
# z.dat <- as.matrix(data[,(length(beta)+2):(length(beta)+2+length(delta)-1)])
#
# ret.list <- SFM.within(optim=F, N = N.input, Time = Time.input, xv = x.dat ,y = y.dat ,z = z.dat, par = c(1,2,1,2,3,4))
#
# # Hier werden
# test <- SFM.inindex(h=ret.list$h, sigma2star=ret.list$sigma_2star, mu2star=ret.list$mu_2star, N=N.input, Time=Time.input) # h is not within transformed
# m.test <- matrix(test,ncol=2)
#
#
# SFM.alpha(y = y.dat, x = x.dat, beta = c(1,2), mu = 0, sigma_u = 1, sigma_v <- 2, h = ret.list$h, epsilon = ret.list$epsilon, N=N.input, Time=Time.input)
#
#
#
# ############### RICE DATA SET
# data(riceProdPhil)
# rice <- (riceProdPhil)
# head(rice)
# table(rice$FMERCODE)
# SFM.caller(PRICE ~ LABORP+ NPKP + BANRAT, R = 1, Time = 8, N = 43, data =rice )
#
#
# # Some Information regarding optimizer
# # http://www.stat.colostate.edu/~jah/Computing_Hints/optimization.pdf
# # There are a number of general-purpose optimization routines in base R that I'm aware of: optim, nlminb, nlm and constrOptim (which handles linear inequality constraints, and calls optim under the hood). Here are some things that you might want to consider in choosing which one to use.
# #
# # optim can use a number of different algorithms including conjugate gradient, Newton, quasi-Newton, Nelder-Mead and simulated annealing. The last two don't need gradient information and so can be useful if gradients aren't available or not feasible to calculate (but are likely to be slower and require more parameter fine-tuning, respectively). It also has an option to return the computed Hessian at the solution, which you would need if you want standard errors along with the solution itself.
# #
# # nlminb uses a quasi-Newton algorithm that fills the same niche as the "L-BFGS-B" method in optim. In my experience it seems a bit more robust than optim in that it's more likely to return a solution in marginal cases where optim will fail to converge, although that's likely problem-dependent. It has the nice feature, if you provide an explicit gradient function, of doing a numerical check of its values at the solution. If these values don't match those obtained from numerical differencing, nlminb will give a warning; this helps to ensure you haven't made a mistake in specifying the gradient (easy to do with complicated likelihoods).
# #
# # nlm only uses a Newton algorithm. This can be faster than other algorithms in the sense of needing fewer iterations to reach convergence, but has its own drawbacks. It's more sensitive to the shape of the likelihood, so if it's strongly non-quadratic, it may be slower or you may get convergence to a false solution. The Newton algorithm also uses the Hessian, and computing that can be slow enough in practice that it more than cancels out any theoretical speedup.
# # When to use and not to use any particular method of maximization depends to a great extent on the type of data you have. nlm will work just fine if the likelihood surface isn't particularly "rough" and is everywhere differentiable. nlminb provides a way to constrain parameter values to particular bounding boxes. optim, which is probably the most-used optimizer, provides a few different optimization routines; for example, BFGS, L-BFGS-B, and simulated annealing (via the SANN option), the latter of which might be handy if you have a difficult optimizing problem. There are also a number of optimizers available on CRAN. rgenoud, for instance, provides a genetic algorithm for optimization. DEoptim uses a different genetic optimization routine. Genetic algorithms can be slow to converge, but are usually guaranteed to converge (in time) even when there are discontinuities in the likelihood. I don't know about DEoptim, but rgenoud is set up to use snow for parallel processing, which helps somewhat.
# #
# # So, a probably somewhat unsatisfactory answer is that you should use nlm or any other optimizer if it works for the data you have. If you have a well-behaved likelihood, any of the routines provided by  optim or nlm will give you the same result. Some may be faster than others, which may or may not matter, depending on the size of the dataset, etc. As for the number of parameters these routines can handle, I don't know, though it's probably quite a few. Of course, the more parameters you have, the more likely you are to run into problems with convergence.

#
# x <- sample(1000)

