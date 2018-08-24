# Unbalanced Panels:
# TODO(Clemens) StartSigma -> affected in sfmfep ------ Done
# TODO(Clemens) Within model ------ Done
# TODO(Clemens) add error handler, if N & T vector doesnt match ------ Done
# TODO(Clemens) First Difference Model ------ Done
# TODO(Oli) SFM.inindex
# TODO(Oli) SFM.alpha
# TODO(Oli) Output

# First Difference Model
# TODO(Clemens) Funktion/Model ------ Done
# TODO(Clemens) sfmfep ------ Done
# TODO(Clemens) unit tester ------ Done
# TODO() output
# TODO(Oli) check BIC/AIC
# TODO(Clemens) check within model & first diff model ------ Done

# TODO(Oli) Alpha & Inefficency Index -> check results

# TODO(Oli) Output - s3 & tester

# Bootstrapping (individual)
# TODO(Clemens) Funktion ------ Done
# TODO(Clemens) sfmfep ------ Done
# TODO(Clemens) unit tester ------ Done
# TODO() output

# TODO() Dokumentation R
# TODO() proper commenting R
# TODO() Dokumentation Abgabe Dok
# TODO() R CMD checker / Travis CI
# TODO() Beispieldatensatz
# TODO(Oli & Clemens) Optimierung & Testläufe inkl. Simulationslauf a la paper




###############################################################################
##############          Test Optim / nlm                          #############
###############################################################################
# library("fepsfrontieR")
# # Specify "your" variables for the model (and lower & upper bounds for optim)
# N <- 100
# Time <- 10
# #
# Time <- rep(10, N)
# cumTime <- c(0, cumsum (Time))
# beta <- c(0.5); K <- length(beta);
# delta=c(0.5); R <- length(delta)
# sigma_u = 0.2; sigma_v = 0.1
# par <-c(sigma_u,sigma_v,beta,delta)
# mu <- 0
# lower <- c(0.05,0.05, rep(-Inf,K), rep(-Inf,R)); upper <-  c(1000,1000,1000,1000,1000,1000)
# # # # # # #
# # # # # # # # Generates data for the optimization & does it "Iteration" times
# Iteration <- 30
# estimateW <- matrix(c(rep(NA, (2+K+R)*Iteration)), nrow= 2 + K + R, ncol=Iteration)
# estimateF <- matrix(c(rep(NA, (2+K+R)*Iteration)), nrow= 2 + K + R, ncol=Iteration)
# profvis(
# for(i in 1:Iteration){
#   # # # # #   # Generate dataset for optim based on inputs above
#   set.seed(sample(99999,1))
#   data <- SFM.generate(N = N, Time = min(Time), beta = beta, delta = delta, mu = mu, sigma_u = sigma_u, sigma_v = sigma_v)
#   x    <- as.matrix(data[, 1:length(beta)]); y <- as.matrix(data[, (length(beta)+1)])
#   z    <- as.matrix(data[, (length(beta)+2):(length(beta)+2+length(delta)-1)])
#
#   # optim.SFMwithin <- nlminb (objective = SFM.within,
#   #                      start = c(sigma_u,sigma_v,beta,delta),
#   #                      lower = c(0.001,0.001,rep(-Inf, R+K)),
#   #                      mu = mu,
#   #                      Time = Time,
#   #                      N = N,
#   #                      R = length(delta), K = length(beta),
#   #                      cumTime = cumTime,
#   #                      xv = x, y = y, z = z,
#   #                      optim = T)
#
#   optim.SFMfirstDiff <- nlminb (objective = SFM.firstDiff,
#                        start = c(sigma_u,sigma_v,beta,delta),
#                        lower = c(0.001,0.001,rep(-Inf, R+K)),
#                        mu = mu,
#                        Time = Time,
#                        N = N,
#                        cumTime = cumTime,
#                        K = K, R = R,
#                        xv = x, y = y, z = z,
#                        optim = T)
#
#    # estimateW[, i] <- optim.SFMwithin$par
#   estimateF[, i] <- optim.SFMfirstDiff$par
# })
#
# apply(estimateF, 1, median); apply(estimateF, 1, mean)
#
# t.formula <- y ~ x + (z)
# ttt <- sfmfep(formula = t.formula, data = data, N = 100, Time = 5, mu = 0.5,
#        sigmaCI = 0.05, estimate = T, method = "firstdiff", bootstrap = T, B = 100,
#                    myPar = c(sigma_u = sigma_u, sigma_v = sigma_v, beta = beta, delta = delta))
# hist(ttt$estimatesMat[,1])
#
# # #
# #
# #   hist(estimateW[1,])
# apply(estimateW,1, median) - apply(estimateF,1, median)
# apply(estimateW,1, mean) - c(sigma_u = sigma_u, sigma_v = sigma_v, beta = beta, delta = delta)
# apply(estimateW,1, mean)
# apply(estimateW,1, var)
# apply(estimateF,1, mean)
# apply(estimateF,1, var)

# for(i in 1:Iteration)
# # # # typeof(dplyr::tbl_df(data))
# # # # #   typeof(data)
# # # # # x <- as.matrix(data[, 1:2], ncol=2)
# # # # # z <- as.matrix(data[, 4:5], ncol=2)
# # # # # y <- as.matrix(data[,3], ncol=1)
# # # # # solve(t(x)%*%x) %*% t(x) %*% y
# # # # # solve(t(z)%*%z) %*% t(z) %*% y
# # # # #
# # # # # apply(data,2,var)
# # # # # apply(estimate,1,mean)
# # # # #
# #

# # # # #
# # # # # ############### RICE DATA SET
# # # # # data(riceProdPhil)
# # # # # rice <- (riceProdPhil)
# # # # # head(rice)
# # # # # table(rice$FMERCODE)
# # # # # SFM.caller(PRICE ~ LABORP+ NPKP + BANRAT, R = 1, Time = 8, N = 43, data =rice )
# # # # #
# # # # #
# # # # # # Some Information regarding optimizer
# # # # # # http://www.stat.colostate.edu/~jah/Computing_Hints/optimization.pdf
# # # # # # There are a number of general-purpose optimization routines in base R that I'm aware of: optim, nlminb, nlm and constrOptim (which handles linear inequality constraints, and calls optim under the hood). Here are some things that you might want to consider in choosing which one to use.
# # # # # #
# # # # # # optim can use a number of different algorithms including conjugate gradient, Newton, quasi-Newton, Nelder-Mead and simulated annealing. The last two don't need gradient information and so can be useful if gradients aren't available or not feasible to calculate (but are likely to be slower and require more parameter fine-tuning, respectively). It also has an option to return the computed Hessian at the solution, which you would need if you want standard errors along with the solution itself.
# # # # # #
# # # # # # nlminb uses a quasi-Newton algorithm that fills the same niche as the "L-BFGS-B" method in optim. In my experience it seems a bit more robust than optim in that it's more likely to return a solution in marginal cases where optim will fail to converge, although that's likely problem-dependent. It has the nice feature, if you provide an explicit gradient function, of doing a numerical check of its values at the solution. If these values don't match those obtained from numerical differencing, nlminb will give a warning; this helps to ensure you haven't made a mistake in specifying the gradient (easy to do with complicated likelihoods).
# # # # # #
# # # # # # nlm only uses a Newton algorithm. This can be faster than other algorithms in the sense of needing fewer iterations to reach convergence, but has its own drawbacks. It's more sensitive to the shape of the likelihood, so if it's strongly non-quadratic, it may be slower or you may get convergence to a false solution. The Newton algorithm also uses the Hessian, and computing that can be slow enough in practice that it more than cancels out any theoretical speedup.
# # # # # # When to use and not to use any particular method of maximization depends to a great extent on the type of data you have. nlm will work just fine if the likelihood surface isn't particularly "rough" and is everywhere differentiable. nlminb provides a way to constrain parameter values to particular bounding boxes. optim, which is probably the most-used optimizer, provides a few different optimization routines; for example, BFGS, L-BFGS-B, and simulated annealing (via the SANN option), the latter of which might be handy if you have a difficult optimizing problem. There are also a number of optimizers available on CRAN. rgenoud, for instance, provides a genetic algorithm for optimization. DEoptim uses a different genetic optimization routine. Genetic algorithms can be slow to converge, but are usually guaranteed to converge (in time) even when there are discontinuities in the likelihood. I don't know about DEoptim, but rgenoud is set up to use snow for parallel processing, which helps somewhat.
# # # # # #
# # # # # # So, a probably somewhat unsatisfactory answer is that you should use nlm or any other optimizer if it works for the data you have. If you have a well-behaved likelihood, any of the routines provided by  optim or nlm will give you the same result. Some may be faster than others, which may or may not matter, depending on the size of the dataset, etc. As for the number of parameters these routines can handle, I don't know, though it's probably quite a few. Of course, the more parameters you have, the more likely you are to run into problems with convergence.
