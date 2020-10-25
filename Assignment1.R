# Question 3

# a)

set.seed(1)
n <- 500
# create matrix Z with columns representing Z1, Z2, Z3
Z <- matrix(0, nrow = n, ncol =3)
for (l in 1:3) {
  Z[,l] <- rnorm(n, mean=0, sd=1)
}
# create vectors for Y1 and Y2 and store them in a data frame for Y
Y1 <- 1 + Z[,1] 
Y2 <- 5 + 2*Z[,1] + Z[,2]
Y <- data.frame(Y1,Y2)

# calculate the index's of the missing values under MAR mechanism.
a <- 2; b <- 0
ind_mar <- which((a*(Y[,1]-1)+b*(Y[,2]-5)+Z[,3]) < 0)

# create vectors which store the observed values of Y2 and the missing values of Y2.
Y2_MAR_obs <- Y[,2][-ind_mar]
Y2_MAR_mis <- Y[,2][ind_mar]

# plot the marginal densities of Y2 (complete) and Y2 (observed, after imposing missingness)
plot(density(Y[,2]), lwd = 2, col = "blue", main = "MAR", 
     xlab = expression(Y[2]), ylim = c(0,0.3))
lines(density(Y2_MAR_obs), lwd = 2, col = "red")
legend(8,0.3, legend = c("Complete data", "Observed data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")

# b)

set.seed(1)
# impose missingness on Y2 under the MAR mechanism.
Y_mis <- Y
Y_mis[,2][ind_mar] <- NA
fitY2 <- lm(Y_mis[,2] ~ Y[,1], data = Y_mis)
summary(fitY2)
predsri <- predict(fitY2, newdata = Y_mis) + rnorm(n, 0, sigma(fitY2))
# impute values of Y2
Y2_sri <- ifelse(is.na(Y_mis[,2]) == TRUE, predsri, Y_mis[,2])

# plot the marginal densities of Y2 (complete) and Y2 (completed 
# after imputation)
plot(density(Y[,2]), lwd = 2, col = "blue", 
     main = "MAR", xlab = expression(Y[2]), ylim = c(0,0.3))
lines(density(Y2_sri, na.rm = TRUE), lwd = 2, col = "red")
legend(8,0.3, legend = c("Complete data", "Observed data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")

# c)

# calculate the index's of the missing values under MNAR mechanism.
a <- 0; b <- 2
ind_mnar <- which((a*(Y[,1]-1)+b*(Y[,2]-5)+Z[,3]) < 0)

# create vectors which store the observed values of Y2 and the missing values of Y2.
Y2_MNAR_obs <- Y[,2][-ind_mnar]
Y2_MNAR_mis <- Y[,2][ind_mnar]

# plot the marginal densities of Y2 (complete) and Y2 (observed, after imposing missingness)
plot(density(Y[,2]), lwd = 2, col = "blue", main = "MNAR", xlab = expression(Y[2]), ylim = c(0,0.3))
lines(density(Y2_MNAR_obs), lwd = 2, col = "red")
legend(8,0.3, legend = c("Complete data", "Observed data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")

# d)

set.seed(1)
# impose missingness on Y2 under the MNAR mechanism.
Y_mis <- Y
Y_mis[,2][ind_mnar] <- NA
# predict values for Y2 
fitY2 <- lm(Y_mis[,2] ~ Y[,1], data = Y_mis)
summary(fitY2)
predsri <- predict(fitY2, newdata = Y_mis) + rnorm(n, 0, sigma(fitY2))
# impute values of Y2
Y2_sri <- ifelse(is.na(Y_mis[,2]) == TRUE, predsri, Y_mis[,2])

# plot the marginal densities of Y2 (complete) and Y2 (completed after imputation)
plot(density(Y[,2]), lwd = 2, col = "blue", main = "MNAR", xlab = expression(Y[2]), ylim = c(0,0.3))
lines(density(Y2_sri), lwd = 2, col = "red")
legend(8,0.3, legend = c("Complete data", "Observed data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")


# Question 4

# a) complete case analysis

load("databp.Rdata")

ind <- which(is.na(databp$recovtime) == FALSE)
mccoverall <- mean(databp$recovtime, na.rm = TRUE)
seccoverall <- sd(databp$recovtime, na.rm = TRUE)/sqrt(length(ind))
mccoverall; seccoverall

rlogdose <- cor(databp$recovtime, databp$logdose, method = 'pearson', use = "complete.obs")
rbloodp <- cor(databp$recovtime, databp$bloodp, method = 'pearson', use = "complete.obs")
rlogdose; rbloodp

# b) mean imputation

recovtime_mi <- ifelse(is.na(databp$recovtime) == TRUE, 
                       mean(databp$recovtime, na.rm = TRUE), databp$recovtime)
n <- length(recovtime_mi)
mmi <- mean(recovtime_mi)
semi <- sd(recovtime_mi)/sqrt(n)
mmi; semi

rmilogdose <- cor(recovtime_mi, databp$logdose, method = 'pearson', use = "everything")
rmibloodp <- cor(recovtime_mi, databp$bloodp, method = 'pearson', use = "everything")
rmilogdose; rmibloodp

# c) mean regression imputation

fitrecovtime <- lm(recovtime ~ logdose + bloodp, data = databp)
summary(fitrecovtime)
predri <- predict(fitrecovtime, newdata = databp)
recovtime_ri <- ifelse(is.na(databp$recovtime) == TRUE, predri, databp$recovtime)
mri <- mean(recovtime_ri)
seri <- sd(recovtime_ri)/sqrt(n)
mri; seri

rmrilogdose <- cor(recovtime_ri, databp$logdose, method = 'pearson', use = "everything")
rmribloodp <- cor(recovtime_ri, databp$bloodp, method='pearson', use = "everything")
rmrilogdose; rmribloodp

# d) stochastic regression imputation

set.seed(1)
predsri <- predict(fitrecovtime, newdata = databp) + rnorm(n, 0, sigma(fitrecovtime))
recovtime_sri <- ifelse(is.na(databp$recovtime) == TRUE, predsri, databp$recovtime)
msri <- mean(recovtime_sri)
sesri <- sd(recovtime_sri)/sqrt(n)
msri; sesri

rsrilogdose <- cor(recovtime_sri, databp$logdose, method = 'pearson', use = "everything")
rsribloodp <- cor(recovtime_sri, databp$bloodp, method='pearson', use = "everything")
rsrilogdose; rsribloodp

# e) predictive mean matching

set.seed(1)
fitrecovtime <- lm(recovtime ~ logdose + bloodp, data = databp)
predhdi <- predict(fitrecovtime, newdata = databp)

donor <- c()

sqerr <- ((predhdi - predhdi[4])^2)
sqerr[4] <- Inf
ind4 <- which.min(sqerr)
donor <- c(donor, databp$recovtime[ind4])

sqerr <- ((predhdi - predhdi[10])^2)
sqerr[10] <- Inf
ind10 <- which.min(sqerr)
donor <- c(donor, databp$recovtime[ind10])

sqerr <- ((predhdi - predhdi[22])^2)
sqerr[22] <- Inf
ind22 <- which.min(sqerr)
donor <- c(donor, databp$recovtime[ind22])

recovtime_pmm <- databp$recovtime
recovtime_pmm[c(4,10,22)] <- donor

mpmm <- mean(recovtime_pmm)
sepmm <- sd(recovtime_pmm)/sqrt(n)
mpmm; sepmm

rpmmlogdose <- cor(recovtime_pmm, databp$logdose, method = 'pearson', use = "everything")
rpmmbloodp <- cor(recovtime_pmm, databp$bloodp, method='pearson', use = "everything")
rpmmlogdose; rpmmbloodp



