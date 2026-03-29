# ==============================================================================
# Global variables declaration for R CMD check
# These include JAGS-specific syntax and internal modeling variables
# ==============================================================================
utils::globalVariables(c(
  "N", "theta", "bx0", "r", "a", "b", "c", "d", "vx0", "mx0",
  "K", "Tx", "mx", "sx0", "A", "B", "R", "Weights", "mu", "bx", "rho", "phi",
  "by", "sy", "sx", "inprod", "T", "kappa", "mod_noOutliers"
))

# ==============================================================================
# The following functions are adapted/copied from the TwoSampleMR package
# Original Authors: Gibran Hemani, Jie Zheng, et al.
# Source: https://github.com/MRCIEU/TwoSampleMR
# ==============================================================================

default_parameters_mr_grip <- function() {
  list(
    test_dist = "z", nboot = 1000, Cov = 0, penk = 20, phi = 1,
    alpha = 0.05, Qthresh = 0.05, over.dispersion = TRUE,
    loss.function = "huber", shrinkage = FALSE
  )
}

mr_grip <- function(b_exp, b_out, se_exp, se_out, parameters) {
  if (length(b_exp) != length(b_out)) stop("The lengths of b_exp and b_out are not equal.", call. = FALSE)
  if (length(se_exp) != length(se_out)) stop("The lengths of se_exp and se_out are not equal.", call. = FALSE)

  nulllist <- list(
    b = NA, se = NA, pval = NA, b_i = NA, se_i = NA, pval_i = NA,
    b.adj = NA, se.adj = NA, pval.adj = NA, nsnp = NA, mod = NA, smod = NA, dat = NA
  )
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3) return(nulllist)

  dat <- data.frame(b_out = b_out, b_exp = b_exp, se_exp = se_exp, se_out = se_out)
  grip_out <- b_out * b_exp
  grip_exp <- b_exp^2
  grip_weights <- 1 / (b_exp^2 * se_out^2)

  mod <- stats::lm(grip_out ~ grip_exp, weights = grip_weights)
  smod <- summary(mod)
  b <- stats::coefficients(smod)[2, 1]
  se <- stats::coefficients(smod)[2, 2]
  b_i <- stats::coefficients(smod)[1, 1]
  se_i <- stats::coefficients(smod)[1, 2]

  grip_weights_adj <- 1 / (se_out^2)
  numer <- sum(grip_weights_adj) * sum(grip_weights_adj * b_out * b_exp * (b_exp^2 - 3 * se_exp^2)) -
    sum(grip_weights_adj * b_out * b_exp) * sum(grip_weights_adj * (b_exp^2 - se_exp^2))
  denom <- sum(grip_weights_adj) * sum(grip_weights_adj * (b_exp^4 - 6 * b_exp^2 * se_exp^2 + 3 * se_exp^4)) -
    (sum(grip_weights_adj * (b_exp^2 - se_exp^2)))^2
  b.adj <- numer / denom

  var_out <- mean((b_out - (b * grip_exp + b_i) / b_exp)^2) * b_exp^2
  numer_se <- sum(grip_weights_adj)^2 * sum(grip_weights_adj^2 * var_out * (b_exp^2 - 3 * se_exp^2)^2) +
    sum(grip_weights_adj^2 * var_out) * sum(grip_weights_adj * (b_exp^2 - se_exp^2))^2 -
    2 * sum(grip_weights_adj) * sum(grip_weights_adj * (b_exp^2 - se_exp^2)) *
    sum(grip_weights_adj^2 * (b_exp^2 - se_exp^2) * var_out)
  se.adj <- sqrt(numer_se) / denom

  pval <- 2 * stats::pt(abs(b / se), length(b_exp) - 2, lower.tail = FALSE)
  pval_i <- 2 * stats::pt(abs(b_i / se_i), length(b_exp) - 2, lower.tail = FALSE)
  pval.adj <- 2 * stats::pt(abs(b.adj / se.adj), length(b_exp) - 2, lower.tail = FALSE)

  list(b = b, se = se, pval = pval, b_i = b_i, se_i = se_i, pval_i = pval_i,
       b.adj = b.adj, se.adj = se.adj, pval.adj = pval.adj, nsnp = length(b_exp),
       mod = mod, smod = smod, dat = dat)
}

# ==============================================================================
# The following functions are adapted/copied from the MR-Horse implementation
# Original Authors: A. J. Grant, Stephen Burgess, et al.
# ==============================================================================

mr_horse_model <- function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1/(sy[i] * sy[i]))
    mu[i] = theta * bx0[i] + alpha[i]
    bx[i] ~ dnorm(bx0[i], 1 / (sx[i] * sx[i]))
    bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0))
    r[i] ~ dbeta(10, 10);T(, 1)
    rho[i] = 2*r[i] - 1
    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
    phi[i] = a[i] / sqrt(b[i])
    a[i] ~ dnorm(0, 1);T(0, )
    b[i] ~ dgamma(0.5, 0.5)
  }
  c ~ dnorm(0, 1);T(0, )
  d ~ dgamma(0.5, 0.5)
  tau = c / sqrt(d)
  vx0 ~ dnorm(0, 1);T(0, )
  mx0 ~ dnorm(0, 1)
  theta ~ dunif(-10, 10)
}

mvmr_horse_model <- function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1 / (sy[i] * sy[i]))
    mu[i] = inprod(bx0[i, 1:K], theta) + alpha[i]
    bx[i,1:K] ~ dmnorm(bx0[i,1:K], Tx[1:K, ((i-1)*K+1):(i*K)])
    kappa[i] = (rho[i]^2 / (1 + K*rho[i]^2))
    bx0[i,1:K] ~ dmnorm(mx + sx0 * rho[i] * alpha[i] / (phi[i] * tau), A - kappa[i] * B)
    r[i] ~ dbeta(10, 10);T(, 1)
    rho[i] = 2*r[i] - 1
    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
    phi[i] = a[i] / sqrt(b[i])
    a[i] ~ dnorm(0, 1);T(0, )
    b[i] ~ dgamma(0.5, 0.5)
  }
  c ~ dnorm(0, 1);T(0, )
  d ~ dgamma(0.5, 0.5)
  tau = c / sqrt(d)
  mx ~ dmnorm(rep(0, K), R[,])
  for (k in 1:K){
    vx0[k] ~ dnorm(0, 1);T(0, )
    sx0[k] = sqrt(vx0[k])
    theta[k] ~ dunif(-10, 10)
    for (j in 1:K){
      A[j, k] = ifelse(j==k, 1/vx0[j], 0)
      B[j, k] = 1 / (sx0[j] * sx0[k])
    }
  }
}

mr_horse <- function(D, no_ini = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000){
  if(!("theta" %in% variable.names)) variable.names <- c("theta", variable.names)
  jags_fit <- R2jags::jags(
    data = list(by = D$betaY, bx = D$betaX, sy = D$betaYse, sx = D$betaXse, N = length(D$betaY)),
    parameters.to.save = variable.names, n.chains = no_ini, n.iter = n.burnin + n.iter, n.burnin = n.burnin, model.file = mr_horse_model
  )
  mr.coda <- coda::as.mcmc(jags_fit)
  summ <- summary(mr.coda[, "theta"])
  mr_estimate <- data.frame(
    "Estimate" = round(unname(summ$statistics[1]), 3), "SD" = round(unname(summ$statistics[2]), 3),
    "2.5% quantile" = round(unname(summ$quantiles[1]), 3), "97.5% quantile" = round(unname(summ$quantiles[5]), 3),
    "Rhat" = round(unname(coda::gelman.diag(mr.coda)$psrf[1]), 3)
  )
  names(mr_estimate) <- c("Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda)
}

mvmr_horse <- function(D, no_ini = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000){
  if(!("theta" %in% variable.names)) variable.names <- c("theta", variable.names)
  p <- dim(D)[1]
  K_val <- sum(sapply(1:dim(D)[2], function(j) substr(names(D)[j], 1, 5) == "betaX"))/2
  Bx <- D[, sprintf("betaX%i", 1:K_val)]
  Sx <- D[, sprintf("betaX%ise", 1:K_val)]
  Tx_mat <- matrix(nrow = K_val, ncol = p*K_val)
  for (j in 1:p) Tx_mat[, ((j-1)*K_val+1):(j*K_val)] <- diag(1 / Sx[j, ]^2)
  jags_fit <- R2jags::jags(
    data = list(by = D$betaY, bx = Bx, sy = D$betaYse, Tx = Tx_mat, N = p, K = K_val, R = diag(K_val)),
    parameters.to.save = variable.names, n.chains = no_ini, n.iter = n.burnin + n.iter, n.burnin = n.burnin, model.file = mvmr_horse_model
  )
  mr.coda <- coda::as.mcmc(jags_fit)
  s <- summary(mr.coda)
  mr_estimate <- data.frame(
    "Parameter" = sprintf("theta[%i]", 1:K_val),
    "Estimate" = round(unname(s$statistics[sprintf("theta[%i]", 1:K_val), 1]), 3),
    "SD" = round(unname(s$statistics[sprintf("theta[%i]", 1:K_val), 2]), 3),
    "2.5% quantile" = round(unname(s$quantiles[sprintf("theta[%i]", 1:K_val), 1]), 3),
    "97.5% quantile" = round(unname(s$quantiles[sprintf("theta[%i]", 1:K_val), 5]), 3),
    "Rhat" = round(unname(coda::gelman.diag(mr.coda)$psrf[sprintf("theta[%i]", 1:K_val), 1]), 3)
  )
  names(mr_estimate) <- c("Parameter", "Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda)
}

# ==============================================================================
# The following function(s) are adapted/copied from the MRPRESSO package
# Original Author: Marie Verbanck (GPL-3)
# ==============================================================================

mr_presso <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data,
                      OUTLIERtest = FALSE, DISTORTIONtest = FALSE,
                      SignifThreshold = 0.05, NbDistribution = 1000, seed = NULL){

  if(!is.null(seed)) set.seed(seed)
  if(!is.data.frame(data)) stop("data must be an object of class data.frame", call. = FALSE)

  "%^%" <- function(x, n) with(eigen(x), vectors %*% (values^n * t(vectors)))

  getRandomData <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data){
    mod_IVW <- lapply(1:nrow(data), function(i)
      stats::lm(stats::as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse=" + "))),
                weights = Weights, data = data[-i, ]))
    dataRandom <- cbind(eval(parse(text = paste0(
      "cbind(", paste0("rnorm(nrow(data), data[, '", BetaExposure, "'], data[, '", SdExposure, "'])", collapse = ", "),
      ", sapply(1:nrow(data), function(i) rnorm(1, predict(mod_IVW[[i]], newdata = data[i, ]), data[i ,'", SdOutcome,"'])))"
    ))), data$Weights)
    colnames(dataRandom) <- c(BetaExposure, BetaOutcome, "Weights")
    dataRandom
  }

  data <- data[, c(BetaOutcome, BetaExposure, SdOutcome, SdExposure)]
  data <- data[rowSums(is.na(data)) == 0, ]
  data[, c(BetaOutcome, BetaExposure)] <- data[, c(BetaOutcome, BetaExposure)] * sign(data[, BetaExposure[1]])
  data$Weights <- 1 / data[, SdOutcome]^2

  if(nrow(data) <= length(BetaExposure) + 2) stop("Not enough instrumental variables", call. = FALSE)
  if(nrow(data) >= NbDistribution) stop("Not enough elements to compute empirical P-values, increase NbDistribution", call. = FALSE)

  getRSS_LOO <- function(BetaOutcome, BetaExposure, data, returnIV){
    dataW <- data[, c(BetaOutcome, BetaExposure)] * sqrt(data[, "Weights"])
    X <- as.matrix(dataW[ , BetaExposure]); Y <- as.matrix(dataW[ , BetaOutcome])
    CausalEstimate_LOO <- sapply(1:nrow(dataW), function(i) {
      (t(X[-i, ]) %*% X[-i, ])%^%(-1) %*% t(X[-i, ]) %*% Y[-i, ]
    })
    RSS <- if(length(BetaExposure) == 1) sum((Y - CausalEstimate_LOO * X)^2, na.rm = TRUE) else sum((Y - rowSums(t(CausalEstimate_LOO) * X))^2, na.rm = TRUE)
    if(returnIV) list(RSS, CausalEstimate_LOO) else RSS
  }

  RSSobs <- getRSS_LOO(BetaOutcome, BetaExposure, data, returnIV = OUTLIERtest)
  randomData <- replicate(NbDistribution, getRandomData(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data), simplify = FALSE)
  RSSexp <- sapply(randomData, getRSS_LOO, BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, returnIV = OUTLIERtest)

  if(OUTLIERtest) GlobalTest <- list(RSSobs = RSSobs[[1]], Pvalue = sum(RSSexp[1, ] > RSSobs[[1]])/NbDistribution)
  else GlobalTest <- list(RSSobs = RSSobs, Pvalue = sum(RSSexp > RSSobs)/NbDistribution)

  if(GlobalTest$Pvalue < SignifThreshold & OUTLIERtest){
    OutlierTest <- do.call("rbind", lapply(1:nrow(data), function(SNV){
      randomSNP <- do.call("rbind", lapply(randomData, function(mat) mat[SNV, ]))
      if(length(BetaExposure) == 1){
        Dif <- data[SNV, BetaOutcome] - data[SNV, BetaExposure] * RSSobs[[2]][SNV]
        Exp <- randomSNP[, BetaOutcome] - randomSNP[, BetaExposure] * RSSobs[[2]][SNV]
      } else {
        Dif <- data[SNV, BetaOutcome] - sum(data[SNV, BetaExposure] * RSSobs[[2]][, SNV])
        Exp <- randomSNP[, BetaOutcome] - rowSums(randomSNP[, BetaExposure] * RSSobs[[2]][, SNV])
      }
      pval <- sum(Exp^2 > Dif^2)/length(randomData)
      cbind.data.frame(RSSobs = Dif^2, Pvalue = pval)
    }))
    OutlierTest$Pvalue <- apply(cbind(OutlierTest$Pvalue*nrow(data), 1), 1, min)
  } else { OUTLIERtest <- FALSE }

  mod_all <- stats::lm(stats::as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse = "+"))), weights = Weights, data = data)

  if(DISTORTIONtest & OUTLIERtest){
    refOutlier <- which(OutlierTest$Pvalue <= SignifThreshold)
    if(length(refOutlier) > 0 && length(refOutlier) < nrow(data)){
      mod_noOutliers <- stats::lm(stats::as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse=" + "))), weights = Weights, data = data[-refOutlier, ])
      BiasObs <- (mod_all$coefficients[BetaExposure] - mod_noOutliers$coefficients[BetaExposure]) / abs(mod_noOutliers$coefficients[BetaExposure])
      BiasTest <- list(`Outliers Indices` = refOutlier, `Distortion Coefficient` = 100*BiasObs, Pvalue = NA)
    } else { BiasTest <- list(`Outliers Indices` = "No valid subsets for distortion test", `Distortion Coefficient` = NA, Pvalue = NA) }
  }

  OriginalMR <- cbind.data.frame(BetaExposure, "Raw", summary(mod_all)$coefficients)
  colnames(OriginalMR) <- c("Exposure", "MR Analysis", "Causal Estimate", "Sd", "T-stat", "P-value")
  MR_final <- OriginalMR
  if(exists("mod_noOutliers") && !is.null(mod_noOutliers)){
    OutlierCorrectedMR <- cbind.data.frame(BetaExposure, "Outlier-corrected", summary(mod_noOutliers)$coefficients, row.names = NULL)
    colnames(OutlierCorrectedMR) <- colnames(OriginalMR)
    MR_final <- rbind.data.frame(OriginalMR, OutlierCorrectedMR)
  }

  list(`Main MR results` = MR_final, `MR-PRESSO results` = list(`Global Test` = GlobalTest))
}
