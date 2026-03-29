# ==============================================================================
# MR-RAPS (vendored) - internal implementation for this package
#
# The functions below are adapted/copied from the mr.raps package (MR-RAPS)
# Original Authors: Qingyuan Zhao, Jingshu Wang, et al.
# License: GPL-3
# Source: https://github.com/qingyuanzhao/mr.raps
#
# Notes:
# - Copied/adapted on: 2026-03-05
# ==============================================================================

# ---- internal helpers ----

mr_raps_check_numeric_vec <- function(x, name) {
  if (!is.numeric(x) || !is.vector(x)) stop(sprintf("%s must be a numeric vector.", name), call. = FALSE)
  if (any(!is.finite(x))) stop(sprintf("%s contains non-finite values.", name), call. = FALSE)
}

mr_raps_check_lengths <- function(...) {
  xs <- list(...)
  n <- length(xs[[1]])
  if (!all(vapply(xs, length, integer(1)) == n)) stop("All inputs must have the same length.", call. = FALSE)
  if (n < 2) stop("Need at least 2 variants (SNPs).", call. = FALSE)
}

# ---- robust loss functions ----

mr_raps_rho_huber <- function(r, k = 1.345, deriv = 0) {
  if (deriv == 0) {
    ifelse(abs(r) <= k, r^2 / 2, k * (abs(r) - k / 2))
  } else if (deriv == 1) {
    ifelse(abs(r) <= k, r, k * sign(r))
  } else if (deriv == 2) {
    ifelse(abs(r) <= k, 1, 0)
  } else {
    stop("deriv must be 0, 1, or 2.", call. = FALSE)
  }
}

mr_raps_rho_tukey <- function(r, k = 4.685, deriv = 0) {
  if (deriv == 0) {
    pmin(1 - (1 - (r / k)^2)^3, 1)
  } else if (deriv == 1) {
    r * (1 - (r / k)^2)^2 * (abs(r) <= k)
  } else if (deriv == 2) {
    t <- (r / k)^2
    ifelse(t < 1, (1 - t) * (1 - 5 * t), 0)
  } else {
    stop("deriv must be 0, 1, or 2.", call. = FALSE)
  }
}

# ---- main exported-by-name (but keep it INTERNAL: do NOT @export) ----
# IMPORTANT: do NOT add roxygen @export tag; then it won't be exported to users.
mr.raps <- function(b_exp, b_out, se_exp, se_out,
                    over.dispersion = FALSE,
                    loss.function = c("l2", "huber", "tukey"),
                    diagnosis = FALSE,
                    se.method = c("sandwich", "bootstrap"),
                    k = switch(loss.function[1], l2 = NULL, huber = 1.345, tukey = 4.685),
                    B = 1000,
                    suppress.warning = FALSE) {

  loss.function <- match.arg(loss.function, c("l2", "huber", "tukey"))
  se.method <- match.arg(se.method, c("sandwich", "bootstrap"))

  mr_raps_check_numeric_vec(b_exp, "b_exp")
  mr_raps_check_numeric_vec(b_out, "b_out")
  mr_raps_check_numeric_vec(se_exp, "se_exp")
  mr_raps_check_numeric_vec(se_out, "se_out")
  mr_raps_check_lengths(b_exp, b_out, se_exp, se_out)

  if (any(se_exp <= 0) || any(se_out <= 0)) stop("Standard errors must be > 0.", call. = FALSE)

  if (loss.function == "l2") {
    fit <- if (!over.dispersion) {
      mr_raps_simple(b_exp, b_out, se_exp, se_out, diagnosis = diagnosis)
    } else {
      mr_raps_overdispersed(
        b_exp, b_out, se_exp, se_out,
        diagnosis = diagnosis,
        suppress.warning = suppress.warning
      )
    }
  } else {
    fit <- if (!over.dispersion) {
      mr_raps_simple_robust(
        b_exp, b_out, se_exp, se_out,
        loss.function = loss.function,
        k = k,
        diagnosis = diagnosis
      )
    } else {
      mr_raps_overdispersed_robust(
        b_exp, b_out, se_exp, se_out,
        loss.function = loss.function,
        k = k,
        suppress.warning = suppress.warning,
        diagnosis = diagnosis
      )
    }
  }

  if (se.method == "bootstrap") {
    fit.bootstrap <- vector("list", B)
    for (b in seq_len(B)) {
      s <- sample.int(length(b_exp), replace = TRUE)
      fit.bootstrap[[b]] <- tryCatch(
        unlist(mr.raps(
          b_exp[s], b_out[s], se_exp[s], se_out[s],
          over.dispersion = over.dispersion,
          loss.function = loss.function,
          diagnosis = FALSE,
          se.method = "sandwich",
          k = k,
          suppress.warning = TRUE
        )),
        error = function(e) NA
      )
    }
    fit.bootstrap <- data.frame(do.call(rbind, fit.bootstrap))
    fit <- c(
      fit,
      list(
        beta.hat.bootstrap = stats::median(fit.bootstrap$beta.hat),
        beta.se.bootstrap = stats::mad(fit.bootstrap$beta.hat)
      )
    )
  }

  fit
}

# ---- l2, no overdispersion ----

mr_raps_simple <- function(b_exp, b_out, se_exp, se_out, diagnosis = FALSE) {
  profile.loglike <- function(beta) {
    -0.5 * sum((b_out - b_exp * beta)^2 / (se_out^2 + se_exp^2 * beta^2))
  }

  bound <- stats::quantile(abs(b_out / b_exp), 0.95) * 2
  beta.hat <- stats::optimize(profile.loglike, bound * c(-1, 1), maximum = TRUE,
                              tol = .Machine$double.eps^0.5)$maximum
  while (abs(beta.hat) > 0.95 * bound) {
    bound <- bound * 2
    beta.hat <- stats::optimize(profile.loglike, bound * c(-1, 1), maximum = TRUE,
                                tol = .Machine$double.eps^0.5)$maximum
  }

  score.var <- sum(((b_exp^2 - se_exp^2) * se_out^2 +
                      (b_out^2 - se_out^2) * se_exp^2 +
                      se_exp^2 * se_out^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)

  I <- sum(((b_exp^2 - se_exp^2) * se_out^2 +
              (b_out^2 - se_out^2) * se_exp^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)

  dif <- b_out - beta.hat * b_exp
  dif.var <- se_out^2 + beta.hat^2 * se_exp^2
  chi.sq.test <- sum((dif / sqrt(dif.var))^2)

  if (diagnosis) {
    std.resid <- (b_out - b_exp * beta.hat) / sqrt(se_out^2 + beta.hat^2 * se_exp^2)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = c(1, 2))
    stats::qqnorm(std.resid); graphics::abline(0, 1)

    if (requireNamespace("nortest", quietly = TRUE)) {
      ad <- nortest::ad.test(std.resid)
      message("Anderson-Darling test: statistic = ", round(ad$statistic, 4),
              ", p-value = ", round(ad$p.value, 4))
    }
    sw <- stats::shapiro.test(std.resid)
    message("Shapiro-Wilk test: statistic = ", round(sw$statistic, 4),
            ", p-value = ", round(sw$p.value, 4))
  }

  list(
    beta.hat = beta.hat,
    beta.se = sqrt(score.var / I^2),
    beta.p.value = min(1, 2 * (1 - stats::pnorm(abs(beta.hat) / sqrt(score.var / I^2)))),
    naive.se = sqrt(1 / I),
    chi.sq.test = chi.sq.test
  )
}

# ---- l2, with overdispersion ----

mr_raps_overdispersed <- function(b_exp, b_out, se_exp, se_out,
                                  initialization = c("simple", "mode"),
                                  suppress.warning = FALSE,
                                  diagnosis = FALSE,
                                  niter = 20,
                                  tol = .Machine$double.eps^0.5) {

  initialization <- match.arg(initialization, c("simple", "mode"))

  profile.loglike.fixbeta <- function(beta, tau2) {
    alpha.hat <- 0
    -0.5 * sum(se_exp^2 * (log(tau2 + se_out^2 + se_exp^2 * beta^2))) -
      0.5 * sum(se_exp^2 * (b_out - alpha.hat - b_exp * beta)^2 / (tau2 + se_out^2 + se_exp^2 * beta^2))
  }

  profile.loglike.fixtau <- function(beta, tau2) {
    alpha.hat <- 0
    -0.5 * sum((b_out - alpha.hat - b_exp * beta)^2 / (tau2 + se_out^2 + se_exp^2 * beta^2))
  }

  bound.beta <- stats::quantile(abs(b_out / b_exp), 0.95) * 10
  bound.tau2 <- stats::quantile(se_out^2, 0.95) * 10

  if (initialization == "mode") {
    stop("Initialization by mode estimator is currently not supported.", call. = FALSE)
  } else {
    fit <- mr_raps_simple(b_exp, b_out, se_exp, se_out)
    beta.hat <- fit$beta.hat
    tau2.hat <- 0
  }

  for (iter in 1:niter) {
    beta.hat.old <- beta.hat
    tau2.hat.old <- tau2.hat

    tau2.hat <- stats::optimize(
      function(tau2) profile.loglike.fixbeta(beta.hat, tau2),
      bound.tau2 * c(0, 1),
      maximum = TRUE,
      tol = .Machine$double.eps^0.5
    )$maximum
    if (tau2.hat < 0) tau2.hat <- 0
    if (tau2.hat > 0.95 * bound.tau2) warning("Estimated overdispersion seems abnormally large.", call. = FALSE)

    beta.hat <- stats::optimize(
      function(beta) profile.loglike.fixtau(beta, tau2.hat),
      bound.beta * c(-1, 1),
      maximum = TRUE,
      tol = .Machine$double.eps^0.5
    )$maximum

    if (abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10) +
        abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) <= tol) {
      break
    }
  }

  if ((tau2.hat <= min(se_out^2) / 5) && (!suppress.warning)) {
    warning("The estimated overdispersion parameter is very small. Consider using the simple model.", call. = FALSE)
  }

  score.var <- diag(c(
    sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) +
           (b_out^2 - tau2.hat - se_out^2) * se_exp^2 +
           se_exp^2 * (tau2.hat + se_out^2)) / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
    sum(2 * se_exp^4 / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2)
  ))

  I <- matrix(c(
    -sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) +
            (b_out^2 - tau2.hat - se_out^2) * se_exp^2) / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
    0,
    -sum(se_exp^2 * beta.hat / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
    -sum(se_exp^2 / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2)
  ), 2, 2)

  asymp.var <- solve(I) %*% score.var %*% t(solve(I))

  out <- list(
    beta.hat = beta.hat,
    tau2.hat = tau2.hat,
    beta.se = sqrt(asymp.var[1, 1]),
    tau2.se = sqrt(asymp.var[2, 2]),
    beta.p.value = min(1, 2 * (1 - stats::pnorm(abs(beta.hat) / sqrt(asymp.var[1, 1]))))
  )

  if (diagnosis) {
    std.resid <- (b_out - b_exp * beta.hat) / sqrt(tau2.hat + se_out^2 + beta.hat^2 * se_exp^2)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = c(1, 2))
    stats::qqnorm(std.resid); graphics::abline(0, 1)
    if (requireNamespace("nortest", quietly = TRUE)) {
      ad <- nortest::ad.test(std.resid)
      message("Anderson-Darling test: statistic = ", round(ad$statistic, 4),
              ", p-value = ", round(ad$p.value, 4))
    }
    sw <- stats::shapiro.test(std.resid)
    message("Shapiro-Wilk test: statistic = ", round(sw$statistic, 4),
            ", p-value = ", round(sw$p.value, 4))
    out$std.resid <- std.resid
  }

  out
}

# ---- robust, no overdispersion ----

mr_raps_simple_robust <- function(b_exp, b_out, se_exp, se_out,
                                  loss.function = c("huber", "tukey"),
                                  k = switch(loss.function[1], huber = 1.345, tukey = 4.685),
                                  diagnosis = FALSE) {

  loss.function <- match.arg(loss.function, c("huber", "tukey"))
  rho <- switch(
    loss.function,
    huber = function(r, ...) mr_raps_rho_huber(r, k, ...),
    tukey = function(r, ...) mr_raps_rho_tukey(r, k, ...)
  )

  delta <- stats::integrate(function(x) x * rho(x, deriv = 1) * stats::dnorm(x), -Inf, Inf)$value
  c1 <- stats::integrate(function(x) rho(x, deriv = 1)^2 * stats::dnorm(x), -Inf, Inf)$value

  robust.loglike <- function(beta) {
    -sum(rho((b_out - b_exp * beta) / sqrt(se_out^2 + se_exp^2 * beta^2)))
  }

  bound <- stats::quantile(abs(b_out / b_exp), 0.95) * 2
  beta.hat <- stats::optimize(robust.loglike, bound * c(-1, 1), maximum = TRUE,
                              tol = .Machine$double.eps^0.5)$maximum
  while (abs(beta.hat) > 0.95 * bound) {
    bound <- bound * 2
    beta.hat <- stats::optimize(robust.loglike, bound * c(-1, 1), maximum = TRUE,
                                tol = .Machine$double.eps^0.5)$maximum
  }

  score.var <- c1 * sum(((b_exp^2 - se_exp^2) * se_out^2 +
                           (b_out^2 - se_out^2) * se_exp^2 +
                           se_exp^2 * se_out^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)

  I <- delta * sum(((b_exp^2 - se_exp^2) * se_out^2 +
                      (b_out^2 - se_out^2) * se_exp^2) / (se_out^2 + beta.hat^2 * se_exp^2)^2)

  dif <- b_out - beta.hat * b_exp
  dif.var <- se_out^2 + beta.hat^2 * se_exp^2
  chi.sq.test <- sum((dif / sqrt(dif.var))^2)

  list(
    beta.hat = beta.hat,
    beta.se = sqrt(score.var / I^2),
    naive.se = sqrt(1 / I),
    chi.sq.test = chi.sq.test,
    beta.p.value = min(1, 2 * (1 - stats::pnorm(abs(beta.hat) / sqrt(score.var / I^2))))
  )
}

# ---- robust, with overdispersion ----

mr_raps_overdispersed_robust <- function(b_exp, b_out, se_exp, se_out,
                                         loss.function = c("huber", "tukey"),
                                         k = switch(loss.function[1], huber = 1.345, tukey = 4.685),
                                         initialization = c("l2", "mode"),
                                         suppress.warning = FALSE,
                                         diagnosis = FALSE,
                                         niter = 20,
                                         tol = .Machine$double.eps^0.5) {

  loss.function <- match.arg(loss.function, c("huber", "tukey"))
  initialization <- match.arg(initialization, c("l2", "mode"))

  rho <- switch(
    loss.function,
    huber = function(r, ...) mr_raps_rho_huber(r, k, ...),
    tukey = function(r, ...) mr_raps_rho_tukey(r, k, ...)
  )

  delta <- stats::integrate(function(x) x * rho(x, deriv = 1) * stats::dnorm(x), -Inf, Inf)$value
  c1 <- stats::integrate(function(x) rho(x, deriv = 1)^2 * stats::dnorm(x), -Inf, Inf)$value
  c2 <- stats::integrate(function(x) x^2 * rho(x, deriv = 1)^2 * stats::dnorm(x), -Inf, Inf)$value - delta^2
  c3 <- stats::integrate(function(x) x^2 * rho(x, deriv = 2) * stats::dnorm(x), -Inf, Inf)$value

  robust.E <- function(beta, tau2) {
    t <- (b_out - beta * b_exp) / sqrt(tau2 + se_out^2 + se_exp^2 * beta^2)
    se_exp^2 * (t * rho(t, deriv = 1) - delta) / (tau2 + se_out^2 + se_exp^2 * beta^2)
  }

  robust.loglike.fixtau <- function(beta, tau2) {
    alpha.hat <- 0
    -0.5 * sum(rho((b_out - alpha.hat - b_exp * beta) / sqrt(tau2 + se_out^2 + se_exp^2 * beta^2)))
  }

  bound.beta <- stats::quantile(abs(b_out / b_exp), 0.95) * 10
  bound.tau2 <- stats::quantile(se_out^2, 0.95) * 10

  if (initialization == "mode") {
    stop("Initialization by mode estimator is currently not supported.", call. = FALSE)
  } else {
    fit <- mr_raps_overdispersed(b_exp, b_out, se_exp, se_out, suppress.warning = TRUE)
    beta.hat <- fit$beta.hat
    tau2.hat <- fit$tau2.hat
  }

  for (iter in 1:niter) {
    beta.hat.old <- beta.hat
    tau2.hat.old <- tau2.hat

    tau2.hat <- tryCatch(
      stats::uniroot(
        function(tau2) sum(robust.E(beta.hat, tau2)),
        bound.tau2 * c(0, 1),
        extendInt = "yes",
        tol = bound.tau2 * .Machine$double.eps^0.25
      )$root,
      error = function(e) {
        warning("Did not find a solution for tau2; setting tau2=0.", call. = FALSE)
        0
      }
    )
    if (tau2.hat < 0) tau2.hat <- 0
    if (tau2.hat > bound.tau2 * 0.95) warning("Estimated overdispersion seems abnormally large.", call. = FALSE)

    beta.hat <- stats::optim(
      beta.hat,
      function(beta) robust.loglike.fixtau(beta, tau2.hat),
      method = "L-BFGS-B",
      lower = -bound.beta,
      upper = bound.beta,
      control = list(fnscale = -1)
    )$par

    if (abs(beta.hat.old - beta.hat) / abs(beta.hat + 1e-10) +
        abs(tau2.hat.old - tau2.hat) / abs(tau2.hat + 1e-10) <= tol) {
      break
    }
  }

  if ((tau2.hat <= min(se_out^2) / 5) && (!suppress.warning)) {
    warning("Estimated overdispersion is very small; consider no-overdispersion model.", call. = FALSE)
  }

  score.var <- diag(c(
    c1 * sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) +
                (b_out^2 - tau2.hat - se_out^2) * se_exp^2 +
                se_exp^2 * (tau2.hat + se_out^2)) / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
    (c2 / 2) * sum(2 * se_exp^4 / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2)
  ))

  I <- matrix(c(
    -delta * sum(((b_exp^2 - se_exp^2) * (tau2.hat + se_out^2) +
                    (b_out^2 - tau2.hat - se_out^2) * se_exp^2) / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
    0,
    -delta * sum(se_exp^2 * beta.hat / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2),
    -(delta + c3) / 2 * sum(se_exp^2 / (tau2.hat + se_out^2 + se_exp^2 * beta.hat^2)^2)
  ), 2, 2)

  asymp.var <- solve(I) %*% score.var %*% t(solve(I))

  out <- list(
    beta.hat = beta.hat,
    tau2.hat = tau2.hat,
    beta.se = sqrt(asymp.var[1, 1]),
    tau2.se = sqrt(asymp.var[2, 2]),
    beta.p.value = min(1, 2 * (1 - stats::pnorm(abs(beta.hat) / sqrt(asymp.var[1, 1]))))
  )

  if (diagnosis) {
    std.resid <- (b_out - b_exp * beta.hat) / sqrt(tau2.hat + se_out^2 + beta.hat^2 * se_exp^2)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = c(1, 2))
    stats::qqnorm(std.resid); graphics::abline(0, 1)
    if (requireNamespace("nortest", quietly = TRUE)) {
      ad <- nortest::ad.test(std.resid)
      message("Anderson-Darling test: statistic = ", round(ad$statistic, 4),
              ", p-value = ", round(ad$p.value, 4))
    }
    sw <- stats::shapiro.test(std.resid)
    message("Shapiro-Wilk test: statistic = ", round(sw$statistic, 4),
            ", p-value = ", round(sw$p.value, 4))
    out$std.resid <- std.resid
  }

  out
}
