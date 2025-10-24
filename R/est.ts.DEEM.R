#' Calculate the DEEM estimator and its estimated standard error in the two-sample setting
#'
#' @description
#' This function implements the Debiased Estimating Equation Method (DEEM) for Mendelian
#' randomization in the two-sample setting, producing consistent causal effect estimates even when instruments are weak
#' and/or correlated.
#' Inputs of this function are the outputs of the function "prepare.DEEM".
#' Examples for the required summary statistics can be downloaded from https://osf.io/j25rc/files/osfstorage
#' (the LD matrix is approximately 6GB after decompression).
#'
#' @param index Number of LD blocks included.
#' @param Gamma_h_lst A list of vectors of estimated SNP–outcome associations within different LD blocks from the outcome sample.
#' @param gamma_h_lst A list of vectors of estimated SNP–exposure associations within different LD blocks from the exposure sample.
#' @param gamma_s_lst A list of vectors of estimated SNP–exposure associations within different LD blocks from the supplemental exposure sample (used to select SNPs).
#' @param Sig_y A list of estimated variance-covariance matrix for the estimated SNP–outcome associations within different LD blocks from the outcome sample.
#' @param Sig_x A list of estimated variance-covariance matrix for the estimated SNP–exposure associations within different LD blocks from the exposure sample.
#' @param Sig_alpha A list of estimated variance-covariance matrix for the pleiotropic effects within different LD blocks.
#' @param R_sel A list of estimated LD matrix within different LD blocks.
#'
#' @return
#' A vector containing:
#' \itemize{
#'   \item \code{beta_hat}: the estimated causal effect.
#'   \item \code{se}: the estimated standard error.
#' }
#'
#' @export
est.ts.DEEM <- function(index, Gamma_h_lst, gamma_h_lst, gamma_s_lst, Sig_y, Sig_x, Sig_alpha, R_sel) {
  D_x <- list()
  D_y <- list()
  D_alpha <- list()

  for (i in 1:index) {
    if (length(gamma_s_lst[[i]]) == 1) {
      D_x[[i]] <- Sig_x[[i]]
      D_y[[i]] <- Sig_y[[i]]
      D_alpha[[i]] <- Sig_alpha[[i]]
    } else {
      D_x[[i]] <- diag(diag(Sig_x[[i]]))
      D_y[[i]] <- diag(diag(Sig_y[[i]]))
      D_alpha[[i]] <- diag(diag(Sig_alpha[[i]]))
    }
  }

  q <- 0
  for (i in 1:index) {
    q <- q + length(gamma_s_lst[[i]])
  }

  W_or <- list()

  for (i in 1:index) {
    W_or[[i]] <- D_y[[i]]
    W_or[[i]] <- solve(sqrt(W_or[[i]]) %*% (1 * diag(1, length(gamma_s_lst[[i]])) + 1 * R_sel[[i]]) %*% sqrt(W_or[[i]]))
  }


  ####################################calculate initial estimator using optimal weighting
  num <- 0
  den <- 0

  for (i in 1:index) {
    num <- num + t(gamma_s_lst[[i]]) %*% W_or[[i]] %*% Gamma_h_lst[[i]]
    den <- den + t(gamma_s_lst[[i]]) %*% W_or[[i]] %*% gamma_h_lst[[i]]
  }
  init <- c(num / den)

  D_x <- list()
  D_y <- list()
  D_alpha <- list()
  for (i in 1:index) {
    if (length(gamma_s_lst[[i]]) == 1) {
      D_x[[i]] <- Sig_x[[i]]
      D_y[[i]] <- Sig_y[[i]]
      D_alpha[[i]] <- Sig_alpha[[i]]
    } else {
      D_x[[i]] <- diag(diag(W_or[[i]] %*% Sig_x[[i]]) / diag(W_or[[i]]))
      D_y[[i]] <- diag(diag(W_or[[i]] %*% Sig_y[[i]]) / diag(W_or[[i]]))
      D_alpha[[i]] <- diag(diag(W_or[[i]] %*% Sig_alpha[[i]]) / diag(W_or[[i]]))
    }
  }

  r_tilde <- list()
  for (i in 1:index) {
    r_tilde[[i]] <- Gamma_h_lst[[i]] - gamma_h_lst[[i]] * init
  }

  W_tau <- list()
  for (i in 1:index) {
    W_tau[[i]] <- W_or[[i]]
  }

  num <- 0
  den <- 0

  for (i in 1:index) {
    num <- num + t(r_tilde[[i]]) %*% W_tau[[i]] %*% r_tilde[[i]] - sum(diag(W_tau[[i]] %*% (D_y[[i]] + init^2 * D_x[[i]])))
    den <- den + sum(diag(W_tau[[i]] %*% D_alpha[[i]]))
  }

  tau2_h <- c(num / den)

  D_x.inv <- list()
  D_r.inv <- list()
  Sig_r <- list()
  for (i in 1:index) {
    D_x.inv[[i]] <- solve(D_x[[i]])
    D_r.inv[[i]] <- solve(D_y[[i]] + init^2 * D_x[[i]] + max(tau2_h, 0) * D_alpha[[i]])
    Sig_r[[i]] <- Sig_y[[i]] + init^2  * Sig_x[[i]] + max(tau2_h, 0) * Sig_alpha[[i]]
  }

  num <- 0
  den <- 0

  for (i in 1:index) {
    num <- num + t(gamma_s_lst[[i]]) %*% W_or[[i]] %*% Gamma_h_lst[[i]]
    den <- den + t(gamma_s_lst[[i]]) %*% W_or[[i]] %*% gamma_h_lst[[i]]
  }

  orDEEM3 <- c(num / den)


  ################################################################calculate matrix required to debias with non-diagonal weighting matrix

  orEE <- function(beta) {
    g <- 0
    for (i in 1:index) {
      g <- g + t(gamma_h_lst[[i]] + beta * D_x[[i]] %*% solve(D_y[[i]] + tau2_h * D_alpha[[i]] + beta^2 * D_x[[i]]) %*% (Gamma_h_lst[[i]] -
                                                                                                                           beta * gamma_h_lst[[i]])) %*% W_or[[i]] %*% (Gamma_h_lst[[i]] - beta * gamma_h_lst[[i]])  #computation can be boosted by modifying this formulation
    }
    g^2
  }

  orDEEM <- nlminb(orDEEM3, orEE, lower = orDEEM3 - 0.4, upper = orDEEM3 + 0.4)$par

  ##########################################################calculate meta estimator

  Sig_JT <- list()

  for (i in 1:index) {
    M_tmp1 <- cbind(Sig_x[[i]], - init * Sig_x[[i]])
    M_tmp2 <- cbind(- init * Sig_x[[i]], Sig_r[[i]])
    Sig_JT[[i]] <- rbind(M_tmp1, M_tmp2)
  }

  M1 <- list()

  for (i in 1:index) {
    qi <- length(gamma_s_lst[[i]])
    M_tmp1 <- cbind(matrix(0, qi, qi), W_or[[i]])
    M_tmp2 <- cbind(W_or[[i]], init*(D_r.inv[[i]] %*% D_x[[i]] %*% W_or[[i]] + W_or[[i]] %*% D_x[[i]] %*% D_r.inv[[i]]))
    M1[[i]] <- 1 / 2 * rbind(M_tmp1, M_tmp2) %*% Sig_JT[[i]]
  }

  M_rho <- list()
  for (i in 1:index) {
    qi <- length(gamma_s_lst[[i]])
    M_tmp1 <- cbind(matrix(0, qi, qi), W_or[[i]])
    M_tmp2 <- cbind(W_or[[i]], matrix(0, qi, qi))
    M_rho[[i]] <- 1 / 2 * rbind(M_tmp1, M_tmp2) %*% Sig_JT[[i]]
  }

  M_tau <- list()
  for (i in 1:index) {
    qi <- length(gamma_s_lst[[i]])
    M_tmp1 <- cbind(matrix(0, qi, qi), matrix(0, qi, qi))
    M_tmp2 <- cbind(- init * W_tau[[i]] %*% Sig_x[[i]], W_tau[[i]] %*% Sig_r[[i]])
    M_tau[[i]] <- rbind(M_tmp1, M_tmp2)
  }

  m1 <- list()
  for (i in 1:index) {
    m1[[i]] <- W_or[[i]] %*% gamma_h_lst[[i]]
  }

  m2 <- list()
  for (i in 1:index) {
    m2[[i]] <- W_or[[i]] %*% gamma_s_lst[[i]]
  }

  m_rho <- m1

  Cov_s <- matrix(0, 4, 4)

  for (i in 1:index) {
    M_tmp <- matrix(0, 4, 4)

    M_tmp[1, 1] <- t(m1[[i]]) %*% Sig_r[[i]] %*% m1[[i]] -
      sum(diag(W_or[[i]] %*% Sig_r[[i]] %*% W_or[[i]] %*% Sig_x[[i]])) +
      2 * sum(as.vector(t(M1[[i]])) * as.vector(M1[[i]]))

    M_tmp[1, 2] <- M_tmp[2, 1] <- t(m1[[i]]) %*% Sig_r[[i]] %*% m2[[i]]

    M_tmp[1, 3] <- M_tmp[3, 1] <-  t(m1[[i]]) %*% Sig_r[[i]] %*% m_rho[[i]] -
      sum(diag(W_or[[i]] %*% Sig_r[[i]] %*% W_or[[i]] %*% Sig_x[[i]])) +
      2 * sum(as.vector(t(M1[[i]])) * as.vector(M_rho[[i]]))

    M_tmp[1, 4] <- M_tmp[4, 1] <- 2 * sum(as.vector(t(M1[[i]])) * as.vector(M_tau[[i]]))

    M_tmp[2, 2] <- t(m2[[i]]) %*% Sig_r[[i]] %*% m2[[i]]

    M_tmp[2, 3] <- M_tmp[3, 2] <- t(m_rho[[i]]) %*% Sig_r[[i]] %*% m2[[i]]

    M_tmp[2, 4] <- M_tmp[4, 2] <- 0

    M_tmp[3, 3] <- t(m_rho[[i]]) %*% Sig_r[[i]] %*% m_rho[[i]] -
      sum(diag(W_or[[i]] %*% Sig_r[[i]] %*% W_or[[i]] %*% Sig_x[[i]])) +
      2 * sum(as.vector(t(M_rho[[i]])) * as.vector(M_rho[[i]]))

    M_tmp[3, 4] <- M_tmp[4, 3] <- 2 * sum(as.vector(t(M_rho[[i]])) * as.vector(M_tau[[i]]))

    M_tmp[4, 4] <- 2 * sum(as.vector(t(M_tau[[i]])) * as.vector(M_tau[[i]]))

    Cov_s <- Cov_s + M_tmp
  }

  a1 <- rep(0, 4)
  a1[1] <- 1

  num <- 0
  den <- 0
  for (i in 1:index) {
    num <- num + sum(diag(W_or[[i]] %*% D_x[[i]] %*% D_r.inv[[i]] %*% D_alpha[[i]]))
    den <- den + sum(diag(W_tau[[i]] %*% D_alpha[[i]]))
  }

  a1[4] <- - init * c(num / den)

  C <- cbind(a1, c(0, 1, 0, 0))

  G <- rep(0, 2)
  for (i in 1:index) {
    G <- G + c(- t(gamma_h_lst[[i]]) %*% W_or[[i]] %*% gamma_h_lst[[i]] + sum(diag(W_or[[i]] %*% D_x[[i]])),
               - t(gamma_s_lst[[i]]) %*% W_or[[i]] %*% gamma_h_lst[[i]])
  }

  v_orDEEM <- sqrt(t(C[, 1]) %*% Cov_s %*% C[, 1] / G[1]^2)
  v_orDEEM3 <- sqrt(t(C[, 2]) %*% Cov_s %*% C[, 2] / G[2]^2)

  V_meta <- t(C) %*% Cov_s %*% C / (G %*% t(G))
  ev <- eigen(V_meta)$values
  V_meta <- V_meta + 0.5 * (min(diag(V_meta)) / max(diag(V_meta)) - ev[2] / ev[1]) * diag(diag(V_meta))
  est <- c(orDEEM, orDEEM3)

  W_meta <- solve(V_meta)
  w_meta <- colSums(W_meta) / sum(W_meta)

  DEEM_meta <- sum(w_meta %*% est)
  v_meta <- sqrt(t(w_meta) %*% V_meta %*% w_meta)

  c(DEEM_meta, v_meta)
}
