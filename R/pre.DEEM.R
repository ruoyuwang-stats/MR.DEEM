#' Prepare the data required by DEEM, outputs of this function are the inputs of "est.ts.DEEM" and "est.os.DEEM"
#'
#' @description
#' The function prepares the data required by DEEM, outputs of this function are the inputs of "est.ts.DEEM"
#' and "est.os.DEEM"
#' Examples for the required summary statistics can be downloaded from https://osf.io/j25rc/files/osfstorage
#' (the LD matrix is approximately 6GB after decompression)
#'
#' @param SumData A data.frame consists of eight variables:
#' \itemize{
#'   \item \code{chr}: the chromosome.
#'   \item \code{pos}: the position.
#'   \item \code{gamma_s}: a vector of estimated SNP–exposure associations from the supplemental exposure sample used to select the SNPs.
#'   \item \code{sigma_s}: a vector of estimated standard error for the estimated SNP–exposure associations from the supplemental exposure sample used to select the SNPs.
#'   \item \code{gamma_h}: a vector of estimated SNP–exposure associations from the exposure sample.
#'   \item \code{sigma_x}: a vector of estimated standard error for the estimated SNP–exposure associations from the exposure sample.
#'   \item \code{Gamma_h}: a vector of estimated SNP–outcome associations from the outcome sample.
#'   \item \code{sigma_y}: a vector of estimated standard error for the estimated SNP–outcome associations from the outcome sample.
#' }
#' SNPs in the three samples should be the same (possible preprocessing required).
#' @param com_list A vector of characters in the form "chr:pos" of the SNPs in the SumData.
#' @param LD_Ref A list of LD matrix in each LD blocks from the reference sample.
#' @param Block data.frame consists of three variables:
#' \itemize{
#'   \item \code{chr}: a vector of character in the form "chr?" where "?" is the chromosome that the LD block locates on.
#'   \item \code{start}: a vector of int indicates where the LD block starts.
#'   \item \code{stop}: a vector of int indicates where the LD block ends.
#' }
#' @param ns sample size of the supplemental exposure sample that generates gamma_s.
#' @param ne sample size of the exposure sample that generates gamma_h.
#' @param no sample size of the outcome sample that generates Gamma_h.
#' @param p_thr p-value threshold, default 0.1.
#' @param r2_thr r2 threshold, default 0.81.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{index}: Number of LD blocks included.
#'   \item \code{Gamma_h_lst}: A list of vectors of estimated SNP–outcome associations within different LD blocks from the outcome sample.
#'   \item \code{gamma_h_lst}: A list of vectors of estimated SNP–exposure associations within different LD blocks from the exposure sample.
#'   \item \code{gamma_s_lst}: A list of vectors of estimated SNP–exposure associations within different LD blocks from the supplemental exposure sample (used to select SNPs).
#'   \item \code{Sig_y}: A list of estimated variance-covariance matrix for the estimated SNP–outcome associations within different LD blocks from the outcome sample.
#'   \item \code{Sig_x}: A list of estimated variance-covariance matrix for the estimated SNP–exposure associations within different LD blocks from the exposure sample.
#'   \item \code{Sig_alpha}: A list of estimated variance-covariance matrix for the pleiotropic effects within different LD blocks.
#'   \item \code{R_sel}: A list of estimated LD matrix within different LD blocks.
#' }
#'
#' @export
prepare.DEEM <- function(SumData, com_list, LD_Ref, Blocks, ns, ne, no, p_thr = 0.1, r2_thr = 0.81) {
  dir_list <- matrix(as.numeric(unlist(strsplit(com_list, ":"))), 2)
  dir <- dir_list[1, ]
  sub_dir <- list()
  for (j in 1:22) {
    sub_dir[[j]] <- dir_list[2, ][dir == j]
  }
  Blocks$chr <- as.numeric(gsub("chr", "", Blocks$chr))
  nLDB <- nrow(Blocks)

  for (i in 1:length(LD_Ref)) {  ##shrink estimator to avoid negative eigenvalues
    for (j in 1:length(LD_Ref[[i]])) {
      if (length(LD_Ref[[i]][[j]]) == 0) {
        next
      } else if (length(LD_Ref[[i]][[j]]) == 1) {
        LD_Ref[[i]][[j]] <- 1
      } else {
        LD_Ref[[i]][[j]] <- (1 - 1e-3) * LD_Ref[[i]][[j]] + diag(1e-3, nrow(LD_Ref[[i]][[j]]))
      }
    }
  }
  sigx_max <- median(SumData$sigma_x) + 10 * sd(SumData$sigma_x)
  sigy_max <- median(SumData$sigma_y) + 10 * sd(SumData$sigma_y)

  gamma_s_lst <- list()
  gamma_h_lst <- list()
  Gamma_h_lst <- list()
  Sig_x <- list()
  Sig_y <- list()
  Sig_alpha <- list()
  R_sel <- list()
  index <- 1
  j <- 0
  k_j <- 0
  for (k in 1:nLDB) {
    k_j <- k_j + 1
    if (j != Blocks$chr[k]) {
      j <- Blocks$chr[k]
      dat <- SumData[SumData$chr == j, ]
      gamma_s <- dat$gamma_s
      sigma_s <- dat$sigma_s
      gamma_h <- dat$gamma_h
      sigma_x <- dat$sigma_x
      Gamma_h <- dat$Gamma_h
      sigma_y <- dat$sigma_y
      k_j <- 1  ##kj-th block in chr j
    }
    LDB <- Blocks$start[k] <= sub_dir[[j]] & sub_dir[[j]] < Blocks$stop[k]
    if (sum(LDB) == 0) {
      next
    }
    gamma_s_B <- gamma_s[LDB]
    gamma_h_B <- gamma_h[LDB]
    Gamma_h_B <- Gamma_h[LDB]
    sigma_s_B <- sigma_s[LDB]
    sigma_x_B <- sigma_x[LDB]
    sigma_y_B <- sigma_y[LDB]

    R <- LD_Ref[[j]][[k_j]]
    p_k <- nrow(R)

    order <- order(abs(gamma_s_B / sigma_s_B))
    C <- R - diag(diag(R))
    sel <- rep(T, p_k)

    pval <- 2 * (1 - pnorm(abs(gamma_s_B / sigma_s_B)))

    for (i in order) {
      sel[i] <- max(abs(C[i, sel])) < sqrt(r2_thr)
    }

    sel <- sel & pval < p_thr & sigma_x_B < sigx_max & sigma_y_B < sigy_max
    if (sum(sel) < 1) {
      next
    }

    q[2] <- q[2] + sum(sel)
    Phi_y <- as.vector(sqrt(sigma_y_B^2 + Gamma_h_B^2 / (no - 2)))
    M_tmp <- Phi_y^{-1} * R
    M_tmp <- Phi_y * t(M_tmp)
    R_sel[[index]] <- R[sel, sel]

    gamma_s_B <- gamma_s_B[sel]
    gamma_h_B <- gamma_h_B[sel]
    Gamma_h_B <- Gamma_h_B[sel]
    sigma_s_B <- sigma_s_B[sel]
    sigma_x_B <- sigma_x_B[sel]
    sigma_y_B <- sigma_y_B[sel]

    Phi_x <- as.vector(sqrt(sigma_x_B^2 + gamma_h_B^2 / (ne - 2)))
    Phi_y <- as.vector(sqrt(sigma_y_B^2 + Gamma_h_B^2 / (no - 2)))

    if (sum(sel) == 1) {
      Sig_x[[index]] <- Phi_x %*% R_sel[[index]] %*% Phi_x
      Sig_y[[index]] <- Phi_y %*% R_sel[[index]] %*% Phi_y
      Sig_alpha[[index]] <- sum(M_tmp[sel, ]^2)

    } else {
      Sig_x[[index]] <- diag(Phi_x) %*% R_sel[[index]] %*% diag(Phi_x)
      Sig_y[[index]] <- diag(Phi_y) %*% R_sel[[index]] %*% diag(Phi_y)
      Sig_alpha[[index]] <- M_tmp[sel, ] %*% t(M_tmp[sel, ])
    }
    gamma_s_lst[[index]] <- gamma_s_B
    gamma_h_lst[[index]] <- gamma_h_B
    Gamma_h_lst[[index]] <- Gamma_h_B
    index <- index + 1
  }

  index <- index - 1

  return(list(index, Gamma_h_lst, gamma_h_lst, gamma_s_lst, Sig_y, Sig_x, Sig_alpha, R_sel))
}
