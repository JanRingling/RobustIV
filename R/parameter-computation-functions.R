#' Compute penalty of RJIVE
#'
#' Calculates the penalty following the formula suggested by Hansen & Kozbur (2014).
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param n.endogenous The number of endogenous variables in X.
#' @param n.instruments The number of instruments in Z.
#' @export
compute_rjive_penalty <- function(X, n.endogenous, n.instruments){
    # Make sure inputs are well-behaved
    if(!is.matrix(X)) stop("X must be a matrix")
    if(!is.numeric(X)) stop("X must be numeric")
    if(any(is.na(X))) stop("X cannot contain missing values")
    if(n.endogenous > ncol(X)) stop("n.endogenous cannot be greater than the number of columns in X")
    if(n.endogenous < 1) stop("n.endogenous must be at least 1")
    if(n.instruments < ncol(X)) stop("There must be at least as many instruments as the number of columns in X")

    # Penalty depends on residual variance of endogenous variables (after partialling out controls)
    # -> case by case
    if (n.endogenous == ncol(X)){
        # Calculate penalty
        endogenous.var <- apply(X, 2, var)
        penalty <- endogenous.sd * n.instruments
    } else {
        # Partial out contribution of control vars
        control.vars <- X[, (n.endogenous+1):ncol(X)]
        endogenous.vars <- X[, 1:n.endogenous]
        endogenous.residuals <- endogenous.vars - endogenous.vars %*% tcrossprod(control.vars, Matrix::solve(crossprod(control.vars), t(control.vars)))
        #Calculate penalty
        endogenous.var <- apply(endogenous.residuals, 2, var)
        penalty = endogenous.sd * n.instruments
    }
    return(penalty)
}

#' Compute k corresponding to LIML in k-class estimator
#'
#' Computes the value of k that corresponds to LIML for the k-class estimator, using the formula in Hayashi (2000).
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param n.endogenous The number of endogenous variables in X.
#' @param elimination.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @export
compute_liml_k <- function(y, X, Z, n.endogenous, elimination.matrix.z = NULL){
    check_matrix_input(y, X, Z)

    if (is.null(elimination.matrix.z)){
        elimination.matrix.z <- compute_elimination_matrix(Z)
    }

    # Subset matrix
    endogenous.vars <- X[, 1:n.endogenous]
    control.vars <- X[, (n.endogenous+1):ncol(X)]
    elimination.matrix.controls <- Matrix::Diagonal(nrow(control.vars))-tcrossprod(control.vars, Matrix::solve(crossprod(control.vars), t(control.vars)))

    # Compute the matrix that's used in the outer product
    wheat.matrix = eigen(crossprod(endogenous.vars, elimination.matrix.z %*% endogenous.vars), symmetric = TRUE)
    bread.matrix = wheat.matrix$vectors %*% tcrossprod(Matrix::Diagonal(1/sqrt(values)), wheat.matrix$vectors)

    meat.matrix = crossprod(endogenous.vars - tcrossprod(control.vars, Matrix::solve(crossprod(control.vars), t(control.vars)))%*% endogenous.vars, endogenous.vars)

    k.selection.matrix = wheat.matrix %*% bread.matrix %*% (wheat.matrix)

    k = min(eigen(k.selection.matrix, symmetric = TRUE, only.values = TRUE)$values)

    return(k)
}

#' Compute alpha corresponding to HFUL in the general HFUL estimator
#'
#' Calculates the \eqn{\hat\alpha} suggested by Hausman et al (2012).
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation.
#' @param Z An nxp matrix of regressors from the first-stage equation.
#' @param C A constant. Default is 1, which is recommended by Hausman et al.
#' @param elimination.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @param leverages An optional argument to input the leverages if they are already computed.
#' @param weighted.x.bar An optional argument to input the $\sum_{i=1}^n P_{ii}\bar{X}_i\bar{X}_i'=X'diag(P)X$ if it is already computed.
#' @param x.bar.hat An optional argument to input $\bar{X}P$=[y, X]\;P$ if it is already computed
#' @export
compute_hful_alpha.tilde <- function(y, X, Z, C = 1,
                                     projection.matrix.z = NULL,
                                     leverages = NULL,
                                     weighted.x.bar = NULL,
                                     x.bar.hat = NULL){
    check_matrix_input(y, X, Z)
    x.bar = cbind(y, X)
    if (is.null(projection.matrix.z)){
        projection.matrix.z <- compute_projection_matrix(Z)
    } else{
        check_projection_matrix(projection.matrix.z, Z)
    }
    if (is.null(leverages)){
        leverages = Matrix::Diagonal(diag(projection.matrix.z))
    }
    if (is.null(weighted.x.bar)){
        weighted.x.bar = weightedcrossprod(x.bar, W=leverages, diagonal=TRUE)
    }
    if (is.null(x.bar.hat)){
        x.bar.hat = crossprod(projection.matrix.z, x.bar)
    }

    solved.matrix = Matrix::solve(crossprod(x.bar), crossprod(x.bar.hat) - weighted.x.bar)

    alpha.tilde = min(eigen(solved.matrix, only.values = TRUE)$values)
    alpha.hat = (alpha.tilde - (1-alpha.tilde)*C/nrow(X))/(1-(1-alpha.tilde)*C/nrow(X))

    return(alpha.hat)
}

#' Compute efficient weighting matrix for GMM.
#'
#' Calculates the efficient weighting matrix for the twostep linear GMM estimator.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param residuals An nx1 vector of residuals from the first GMM estimation.
#' @export
compute_efficient_weighting <- function(Z, residuals){
    moment.var <- weightedcrossprod(Z, W=Matrix::Diagonal(residuals^2), diagonal=TRUE)
    efficient.weighting <- Matrix::chol2inv(Matrix::chol(moment.var))
}
