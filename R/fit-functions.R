#' Fit 2SLS estimator
#'
#' Fits the 2SLS estimator given the outcomes, instruments and structural regressors.
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @export
tsls.fit <- function(y, X, Z){
    check_matrix_input(y, X, Z)

    projection.matrix.z <- compute_projection_matrix(Z)
    x.hat <- crossprod(projection.matrix.z, X)
    result <- Matrix::solve(crossprod(x.hat, X), crossprod(x.hat, y))
    names(result) <- colnames(X)
    return(result)
}

#' Fit JIVE1 estimator
#'
#' Fits the JIVE1 estimator given the outcomes, instruments and structural regressors.
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @export
jive1.fit <- function(y, X, Z){
    check_matrix_input(y, X, Z)

    projection.matrix.z <- compute_projection_matrix(Z)
    leverages <- tr(projection.matrix.z)
    x.hat <- (projection.matrix.z %*% X - leverages * X) / (1 - leverages)

    result <- Matrix::solve(crossprod(x.hat, x), crossprod(x.hat, y))
    names(result) <- colnames(X)
    return(result)
}

#' Fit JIVE2 estimator
#'
#' Fits the JIVE2 estimator given the outcomes, instruments and structural regressors.
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @export
jive2.fit <- function(y, X, Z){
    check_matrix_input(y, X, Z)

    projection.matrix.z <- compute_projection_matrix(Z)
    leverages <- tr(projection.matrix.z)
    x.hat <- (projection.matrix.z %*% X - leverages * X) / (1 - 1/nrow(X))

    result <- Matrix::solve(crossprod(x.hat, x), crossprod(x.hat, y))
    names(result) <- colnames(X)
    return(result)
}

#' Fit k-classs estimator given k
#'
#'  Fits the k-class estimator given a formula and a dataset. Given an appropriate choice of k, this function also fits the LIML estimator.
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @param k A constant to be used for the k-class estimator
#' @param elimination.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @export
kclass.fit <- function(y, X, Z, k,
                       elimination.matrix.z = NULL){
    check_matrix_input(y, X, Z)
    # Check k
    if(!is.numeric(k)) stop("k must be numeric")
    n.obs = nrow(X)
    if (is.null(elimination.matrix.z)){
        elimination.matrix.z <- compute_elimination_matrix(Z)
    }
    k.term <- crossprod(X, Matrix::Diagonal(n.obs-k*elimination.matrix.z))

    result <- Matrix::solve(k.term %*% X, k.term %*% y)
    # Add names
    names(result) <- colnames(X)
    return(result)
}

#' Fit RJIVE estimator given a penalty
#'
#' Fits the RJIVE estimator.
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation.
#' @param penalty A constant or vector > 0 to be used as Ridge penalty. Must be a vector for multiple endogenous variables. Corresponds to \gamma in the original notation.
#' @param n.endogenous The number of endogenous variables in X.
#' @export
rjive.fit <- function(y, X, Z, penalty, n.endogenous){
    check_matrix_input(y, X, Z)
    # Check penalty
    if(!is.numeric(penalty)) stop("penalty must be numeric")
    if(penalty <= 0) stop("penalty must be positive")
    n.obs <- nrow(X)
    n.instruments <- ncol(Z)
    n.structural <- ncol(X)

    var.cov.z <- crossprod(Z)
    var.cov.z.inv <- Matrix::solve(var.cov.z)
    # Later use eigen decomposition to speed up inversion.
    eigenval.decomp <- eigen(var.cov.z)
    d.matrix = Matrix::Diagonal(eigenval.decomp$values)
    p.inverse = Matrix::solve(eigenval.decomp$vectors)

    # Cycle through the endogenous variables as first stage regressions
    x.hat.matrix = matrix(NA, nrow = n.obs, ncol = n.endogenous)
    for (i in 1:n.endogenous){
        # Current endogenous var
        current.x <- X[,c(i,(n.endogenous+1):n.structural)]
        # Compute penalty matrix
        penalty.matrix <- Matrix::Diagonal(penalty[i],nrow = n.instruments)
        # Compute penalised inverse quickly using eigenvalue decomposition
        penalized.var.cov.z.inv <- eigenval.decomp$vectors %*% Matrix::solve(d.matrix + penalty.matrix, p.inverse)
        # Cycle through observations to compute jackknife estimates
        for (j in 1:n.obs){
            x.hat.matrix[j,i] = Z %*% (sherman.morrison(penalized.var.cov.z.inv, -Z[j,], t(Z[j,]), inverted = TRUE) %*% (crossprod(Z, current.x) - tcrossprod(Z[j,], current.x[j,])))
        }
    }
    instrumented.matrix = cbind(x.hat.matrix, X[,(n.endogenous+1):n.structural])
    result = Matrix::solve(crossprod(x.hat.matrix, X), crossprod(x.hat.matrix, y))
    # Add names
    names(result) <- colnames(X)
    return(result)
}

#' HFUL formula given alpha hat in the general HFUL estimator
#'
#' Computes  suggested by Hausman et al (2012).
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation.
#' @param Z An nxp matrix of regressors from the first-stage equation.
#' @param alpha.hat The alpha hat computed according to the HFUL procedure.
#' @param projection.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @param leverages An optional argument to input the leverages if they are already computed.
#' @param weighted.x.bar An optional argument to input the $\sum_{i=1}^n P_{ii}\bar{X}_i\bar{X}_i'=X'diag(P)X$ if it is already computed.
#' @param x.bar.hat An optional argument to input $\bar{X}P$=[y, X]\;P$ if it is already computed
#' @export
hful.fit <- function(y, X, Z, alpha.hat,
                     elimination.matrix.z = NULL,
                     leverages = NULL,
                     weighted.x.bar = NULL,
                     x.bar.hat = NULL){
    check_matrix_input(y, X, Z)
    if (is.null(projection.matrix.z)){
        projection.matrix.z <- compute_projection_matrix(Z)
    } else{
        check_input_matrix(projection.matrix.z, "projection.matrix.z")
    }
    if (is.null(leverages)){
        leverages = Matrix::Diagonal(diag(projection.matrix.z))
    }
    if (is.null(weighted.x.bar)){
        weighted.x = crossprod(X, leverages %*% X)
    } else{
        weighted.x = weighted.x.bar[2:ncol(weighted.x.bar),2:ncol(weighted.x.bar)]
    }
    if (is.null(x.bar.hat)){
        x.hat = crossprod(projection.matrix.z, X)
    } else{
        x.hat = x.bar.hat[,2:ncol(x.bar.hat)]
    }

    delta.hat = Matrix::solve(crossprod(x.hat, X- weighted.x - alpha.hat *  crossprod(X)), crossprod(x.hat, y)-crossprod(X, leverages %*% y)-alpha.hat * crossprod(X, y))
    # Add names
    names(delta.hat) <- colnames(X)
    return(delta.hat)
}
#' Fit GMM moments estimator given weighting matrix W
#'
#'  Fits the linear GMM estimator given the outcomes, instruments, structural regressors and a weighting matrix.
#' @param y An nx1 vector of outcomes
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @param W A p.s.d. weighting matrix. If no argument is supplied, \eqn{(Z'Z)^{-1}} is used corresponding to 2SLS. Has no effect if the parameters are just-identified.
#' @param cov.zx An optional argument to input \eqn{Z'X} if it is already computed.
#' @export
gmm.fit <- function(y, X, Z, W = NULL, cov.zx = NULL){
    check_matrix_input(y, X, Z)
    if (is.null(cov.zx)){
        cov.zx <- crossprod(Z, X)
    } else{
        check_input_matrix(cov.zx, "cov.zx")
    }
    if (ncol(X)==ncol(Z)){
        # Just identified case
        result <- Matrix::solve(cov.zx, crossprod(Z, y))
    } else{
        # Overidentified case
        if (is.null(W)){
            W <- Matrix::solve(crossprod(Z))
        } else{
            check_individual_matrix(W, "Weighing Matrix W")
        }
        cov.zy <- crossprod(Z, y)
        weighting.matrix.times.cov.zx <- crossprod(W, cov.zx)
        result <- Matrix::solve(crossprod(weighting.matrix.times.cov.zx, cov.zx), crossprod(weighting.matrix.times.cov.zx, cov.zy))
    }
    # Bind names to result vector
    names(result) <- colnames(X)
    return(result)
}

