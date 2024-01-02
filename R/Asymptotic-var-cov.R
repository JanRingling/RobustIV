#' Calculate 2SLS asymptotic var-cov matrix
#'
#' Calculates the 2SLS var-cov matrix given asymptotic theory
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param type The type of var-cov matrix to calculate. Defaults to HC1. Currently supports homoskedastic SEs, HC0 & HC1
#' @export
vcov.tsls <- function(X, Z, e.hat,
                     type = "hc1"){
    check_matrix_input(y, X, Z)
    check_individual_matrix_input(e.hat, "residuals e.hat")
    # Type must be one of "homo", "hc0", "hc1", throw error otherwise
    if (!(type %in% c("homo", "hc0", "hc1"))) stop("Type must be one of homo, hc0, hc1")

    n.obs = nrow(y)
    var.z = crossprod(Z)
    cov.xz = crossprod(X, Z)

    if (type == "homo"){
        sigma.2 <- sum(e.hat^2)/n.obs
        v.hat = Matrix::solve(cov.xz, Matrix::solve(var.z, t(cov.xz))) * sigma.2
    } else{

    }
    bread.matrix <- bread.tsls(X, Z, var.z, cov.xz)
    bread.inv <- Matrix::solve(bread.matrix)
    meat.matrix <- meat.tsls(X, Z, e.hat, var.z, cov.xz, type)
    result <- weightedcrossprod(bread.inv, W=meat.matrix)
    # Assign rownames and colnames to correspond to X
    rownames(result) <- colnames(X)
    colnames(result) <- colnames(X)
    return(result)
}

#' Calculate k-class asymptotic var-cov matrix
#'
#' Calculates the var-cov matrix for the k-class estimator given asymptotic theory. This estimator makes use of the special k-class structure. The LIML variance-covariance matrix can asymptotically also be derived as the same as of the 2SLS estimator, but it is better to make use of the special structure of the LIML estimator.
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param n.endogenous The number of endogenous variables in X.
#' @param k A constant to be used for the k-class estimator.
#' @param type The type of var-cov matrix to calculate. Defaults to HC1. Currently supports HC0 & HC1
#' @export
vcov.kclass <- function(X, Z, e.hat, n.endogenous, k, type = "HC1"){
    elimination.matrix.z <- compute_elimination_matrix(Z)
    x.tilde <- rbind(X[, (n.endogenous+1):ncol(X)], (Matrix::Diagonal(nrow(Z))-k*elimination.matrix.z)%*%X[, 1:n.endogenous])
    bread.matrix <- bread.kclassliml(X, Z, n.endogenous, k, elimination.matrix.z, x.tilde)
    bread.inv <- Matrix::solve(bread.matrix)
    meat.matrix <- meat.kclassliml(X, Z, e.hat, n.endogenous, k, elimination.matrix.z, x.tilde, type)
    result <- weightedcrossprod(bread.inv, W=meat.matrix)
    rownames(result) <- colnames(X)
    colnames(result) <- colnames(X)
    return(result)
}

#' Calculate JIVE asymptotic var-cov matrix
#'
#' Calculates the JIVE var-cov matrix given asymptotic theory
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @export
vcov.jive <- function(X, Z, e.hat){
    check_matrix_input(y, X, Z)
    check_individual_matrix_input(e.hat, "residuals e.hat")

    projection.matrix.z <- Z %*% Matrix::solve(crossprod(Z), t(Z))

    bread.matrix <- bread.jive(X, Z, projection.matrix.z)
    bread.inv <- Matrix::solve(bread.matrix)
    meat.matrix <- meat.jive(X, Z, e.hat, projection.matrix.z)
    result <- weightedcrossprod(bread.inv, W=meat.matrix)
    rownames(result) <- colnames(X)
    colnames(result) <- colnames(X)
    return(result)
}

#' Calculate HFUL asymptotic var-cov matrix
#'
#' Calculates the HFUL var-cov matrix given asymptotic theory. The HFUL variance-covariance matrix is derived in \insertCite{hausman2012hful}{RobustIV}.
#' @importFrom Rdpack reprompt
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param alpha.hat The alpha hat computed according to the HFUL procedure.
#' @param projection.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @param SSR Optional. The sum of squared residuals from the structural equation.
#' @references
#'   \insertRef{hausman2012hful}{RobustIV}
#' @export
vcov.hful <- function(X, Z, e.hat,
                      alpha.hat,
                      projection.matrix.z = NULL,
                      SSR = NULL){
    check_matrix_input(y, X, Z)
    check_individual_matrix_input(e.hat, "residuals e.hat")
    check_individual_matrix_input(projection.matrix.z, "projection.matrix.z")

    if (is.null(projection.matrix.z)){
        projection.matrix.z <- compute_projection_matrix(Z)
    }
    if (is.null(SSR)){
        SSR <- sum(e.hat^2)
    }

    bread.matrix <- bread.hful(X, Z, alpha.hat, projection.matrix.z, SSR)
    bread.inv <- Matrix::solve(bread.matrix)
    meat.matrix <- meat.hful(X, Z, e.hat, alpha.hat, projection.matrix.z, SSR)
    result <- weightedcrossprod(bread.inv, W=meat.matrix)
    rownames(result) <- colnames(X)
    colnames(result) <- colnames(X)
    return(result)
}

#' Calculate efficient GMM asymptotic var-cov matrix
#'
#' Calculates the GMM var-cov matrix given asymptotic theory. The efficient GMM var-cov matrix does not take a standard matrix form (only insofar as the bread is the identity matrix). Computes \eqn{(X'ZS^{-1}ZX')^{-1}}
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param cov.zx Optional. \eqn{X'Z} if it is already computed. If it is provided, X has no more effect. If cov.xz and S are provided, Z also has no effect.
#' @param S Optional. The variance-covariance matrix of the moment conditions. If S is provided, e.hat has no effect. If S and cov.xz are provided, X and Z have no effect.
#' @param S.inv Optional. The inverse of S. If S.inv is provided, S has no effect. Otherwise has the same effect as S on making other arguments irrelevant.
#' @export
vcov.gmm <- function(X, Z, e.hat,
                     cov.zx = NULL,
                     S = NULL,
                     S.inv = NULL){
    check_matrix_input(y, X, Z)
    if (is.null(cov.xz)){
        cov.zx <- crossprod(Z, X)
    } else{
        check_individual_matrix_input(cov.zx, "cov.zx")
    }
    if (is.null(S.inv) & is.null(S)){
        check_individual_matrix_input(e.hat, "residuals e.hat")
        S.inv <- compute_efficient_weighting(Z, e.hat)
    } else if (is.null(S.inv) & !is.null(S)){
        check_individual_matrix_input(S, "S")
        S.inv <- Matrix::solve(S)
    } else{
        check_individual_matrix_input(S.inv, "S.inv")
    }

    result <- Matrix::solve(weightedcrossprod(cov.zx, S.inv))
    rownames(result) <- colnames(X)
    colnames(result) <- colnames(X)
    return(result)
}
