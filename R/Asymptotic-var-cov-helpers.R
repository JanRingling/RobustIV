#' Calculate the bread for 2SLS asymptotic var-cov matrix
#'
#' Calculates \eqn{X'Z(Z'Z)^{-1}Z'X}, the bread matrix for 2SLS
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @param var.z Optional input. $n$ times variance-covariance matrix of Z ($Z'Z$)
#' @param cov.xz Optional input. $n$ times the covariance between X and Z ($X'Z$)
#' @export
bread.tsls <- function(X, Z,
                       var.z = NULL,
                       cov.xz = NULL){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    if (is.null(var.z)){
        var.z = crossprod(Z)
    } else{
        check_individual_matrix_input(var.z, "var.z")
    }
    if (is.null(cov.xz)){
        cov.xz = crossprod(X, Z)
    } else{
        check_individual_matrix_input(cov.xz, "cov.xz")
    }
    bread.matrix = cov.xz %*% Matrix::solve(var.z, t(cov.xz))
    return(bread.matrix)
}

#' Calculate the meat for 2SLS asymptotic var-cov matrix
#'
#' Calculates \eqn{X'Z(Z'Z)^{-1}Z'X}, the bread matrix for 2SLS
#' @param X An nxk matrix of regressors from the structural equation
#' @param Z An nxp matrix of regressors from the first-stage equation
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param var.z Optional input. $n$ times the variance-covariance matrix of Z ($Z'Z$)
#' @param cov.xz Optional input. $n$ times the covariance between X and Z ($X'Z$)
#' @param type The type of var-cov matrix to calculate. Defaults to HC1. Currently supports HC0 & HC1
#' @export
meat.tsls <- function(X, Z, e.hat,
                       var.z = NULL,
                       cov.xz = NULL,
                      type = "hc1"){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    check_individual_matrix_input(e.hat, "residuals")
    if (is.null(var.z)){
        var.z = crossprod(Z)
    } else{
        check_individual_matrix_input(var.z, "var.z")
    }
    if (is.null(cov.xz)){
        cov.xz = crossprod(X, Z)
    } else{
        check_individual_matrix_input(cov.xz, "cov.xz")
    }
    n.obs = nrow(X)
    n.estimated.structural = ncol(X)
    partial.meat = Matrix::solve(var.z, t(cov.xz))
    if (type=="hc0"){
        inner.diagonal = Matrix::Diagonal(e.hat^2)
    } else if (type=="hc1"){
        inner.diagonal = Matrix::Diagonal(e.hat^2 * n.obs / (n.obs-n.estimated.structural))
    } else if (type=="homo"){
        return(sum(e.hat^2))
    }
    omega = Matrix::SymmetricMatrix(weightedcrossprod(Z, W=inner.diagonal, diagonal=TRUE))

    meat.matrix = weightedcrossprod(partial.meat, W=omega)
    return(bread.matrix)
}

#' Calculate the bread for LIML asymptotic var-cov matrix
#'
#' Calculates \eqn{\tilde{X}'X}, the bread matrix for LIML.
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param n.endogenous The number of endogenous variables in X.
#' @param k A constant to be used for the k-class estimator. Must always be provided but has no effect if x.tilde is provided.
#' @param elimination.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @param x.tilde An optional argument to input x.tilde if it is already computed. If x.tilde is input, elimination.matrix.z and k has no effect. If x.tilde is not input, k is required and the computation is accelerated by inputting elimination.matrix.z
#' @export
bread.kclassliml <- function(X, Z,
                       n.endogenous,
                       k,
                       elimination.matrix.z = NULL,
                       x.tilde = NULL
                       ){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    if (is.null(elimination.matrix.z)){
        elimination.matrix.z = Matrix::Diagonal(ncol(Z)) - Matrix::symmetricMatrix(Z %*% Matrix::solve(crossprod(Z), t(Z)))
    } else{
        check_individual_matrix_input(elimination.matrix.z, "elimination matrix of Z")
    }
    if (is.null(x.tilde)){
        x.tilde = rbind(X[, (n.endogenous+1):ncol(X)], (Matrix::Diagonal(nrow(Z))-k*elimination.matrix.z)%*%X[, 1:n.endogenous])
    } else{
        check_individual_matrix_input(x.tilde, "x.tilde")
    }
    bread.matrix = crossprod(X, x.tilde)
    return(bread.matrix)
}

#' Calculate the meat for LIML asymptotic var-cov matrix
#'
#' Calculates \eqn{\tilde{X}'X}, the bread matrix for LIML.
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param n.endogenous The number of endogenous variables in X.
#' @param k A constant to be used for the k-class estimator
#' @param elimination.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @param x.tilde An optional argument to input x.tilde if it is already computed. If x.tilde is input, elimination.matrix.z and k has no effect.
#' @param type The type of var-cov matrix to calculate. Defaults to HC1. Currently supports HC0 & HC1
#' @export
meat.kclassliml <- function(X, Z,
                      e.hat,
                       n.endogenous,
                       k,
                       elimination.matrix.z = NULL,
                       x.tilde = NULL,
                      type = "hc1"){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    check_individual_matrix_input(e.hat, "residuals")
    if (!(type %in% c("hc0", "hc1"))) stop("type must be hc0 or hc1")
    if (is.null(elimination.matrix.z)){
        elimination.matrix.z = Matrix::Diagonal(ncol(Z)) - Matrix::symmetricMatrix(Z %*% Matrix::solve(crossprod(Z), t(Z)))
    } else{
        check_individual_matrix_input(elimination.matrix.z, "elimination matrix of Z")
    }
    if (is.null(x.tilde)){
        x.tilde = rbind(X[, (n.endogenous+1):ncol(X)], (Matrix::Diagonal(nrow(Z))-k*elimination.matrix.z)%*%X[, 1:n.endogenous])
    } else{
        check_individual_matrix_input(x.tilde, "x.tilde")
    }
    if (type=="hc0"){
        adjustment = 1
    } else if (type=="hc1"){
        adjustment = nrow(Z) / (nrow(Z)-ncol(X))
    }

    meat.matrix = weightedcrossprod(x.tilde, W=adjustment * Matrix::Diagonal(x=e.hat^2), diagonal=TRUE)
    return(meat.matrix)
}

#' Calculate the bread for JIVE asymptotic var-cov matrix
#'
#' Calculates \eqn{\hat H=X'PX-\sum_i X_i P_{ii} X_i'}, the bread matrix for JIVE estimators.
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param projection.matrix.z An optional argument to input the projection matrix for Z if it is already computed.
#' @export
bread.jive <- function(X, Z,
                       projection.matrix.z = NULL
){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    if (is.null(elimination.matrix.z)){
        projection.matrix.z = compute_projection_matrix(Z)
    } else{
        check_individual_matrix_input(projection.matrix.z, "projection matrix of Z")
    }
    bread.matrix = crossprod(X, projection.matrix.z %*% X) - weightedcrossprod(X, W=Matrix::Diagonal(x=diag(projection.matrix.z)))
    return(bread.matrix)
}

#' Calculate the meat for JIVE asymptotic var-cov matrix
#'
#' Calculates \eqn{\hat{\Sigma}'X}, the bread matrix for JIVE estimators.
#' @param X An nxk matrix of regressors from the structural equation.
#' @param Z An nxp matrix of regressors from the first-stage equation.
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param projection.matrix.z An optional argument to input the elimination matrix for Z if it is already computed.
#' @export
meat.jive <- function(X, Z,
                      e.hat,
                      projection.matrix.z = NULL){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    check_individual_matrix_input(e.hat, "residuals")
    if (is.null(projection.matrix.z)){
        projection.matrix.z = compute_projection_matrix(Z)
    } else{
        check_individual_matrix_input(projection.matrix.z, "projection matrix of Z")
    }

    x.tilde = X / (1-diag(projection.matrix.z))
    h.tilde = weightedcrossprod(X, x.tilde, W=projection.matrix.z) - weightedcrossprod(X, x.tilde, W=Matrix::Diagonal(x=diag(projection.matrix.z)))
    x.bar = projection.matrix.z %*% X
    z.tilde = t(solve(crossprod(Z), t(Z)))
    e.hat.squared = e.hat^2

    first.term = weightedcrossprod(x.tilde, W=Matrix::Diagonal(x=e.hat.squared))-
        weightedcrossprod(X, x.bar, W=Matrix::Diagonal(x=e.hat.squared*diag(projection.matrix.z)))-
        weightedcrossprod(x.bar, X, W=Matrix::Diagonal(x=e.hat.squared*diag(projection.matrix.z)))
    second.term <- 0
    for (k in 1:ncol(Z)){
        for (l in 1:ncol(Z)){
            for (i in 1:nrow(X)){
                front.term <- z.tilde[i,k] * z.tilde[i,l] * X[i,] * e.hat[i]
                back.term <- Z[i,k] * Z[i,l] * X[i,] * e.hat[i]
            }
            second.term <- second.term + front.term %*% t(back.term)
        }
    }
    meat.matrix <- first.term-second.term
    return(meat.matrix)
}

#' Calculate the bread for RJIVE asymptotic var-cov matrix (todo)
#'
#' Calculates \eqn{\hat H=X'P^A X-\sum_i X_i P_{ii}^A X_i'}, the bread matrix for RJIVE estimators.
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param projection.matrix.z An optional argument to input the augmented projection matrix $Z(Z'Z+\Lambda'\Lambda)^{-1}Z'$ for Z if it is already computed.
#' @param penalty.matrix.z An optional argument to input the penalty matrix $\Lambda$ for Z if it is already computed.
#' @export
bread.rjive <- function(X, Z,
                       projection.matrix.z = NULL,
                       penalty.matrix.z = NULL
){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    if (is.null(elimination.matrix.z)){
        projection.matrix.z = as(Z %*% Matrix::solve(crossprod(Z)+penalty.matrix.z, Z), "symmetricMatrix")
    } else{
        check_individual_matrix_input(projection.matrix.z, "projection matrix of Z")
    }
    bread.matrix = crossprod(X, projection.matrix.z %*% X) - weightedcrossprod(X, W=Matrix::Diagonal(x=diag(projection.matrix.z)))
    return(bread.matrix)
}

#' Calculate the bread for HFUL asymptotic var-cov matrix
#'
#' Calculates \eqn{\hat H=X'P X-\sum_i X_i P_{ii} X_i'-\hat\alpha X'X}, the bread matrix for HFUL estimators.
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param alpha.hat  \eqn{\hat\alpha} used when computing the estimator
#' @param projection.matrix.z An optional argument to input the projection matrix of Z if it is already computed.
#' @export
bread.hful <- function(X, Z, alpha.hat, projection.matrix.z = NULL
){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    if (is.null(elimination.matrix.z)){
        projection.matrix.z = as(Z %*% Matrix::solve(crossprod(Z)+penalty.matrix.z, Z), "symmetricMatrix")
    } else{
        check_individual_matrix_input(projection.matrix.z, "projection matrix of Z")
    }
    bread.matrix = crossprod(X, projection.matrix.z %*% X) - weightedcrossprod(X, W=Matrix::Diagonal(x=diag(projection.matrix.z)))-alpha.hat*crossprod(X)
    return(bread.matrix)
}

#' Calculate the meat for HFUL asymptotic var-cov matrix
#'
#' Calculates the meat matrix for HFUL estimators.
#' @param X An nxk matrix of regressors from the structural equation. The endogenous variables must be listed first.
#' @param Z An nxp matrix of regressors from the first-stage equation. The included instruments (from the structural equation) must be listed first
#' @param e.hat An nx1 vector of residuals from the structural equation
#' @param alpha.hat  \eqn{\hat\alpha} used when computing the estimator
#' @param projection.matrix.z An optional argument to input the projection matrix of Z if it is already computed.
#' @param SSR Optional. The sum of squared residuals from the structural equation
#' @export
meat.hful <- function(X, Z, e.hat, alpha.hat, projection.matrix.z = NULL, SSR
){
    check_individual_matrix_input(X)
    check_individual_matrix_input(Z)
    if (is.null(elimination.matrix.z)){
        projection.matrix.z = as(Z %*% Matrix::solve(crossprod(Z)+penalty.matrix.z, Z), "symmetricMatrix")
    } else{
        check_individual_matrix_input(projection.matrix.z, "projection matrix of Z")
    }
    if (is.null(SSR)){
        SSR = sum(e.hat^2)
    }
    gamma.hat = crossprod(X, e.hat)/SSR
    x.hat <- X - tcrossprod(e.hat, gamma.hat)
    x.dot <- P %*% x.hat
    z.tilde <- t(Matrix::solve(crossprod(Z), t(Z)))
    first.term <- crossprod(x.dot)-
        weightedcrossprod(x.dot, x.hat, W=Matrix::Diagonal(x=diag(projection.matrix.z)), diagonal=TRUE)-
        weightedcrossprod(x.hat, x.dot, W=Matrix::Diagonal(x=diag(projection.matrix.z)), diagonal=TRUE)
    second.term <- 0
    for (k in 1:ncol(Z)){
        for (l in 1:ncol(Z)){
            for (i in 1:nrow(X)){
                front.term <- z.tilde[i,k] * z.tilde[i,l] * x.hat[i,] * e.hat[i]
                back.term <- Z[i,k] * Z[i,l] * x.hat[i,] * e.hat[i]
            }
            second.term <- second.term + front.term %*% t(back.term)
        }
    }
    meat.matrix <- first.term-second.term
    return(meat.matrix)
}
