#' Estimate 2SLS
#'
#' Estimates the 2SLS estimator.
#' @param formula A formula object of the form y ~ X | Z, where y denotes the dependent variable, X denotes the regressors in the structural equation, and Z denotes the regressors in the first-stage equation. Z must include all exogenous variables, including the ones included in the structural equation.
#' @param data A data frame containing the variables in the formula.
#' @param first.stage A boolean indicating whether to return the matrix of first-stage coefficients. Defaults to FALSE.
#' @param na.action A function to specify the action to be taken if NAs are found in the data.
#' @export
robustiv.tsls <- function(formula, data, first.stage = FALSE, na.action = na.omit){
    if(!is.formula(formula)){
        stop("formula must be a formula object")
    }
    separated.data <- extract_data_from_formula(formula, data, na.action)
    y <- separated.data$y
    X <- separated.data$X
    Z <- separated.data$Z
    rm(separated.data)
    coefficients <- tsls.fit(y, X, Z)
    if (first.stage){
        first.stage.coefficients <- Matrix::solve(crossprod(Z), crossprod(Z,X))
        # Set rownames to be the names of the Z variables
        rownames(first.stage.coefficients) <- colnames(Z)
        # Set colnames to be the names of the X variables
        colnames(first.stage.coefficients) <- colnames(X)
    }
    # Compute fitted values
    fitted.values <- X %*% coefficients
    # Compute residuals
    residuals <- y - fitted.values

    # Compute var-cov-matrix
    model.vcov <- vcov.tsls(X, Z, residuals)

    # Create return list including model frame
    return.list <- list(coefficients = coefficients,
                        fitted.values = fitted.values,
                        residuals = residuals,
                        model.frame = data.frame(y = y, X = X, Z = Z),
                        varcov.matrix = model.vcov)
    if (first.stage){
        return.list$first.stage.coefficients <- first.stage.coefficients
    }
    # Set class
    class(return.list) <- "robustiv.tsls"
    return(return.list)
}

#' Estimate JIVE1
#'
#' Estimates the JIVE1 estimator.
#' @param formula A formula object of the form y ~ X | Z, where y denotes the dependent variable, X denotes the regressors in the structural equation, and Z denotes the regressors in the first-stage equation. Z must include all exogenous variables, including the ones included in the structural equation.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function to specify the action to be taken if NAs are found in the data.
#' @export
robustiv.jive1 <- function(formula, data, na.action = na.omit){
    if(!is.formula(formula)){
        stop("formula must be a formula object")
    }
    separated.data <- extract_data_from_formula(formula, data, na.action)
    y <- separated.data$y
    X <- separated.data$X
    Z <- separated.data$Z
    rm(separated.data)
    coefficients <- jive1.fit(y, X, Z)
    # Compute fitted values
    fitted.values <- X %*% coefficients
    # Compute residuals
    residuals <- y - fitted.values

    # Compute var-cov-matrix
    model.vcov <- vcov.jive(X, Z, residuals)

    # Create return list including model frame
    return.list <- list(coefficients = coefficients,
                        fitted.values = fitted.values,
                        residuals = residuals,
                        model.frame = data.frame(y = y, X = X, Z = Z),
                        varcov.matrix = model.vcov)
    # Set class
    class(return.list) <- "robustiv.jive1"
    return(return.list)
}

#' Estimate JIVE2
#'
#' Estimates the JIVE2 estimator.
#' @param formula A formula object of the form y ~ X | Z, where y denotes the dependent variable, X denotes the regressors in the structural equation, and Z denotes the regressors in the first-stage equation. Z must include all exogenous variables, including the ones included in the structural equation.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function to specify the action to be taken if NAs are found in the data.
#' @export
robustiv.jive2 <- function(formula, data, na.action = na.omit){
    if(!is.formula(formula)){
        stop("formula must be a formula object")
    }
    separated.data <- extract_data_from_formula(formula, data, na.action)
    y <- separated.data$y
    X <- separated.data$X
    Z <- separated.data$Z
    rm(separated.data)
    coefficients <- jive2.fit(y, X, Z)
    # Compute fitted values
    fitted.values <- X %*% coefficients
    # Compute residuals
    residuals <- y - fitted.values

    # Compute var-cov-matrix
    model.vcov <- vcov.jive(X, Z, residuals)

    # Create return list including model frame
    return.list <- list(coefficients = coefficients,
                        fitted.values = fitted.values,
                        residuals = residuals,
                        model.frame = data.frame(y = y, X = X, Z = Z),
                        varcov.matrix = model.vcov)
    # Set class
    class(return.list) <- "robustiv.jive2"
    return(return.list)
}

#' Estimate LIML
#'
#' Estimates the LIML estimator.
#' @param formula A formula object of the form y ~ X | Z, where y denotes the dependent variable, X denotes the regressors in the structural equation, and Z denotes the regressors in the first-stage equation. Z must include all exogenous variables, including the ones included in the structural equation.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function to specify the action to be taken if NAs are found in the data.
#' @export
robustiv.liml <- function(formula, data, na.action = na.omit){
    if(!is.formula(formula)){
        stop("formula must be a formula object")
    }
    separated.data <- extract_data_from_formula(formula, data, na.action)
    y <- separated.data$y
    X <- separated.data$X
    Z <- separated.data$Z
    rm(separated.data)

    # n.endogenous are the number of variables that are in X but not in Z
    # Compute number of instruments in the regressors
    n.endogenous <- ncol(X)-sum(as.numeric(colnames(instruments) %in% colnames(regressors)))
    # Compute elimination matrix for z
    elimination.matrix.z <- compute_elimination_matrix(Z)
    # Compute k associated with LIML for k-class estimators
    liml.k <- compute_liml_k(y, X, Z, n.endogenous, elimination.matrix.z)

    coefficients <- kclass.fit(y, X, Z, liml.k, elimination.matrix.z)
    # Compute fitted values
    fitted.values <- X %*% coefficients
    # Compute residuals
    residuals <- y - fitted.values

    # Compute var-cov-matrix
    model.vcov <- vcov.kclass(X, Z, residuals, n.endogenous, liml.k)

    # Create return list including model frame
    return.list <- list(coefficients = coefficients,
                        fitted.values = fitted.values,
                        residuals = residuals,
                        model.frame = data.frame(y = y, X = X, Z = Z),
                        varcov.matrix = model.vcov,
                        k = liml.k)
    # Set class
    class(return.list) <- "robustiv.liml"
    return(return.list)
}

#' Estimate RJIVE
#'
#' Estimates the RJIVE estimator. RJIVE combines a Ridge penalty in the first-stage regression with the JIVE1 estimator.
#' @param formula A formula object of the form y ~ X | Z, where y denotes the dependent variable, X denotes the regressors in the structural equation, and Z denotes the regressors in the first-stage equation. Z must include all exogenous variables, including the ones included in the structural equation.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function to specify the action to be taken if NAs are found in the data.
#' @export
robustiv.rjive <- function(formula, data, na.action = na.omit){
    if(!is.formula(formula)){
        stop("formula must be a formula object")
    }
    separated.data <- extract_data_from_formula(formula, data, na.action)
    y <- separated.data$y
    X <- separated.data$X
    Z <- separated.data$Z
    rm(separated.data)

    n.endogenous <- ncol(X)-sum(as.numeric(colnames(instruments) %in% colnames(regressors)))

    rjive.penalty <- compute_rjive_penalty(X, n.endogenous, ncol(Z))

    coefficients <- rjive.fit(y, X, Z, rjive.penalty, n.endogenous)
    # Compute fitted values
    fitted.values <- X %*% coefficients
    # Compute residuals
    residuals <- y - fitted.values

    # Compute var-cov-matrix
    model.vcov <- vcov.rjive(X, Z, residuals)

    # Create return list including model frame
    return.list <- list(coefficients = coefficients,
                        fitted.values = fitted.values,
                        residuals = residuals,
                        model.frame = data.frame(y = y, X = X, Z = Z),
                        varcov.matrix = model.vcov,
                        penalty = rjive.penalty)
    # Set class
    class(return.list) <- "robustiv.rjive"
    return(return.list)
}

#' Estimate HFUL
#'
#' Estimates the HFUL estimator. HFUL combines the jackknife estimator structure as for JIVE estimators with the LIML estimator.
#' @param formula A formula object of the form y ~ X | Z, where y denotes the dependent variable, X denotes the regressors in the structural equation, and Z denotes the regressors in the first-stage equation. Z must include all exogenous variables, including the ones included in the structural equation.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function to specify the action to be taken if NAs are found in the data.
#' @param C A constant that is used to compute the jackknife weights. The default value is 1 as suggested by Hausman et al (2012).
#' @export
robustiv.hful <- function(formula, data, na.action = na.omit, C = 1){
    if(!is.formula(formula)){
        stop("formula must be a formula object")
    }
    separated.data <- extract_data_from_formula(formula, data, na.action)
    y <- separated.data$y
    X <- separated.data$X
    Z <- separated.data$Z
    rm(separated.data)

    projection.matrix.z <- compute_projection_matrix(Z)
    leverages <- Matrix::Diagonal(diag(projection.matrix.z))
    x.bar = cbind(y, X)
    weighted.x.bar <-  weightedcrossprod(x.bar, W=leverages, diagonal=TRUE)
    x.bar.hat = crossprod(projection.matrix.z, x.bar)

    alpha.hat <- compute_alpha_hat(y, X, Z, C, projection_matrix.z, weighted.x.bar, x.bar.hat)

    coefficients <- hful.fit(y, X, Z, alpha.hat, projection.matrix.z, leverages, weighted.x.bar, x.bar.hat)
    # Compute fitted values
    fitted.values <- X %*% coefficients
    # Compute residuals
    residuals <- y - fitted.values

    # Compute var-cov-matrix
    model.vcov <- vcov.hful(X, Z, residuals, alpha.hat, projection.matrix.z)

    # Create return list including model frame
    return.list <- list(coefficients = coefficients,
                        fitted.values = fitted.values,
                        residuals = residuals,
                        model.frame = data.frame(y = y, X = X, Z = Z),
                        varcov.matrix = model.vcov,
                        alpha.hat = alpha.hat)
    # Set class
    class(return.list) <- "robustiv.hful"
    return(return.list)
}

#' Estimate GMM
#'
#' Estimates the IV model via GMM. Currently only supports the two-step GMM estimator.
#' @param formula A formula object of the form y ~ X | Z, where y denotes the dependent variable, X denotes the regressors in the structural equation, and Z denotes the regressors in the first-stage equation. Z must include all exogenous variables, including the ones included in the structural equation.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function to specify the action to be taken if NAs are found in the data.
#' @export
robustiv.gmm <- function(formula, data, na.action = na.omit){
    if(!is.formula(formula)){
        stop("formula must be a formula object")
    }
    separated.data <- extract_data_from_formula(formula, data, na.action)
    y <- separated.data$y
    X <- separated.data$X
    Z <- separated.data$Z
    rm(separated.data)

    cov.zx <- crossprod(Z,X)
    # First step
    # Compute first-step coefficients
    first.step <- gmm.fit(y, X, Z, cov.zx = cov.zx)
    first.step.residuals <- y - X %*% first.step
    # Compute efficient weighting matrix
    weighting.matrix <- compute_weighting_matrix(Z, first.step.residuals)

    # Second step
    # Compute second-step coefficients
    coefficients <- gmm.fit(y, X, Z, W=weighting.matrix, cov.zx = cov.zx)

    # Compute fitted values
    fitted.values <- X %*% coefficients
    # Compute residuals
    residuals <- y - fitted.values

    # Compute var-cov-matrix
    model.vcov <- vcov.gmm(X, Z, residuals, cov.zx=cov.zx, W=weighting.matrix)

    # Create return list including model frame
    return.list <- list(coefficients = coefficients,
                        fitted.values = fitted.values,
                        residuals = residuals,
                        model.frame = data.frame(y = y, X = X, Z = Z),
                        varcov.matrix = model.vcov,
                        weighting.matrix = weighting.matrix)
    # Set class
    class(return.list) <- "robustiv.gmm"
    return(return.list)
}
