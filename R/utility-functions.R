check__overall_matrix_input <- function(y, X, Z){
    # Check inputs
    if(!is.matrix(X)) stop("X must be a matrix")
    if(!is.matrix(Z)) stop("Z must be a matrix")
    if(!is.vector(y)) stop("y must be a vector")
    if(!is.numeric(y)) stop("y must be numeric")
    if(!is.numeric(X)) stop("X must be numeric")
    if(!is.numeric(Z)) stop("Z must be numeric")

    # Check dimensions
    if(nrow(X) != nrow(Z)) stop("X and Z must have the same number of rows")
    if(length(y) != nrow(X)) stop("y must have the same number of rows as X")

    # Check for missing values
    if(any(is.na(y))) stop("y cannot contain missing values")
    if(any(is.na(X))) stop("X cannot contain missing values")
    if(any(is.na(Z))) stop("Z cannot contain missing values")
}

check_individual_matrix_input <- function(A, name){
    # Check inputs
    if(!is.matrix(A)) stop(paste(name, "must be a matrix"))
    if(!is.numeric(A)) stop(paste(name, "must be numeric"))

    # Check for missing values
    if(any(is.na(A))) stop(paste(name, "cannot contain missing values"))
}

compute_projection_matrix <- function(A){
    # Check inputs
    check_individual_matrix_input(A, "A")

    # Compute projection matrix
    #projection.matrix <- as(A %*% Matrix::solve(crossprod(A), t(A)), "symmetricMatrix") # currently not possible to extract diagonal entries from symmetric matrix
    projection.matrix <- A %*% Matrix::solve(crossprod(A), t(A))
    return(projection.matrix)
}

compute_elimination_matrix <- function(A){
    # Check inputs
    check_individual_matrix_input(A, "A")

    # Compute elimination matrix
    elimination.matrix <- Matrix::symmetricMatrix(Matrix::Diagonal(nrow(A))-A %*% Matrix::solve(crossprod(A), t(A)))
    return(elimination.matrix)
}

#' Utility function for weighted crossproducts
#'
#' Utility function to compute a weighted crossproduct of the form $A'WB$
#' @param A A matrix
#' @param B Another matrix. Defaults to A if no input is supplied.
#' @param W The weighting matrix. Defaults to the identity matrix if no input is supplied.
#' @param diagonal Whether the input is a diagonal matrix. If it is, the weighted crossproduct is computed using $(A'W^{1/2})(W^{1/2}A)$ instead of $A'WA$.
#' @export
weightedcrossprod <- function(A, B = NULL, W = NULL, diagonal = FALSE){
    # Check inputs
    check_individual_matrix_input(A, "A")
    if (!is.null(B)){
        check_individual_matrix_input(B, "B")
    }
    if (is.null(W)){
        if (is.null(B)){
            return(Matrix::crossprod(A))
        } else{
            return(Matrix::crossprod(A,B))
        }
    } else{
        #check_individual_matrix_input(W, "W")
        if (is.null(B)){
            if (diagonal){
                return(Matrix::crossprod(sqrt(W) %*% A))
            } else{
                return(Matrix::crossprod(A, W %*% A))
            }
        } else{
            return(Matrix::crossprod(A, W %*% B))
        }
    }
}

extract_data_from_formula <- function(formula, data, na.action){
    # Check inputs
    if(!is.formula(formula)) stop("formula must be a formula")
    if(!is.data.frame(data)) stop("data must be a data frame")

    # Apply na.action
    data = na.action(na.action, data)
    # Extract data using Formula
    instruments = model.matrix(formula, data, rhs=2)
    regressors = model.matrix(formula, data, lhs=1)
    instruments.in.regressors = (colnames(instruments) %in% colnames(regressors))
    rearranged.regressors = cbind(regressors[, !instruments.in.regressors], regressors[,instruments.in.regressors])
    y = model.response(formula, data)

    return(list(y = y, X = rearranged.regressors, Z = instruments))
}
