# robustIV

## Description

This package intends to implement procedures related to instrumental variables estimation that is robust to the weak & many instruments problem. In its current form, the package entails multiple estimators that exhibit this robustness. In the future, I plan to add robust test procedures, support for the `sandwich` package, and integration with a table generator.

This package is still in alpha development and the estimators have not been rigorously tested yet. *Use at your own risk.* It would be better not to use this package for research papers just yet.

As this package is still in its alpha phase, variable names, function names, and function inputs can still change without ensuring backward-compatibility.

This is the first time I am creating a package. I am grateful for any input.

## Currently supported estimators
The package currently provides functions for the following estimators It calculates the coefficients and a heteroskedasticity-consistent variance-covariance matrix. Currently, no estimators permit any autocorrelation in the error term (i.e. also not clustered observations). 

The supported estimators are as follows. Abbreviations used in functions are provided in brackets:

- Two-stage-least-squares (tsls): Standard 2SLS estimation. Currently has no routines for dealing with weak and/or many instruments.
- Limited Information Maximum Likelihood (liml & kclassliml): LIML and k-class LIML estimation. LIML has similiar asymptotic properties, but usually better finite-sample properties than 2SLS, especially when there are many instruments or weak instruments.
- Jackknife instrumental variables estimators (jive1 & jive2): These estimators implement a jackknife procedure (leave-one-out) in the first stage of 2SLS, thus providing a robust estimator that is approximately unbiased.
- RJIVE (rjive): This estimator combines the JIVE1 estimator with a ridge regression in the first stage. This procedure generally performs very well even with few instruments, but in particular it permits using more instruments than observations. Currently does not have a routine for estimating the variance-covariance matrix.
- HFUL (hful): This estimator combines LIML with the first-stage jackknife procedure. It has attractive properties in terms of robustness to weak & many instruments as well as efficency.
- GMM (gmm): Standard twostep efficient GMM estimation. Currently has no routines for dealing with weak and/or many instruments.

The estimators can generally be called using `robustiv.estimatorname`, with `estimatorname` being replaced by the name of the estimator. For example, to use LIML estimator, one would write `robustiv.liml`.

## Changelog

- v0.1.0: Initial alpha release.

## To-do list

### High priority

- Add Wrapper function for included estimators
- Add support for summary function
- Add var-cov matrix estimation for RJIVE
- Include references in the documentation

### Medium priority
- Add support for sandwich package
- Add support for clustered observations
- Include procedures for 2SLS that are robust to weak & many instruments

### Low priority
- Extend documentation
