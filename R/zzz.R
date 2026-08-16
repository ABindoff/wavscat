# recipes is an optional integration, so its methods are registered at load
# time rather than declared in NAMESPACE. That keeps recipes in Suggests: users
# who never touch tidymodels do not have to install it.

.onLoad <- function(libname, pkgname) {
  s3_register("recipes::prep", "step_scattering")
  s3_register("recipes::bake", "step_scattering")
  s3_register("recipes::tidy", "step_scattering")
  s3_register("recipes::required_pkgs", "step_scattering")
  invisible()
}

#' Register an S3 method for a generic in a package that may not be installed
#'
#' Vendored from rlang, which stopped exporting it at version 1.3.0 and
#' recommends copying it instead. Registration is deferred until the owning
#' package loads, so nothing breaks when it is absent.
#'
#' @param generic Generic as `"pkg::generic"`.
#' @param class Class name.
#' @param method Optional method function; found by name if omitted.
#' @noRd
s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]

  caller <- parent.frame()

  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) asNamespace(environmentName(top)) else caller
  }
  get_method <- function(method) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    } else {
      method
    }
  }

  register <- function(...) {
    envir <- asNamespace(package)
    # Fetched afresh each time, since devtools::load_all() may have replaced it.
    method_fn <- get_method(method)
    stopifnot(is.function(method_fn))
    if (exists(generic, envir)) {
      registerS3method(generic, class, method_fn, envir = envir)
    } else if (identical(Sys.getenv("NOT_CRAN"), "true")) {
      warning(sprintf(
        "Can't find generic `%s` in package %s to register S3 method.",
        generic, package
      ), call. = FALSE)
    }
  }

  # Register on a hook as well, in case the package is unloaded and reloaded.
  setHook(packageEvent(package, "onLoad"), function(...) register())
  if (isNamespaceLoaded(package)) {
    register()
  }
  invisible()
}
