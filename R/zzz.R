.onLoad <- function(libname, pkgname) {
  if (requireNamespace("vctrs", quietly = TRUE) &&
      requireNamespace("pillar", quietly = TRUE)) {
    vctrs::s3_register("pillar::type_sum", "accel", method = function(x) "accel")
  }
}
