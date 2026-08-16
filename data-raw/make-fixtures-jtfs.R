## Convert the Kymatio JTFS reference dump into a test fixture.
##
## Run data-raw/make-fixtures-jtfs.py first, then:
##
##   Rscript data-raw/make-fixtures-jtfs.R <dump-directory>
##
## Writes tests/testthat/fixtures/kymatio-jtfs.rds.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript data-raw/make-fixtures-jtfs.R <dump-directory>")
}
dump_dir <- args[[1L]]
out_dir <- file.path("tests", "testthat", "fixtures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

man <- utils::read.delim(
  file.path(dump_dir, "manifest.tsv"), header = FALSE,
  col.names = c("name", "n"), stringsAsFactors = FALSE
)
lens <- stats::setNames(man$n, man$name)
ref <- function(nm) {
  readBin(file.path(dump_dir, paste0(nm, ".bin")), "double",
          n = lens[[nm]], size = 8, endian = "little")
}

n_cases <- as.integer(ref("n_cases"))
cases <- lapply(seq_len(n_cases) - 1L, function(i) {
  cfg <- ref(sprintf("case%d_cfg", i))
  T_raw <- ref(sprintf("case%d_T", i))
  F_raw <- ref(sprintf("case%d_F", i))
  npaths <- as.integer(ref(sprintf("case%d_npaths", i)))
  keys <- matrix(ref(sprintf("case%d_keys", i)), ncol = 3L, byrow = TRUE)
  keylen <- as.integer(ref(sprintf("case%d_keylen", i)))
  ndim <- as.integer(ref(sprintf("case%d_ndim", i)))
  shape <- matrix(as.integer(ref(sprintf("case%d_shape", i))),
                  ncol = 2L, byrow = TRUE)

  coefs <- lapply(seq_len(npaths), function(p) {
    v <- ref(sprintf("case%d_c%d", i, p - 1L))
    # Kymatio arrays are row-major; a joint path is [band, time].
    if (ndim[p] == 2L) t(matrix(v, nrow = shape[p, 2L])) else v
  })

  list(
    n = as.integer(cfg[1]), J = as.integer(cfg[2]), J_fr = as.integer(cfg[3]),
    Q = as.integer(cfg[4:5]), Q_fr = as.integer(cfg[6]),
    format = if (cfg[7] == 0) "time" else "joint",
    T = if (T_raw == -1) NULL else if (T_raw == -2) "global" else T_raw,
    F = if (F_raw == -1) NULL else if (F_raw == -2) "global" else F_raw,
    x = ref(sprintf("case%d_x", i)),
    # One key per path, as a comma-separated string of zero-based indices.
    keys = vapply(seq_len(npaths), function(p) {
      paste(keys[p, seq_len(keylen[p])], collapse = ",")
    }, ""),
    coefs = coefs
  )
})

saveRDS(cases, file.path(out_dir, "kymatio-jtfs.rds"), compress = "xz")

cat("wrote fixtures to", out_dir, "\n")
for (f in list.files(out_dir, full.names = TRUE)) {
  cat(sprintf("  %-32s %8.1f KB\n", basename(f), file.size(f) / 1024))
}
