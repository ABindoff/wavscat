## Convert the Kymatio reference dump into test fixtures.
##
## Run data-raw/make-fixtures.py first, then point this script at its output:
##
##   Rscript data-raw/make-fixtures.R <dump-directory>
##
## Writes tests/testthat/fixtures/kymatio-filter-bank.rds and
## tests/testthat/fixtures/kymatio-transforms.rds. Neither Python nor Kymatio is
## needed to run the test suite afterwards.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript data-raw/make-fixtures.R <dump-directory>")
}
dump_dir <- args[[1L]]
out_dir <- file.path("tests", "testthat", "fixtures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

man <- utils::read.delim(
  file.path(dump_dir, "manifest.tsv"), header = FALSE,
  col.names = c("name", "n"), stringsAsFactors = FALSE
)
lens <- stats::setNames(man$n, man$name)

ref <- function(name) {
  readBin(file.path(dump_dir, paste0(name, ".bin")), "double",
          n = lens[[name]], size = 8, endian = "little")
}

## ---------------------------------------------------------------- filters --
fb <- list(
  scalars = list(
    sigmas = ref("scalar_sigmas"), P = ref("scalar_P"),
    Q = ref("scalar_Q"), xi_max = ref("scalar_xi_max"),
    sigma_psi = ref("scalar_sigma_psi"),
    pair_xi = ref("scalar_pair_xi"), pair_sigma = ref("scalar_pair_sigma"),
    maxsub = ref("scalar_maxsub")
  ),
  morlet = list(
    meta = matrix(ref("morlet_meta"), ncol = 3, byrow = TRUE),
    values = lapply(seq_len(length(ref("morlet_meta")) / 3) - 1L,
                    function(i) ref(paste0("morlet_", i)))
  ),
  gauss = list(
    meta = matrix(ref("gauss_meta"), ncol = 2, byrow = TRUE),
    values = lapply(seq_len(length(ref("gauss_meta")) / 2) - 1L,
                    function(i) ref(paste0("gauss_", i)))
  ),
  tsupp = list(
    meta = matrix(ref("tsupp_meta"), ncol = 2, byrow = TRUE),
    values = ref("tsupp")
  ),
  generator = list(
    meta = matrix(ref("gen_meta"), ncol = 2, byrow = TRUE),
    xi = lapply(seq_len(length(ref("gen_meta")) / 2) - 1L,
                function(i) ref(paste0("gen_xi_", i))),
    sigma = lapply(seq_len(length(ref("gen_meta")) / 2) - 1L,
                   function(i) ref(paste0("gen_sigma_", i)))
  )
)

fac_meta <- matrix(ref("fac_meta"), ncol = 5, byrow = TRUE)
fb$factory <- lapply(seq_len(nrow(fac_meta)), function(k) {
  i <- k - 1L
  phi_nlev <- as.integer(ref(sprintf("fac%d_phi_nlev", i)))
  entry <- list(
    cfg = fac_meta[k, ],
    phi = list(
      j = ref(sprintf("fac%d_phi_j", i)),
      sigma = ref(sprintf("fac%d_phi_sigma", i)),
      levels = lapply(seq_len(phi_nlev) - 1L,
                      function(lev) ref(sprintf("fac%d_phi_lev%d", i, lev)))
    )
  )
  for (ord in 1:2) {
    nlev <- as.integer(ref(sprintf("fac%d_p%d_nlev", i, ord)))
    entry[[paste0("psi", ord)]] <- list(
      xi = ref(sprintf("fac%d_p%d_xi", i, ord)),
      sigma = ref(sprintf("fac%d_p%d_sigma", i, ord)),
      j = ref(sprintf("fac%d_p%d_j", i, ord)),
      levels = lapply(seq_along(nlev), function(n) {
        lapply(seq_len(nlev[n]) - 1L,
               function(lev) ref(sprintf("fac%d_p%d_%d_lev%d", i, ord, n - 1L, lev)))
      })
    )
  }
  entry
})

saveRDS(fb, file.path(out_dir, "kymatio-filter-bank.rds"), compress = "xz")

## ------------------------------------------------------------- transforms --
n_cases <- as.integer(ref("n_cases"))
transforms <- lapply(seq_len(n_cases) - 1L, function(i) {
  cfg <- ref(sprintf("case%d_cfg", i))
  T_raw <- ref(sprintf("case%d_T", i))
  mo <- as.integer(cfg[5])
  out_type <- if (ref(sprintf("case%d_out_type", i)) == 0) "array" else "list"

  entry <- list(
    n = as.integer(cfg[1]), J = as.integer(cfg[2]),
    Q = as.integer(cfg[3:4]), max_order = mo,
    stride = if (cfg[6] < 0) NULL else cfg[6],
    T = if (T_raw == -1) NULL else if (T_raw == -2) "global" else T_raw,
    out_type = out_type,
    x = ref(sprintf("case%d_x", i)),
    meta = list(
      order = ref(sprintf("case%d_meta_order", i)),
      n = matrix(ref(sprintf("case%d_meta_n", i)), ncol = mo, byrow = TRUE),
      j = matrix(ref(sprintf("case%d_meta_j", i)), ncol = mo, byrow = TRUE),
      xi = matrix(ref(sprintf("case%d_meta_xi", i)), ncol = mo, byrow = TRUE)
    )
  )
  if (out_type == "array") {
    d <- as.integer(ref(sprintf("case%d_dim", i)))
    entry$out <- matrix(ref(sprintf("case%d_out", i)), nrow = d[1], byrow = TRUE)
  } else {
    npaths <- as.integer(ref(sprintf("case%d_npaths", i)))
    entry$out <- lapply(seq_len(npaths) - 1L,
                        function(k) ref(sprintf("case%d_p%d", i, k)))
  }
  entry
})

saveRDS(transforms, file.path(out_dir, "kymatio-transforms.rds"), compress = "xz")

cat("wrote fixtures to", out_dir, "\n")
for (f in list.files(out_dir, full.names = TRUE)) {
  cat(sprintf("  %-32s %8.1f KB\n", basename(f), file.size(f) / 1024))
}
