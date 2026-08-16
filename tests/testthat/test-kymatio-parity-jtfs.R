# Numerical parity with Kymatio for joint time-frequency scattering.
#
# The fixture holds Kymatio's dict output, keyed by path, rather than its
# stacked array. As of 0.4.0.dev0 Kymatio stacks format = "time" coefficients in
# generator order while reporting metadata in sorted order, so its array and its
# meta() disagree with each other. wavscat sorts both, so comparing by key is
# the only meaningful test.

test_that("joint time-frequency transforms match Kymatio", {
  cases <- readRDS(test_path("fixtures", "kymatio-jtfs.rds"))

  for (case in cases) {
    tag <- sprintf("n=%d J=%d J_fr=%d Q=%d/%d Q_fr=%d T=%s F=%s format=%s",
                   case$n, case$J, case$J_fr, case$Q[1], case$Q[2], case$Q_fr,
                   if (is.null(case$T)) "NULL" else as.character(case$T),
                   if (is.null(case$F)) "NULL" else as.character(case$F),
                   case$format)

    op <- scattering_jtfs(
      n = case$n, J = case$J, J_fr = case$J_fr, Q = case$Q, Q_fr = case$Q_fr,
      T = case$T, F = case$F, format = case$format, out_type = "list"
    )
    S <- scat_transform(op, case$x)
    meta <- S$meta

    expect_equal(nrow(meta), length(case$keys), label = paste("path count,", tag))

    # Rebuild Kymatio's key from our metadata: zero-based indices, with the
    # components that apply to this format and order.
    mine <- vapply(seq_len(nrow(meta)), function(r) {
      if (meta$order[r] == 0L) {
        return("")
      }
      paste(c(
        if (case$format == "time") meta$n1[r] - 1L,
        if (meta$order[r] == 2L) meta$n2[r] - 1L,
        meta$n_fr[r] - 1L
      ), collapse = ",")
    }, "")

    expect_setequal(mine, case$keys)

    idx <- match(mine, case$keys)
    scale <- max(abs(unlist(case$coefs, use.names = FALSE)))
    err <- vapply(seq_len(nrow(meta)), function(r) {
      want <- case$coefs[[idx[r]]]
      got <- S$coef[[r]]
      # Drop the trailing single-channel axis.
      got <- array(got, dim = utils::head(dim(got), -1L))
      if (length(got) != length(want)) {
        return(Inf)
      }
      max(abs(as.numeric(got) - as.numeric(want)))
    }, 0)

    expect_lt(max(err) / scale, 1e-10, label = paste("coefficients,", tag))
  }
})
