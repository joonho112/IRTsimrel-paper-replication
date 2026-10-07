# Save, restore, and identify the random-number generator and named streams.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

DEFAULT_RNG_KIND <- c(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

rng_kind_triple <- function() {
  out <- RNGkind()
  if (length(out) != 3L) {
    stop("this R runtime did not return the full RNGkind triple", call. = FALSE)
  }
  stats::setNames(as.character(out), names(DEFAULT_RNG_KIND))
}

restore_rng_kind <- function(kind) {
  if (length(kind) != 3L || anyNA(kind) || any(!nzchar(kind))) {
    stop("an RNG kind state must contain kind, normal.kind, and sample.kind",
      call. = FALSE
    )
  }
  kind <- as.character(kind)
  suppressWarnings(RNGkind(
    kind = kind[1], normal.kind = kind[2],
    sample.kind = kind[3]
  ))
  got <- rng_kind_triple()
  if (!identical(unname(got), unname(kind))) {
    stop("the full RNGkind triple could not be restored", call. = FALSE)
  }
  invisible(got)
}

capture_rng_state <- function() {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(
    kind = rng_kind_triple(), has_seed = has_seed,
    seed = if (has_seed) {
      get(".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    } else {
      NULL
    }
  )
}

restore_rng_state <- function(state) {
  if (!is.list(state) || !all(c("kind", "has_seed", "seed") %in% names(state))) {
    stop("invalid captured RNG state", call. = FALSE)
  }
  restore_rng_kind(state$kind)
  if (isTRUE(state$has_seed)) {
    assign(".Random.seed", state$seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(state$kind)
}

pin_rng <- function(kind = DEFAULT_RNG_KIND[["kind"]],
                    normal.kind = DEFAULT_RNG_KIND[["normal.kind"]],
                    sample.kind = DEFAULT_RNG_KIND[["sample.kind"]]) {
  previous <- rng_kind_triple()
  succeeded <- FALSE
  on.exit(if (!succeeded) restore_rng_kind(previous), add = TRUE)
  suppressWarnings(RNGkind(
    kind = kind, normal.kind = normal.kind,
    sample.kind = sample.kind
  ))
  got <- rng_kind_triple()
  expected <- c(
    kind = kind, normal.kind = normal.kind,
    sample.kind = sample.kind
  )
  if (!identical(unname(got), unname(expected))) {
    stop("the complete RNGkind triple could not be pinned; expected ",
      paste(expected, collapse = " / "), ", got ",
      paste(got, collapse = " / "),
      call. = FALSE
    )
  }
  succeeded <- TRUE
  attr(got, "previous_rng_kind") <- previous
  invisible(got)
}

restore_pinned_rng <- function(pinned) {
  previous <- attr(pinned, "previous_rng_kind", exact = TRUE)
  if (is.null(previous)) {
    stop("pinned RNG record does not contain previous_rng_kind",
      call. = FALSE
    )
  }
  restore_rng_kind(previous)
}

rng_provenance_string <- function(kind = rng_kind_triple()) {
  if (length(kind) != 3L || anyNA(kind)) {
    stop("RNG provenance requires the complete RNGkind triple",
      call. = FALSE
    )
  }
  paste0(
    "rng_kind=", kind[1], "; rng_normal_kind=", kind[2],
    "; rng_sample_kind=", kind[3]
  )
}

with_pinned_rng <- function(code,
                            kind = DEFAULT_RNG_KIND[["kind"]],
                            normal.kind = DEFAULT_RNG_KIND[["normal.kind"]],
                            sample.kind = DEFAULT_RNG_KIND[["sample.kind"]]) {
  previous <- capture_rng_state()
  on.exit(restore_rng_state(previous), add = TRUE)
  pin_rng(kind, normal.kind, sample.kind)
  force(code)
}

.rng_encode_text <- function(x) {
  x <- enc2utf8(as.character(x))
  paste0(nchar(x, type = "bytes"), ":", x)
}

.rng_encode_atomic <- function(x) {
  cls <- paste(class(x), collapse = ",")
  type <- typeof(x)
  values <- vapply(seq_along(x), function(i) {
    z <- x[i]
    value <- if (is.na(z)) {
      if (is.numeric(z) && is.nan(z)) "NaN" else "NA"
    } else if (is.double(z)) {
      if (is.infinite(z)) {
        if (z > 0) "Inf" else "-Inf"
      } else {
        format(z, digits = 17L, scientific = TRUE, trim = TRUE)
      }
    } else if (is.raw(z)) {
      sprintf("%02x", as.integer(z))
    } else {
      as.character(z)
    }
    .rng_encode_text(value)
  }, character(1))
  paste0(
    "atomic<", .rng_encode_text(type), ";", .rng_encode_text(cls), ";",
    length(x), ">[", paste(values, collapse = ";"), "]"
  )
}

rng_canonical_value <- function(x) {
  if (is.null(x)) {
    return("null")
  }
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) nm <- rep("", length(x))
    if (any(nzchar(nm)) && any(!nzchar(nm))) {
      stop("canonical RNG lists must be entirely named or entirely unnamed",
        call. = FALSE
      )
    }
    if (any(nzchar(nm))) {
      if (anyDuplicated(nm)) {
        stop("canonical RNG list names must be unique", call. = FALSE)
      }
      ord <- order(enc2utf8(nm), method = "radix")
      x <- x[ord]
      nm <- nm[ord]
    } else {
      nm <- as.character(seq_along(x))
    }
    pieces <- vapply(
      seq_along(x), function(i) {
        paste0(.rng_encode_text(nm[i]), "=", rng_canonical_value(x[[i]]))
      },
      character(1)
    )
    return(paste0(
      "list<", length(x), ">[", paste(pieces, collapse = "|"),
      "]"
    ))
  }
  if (!is.atomic(x)) {
    stop("RNG identity fields must be atomic values or lists", call. = FALSE)
  }
  .rng_encode_atomic(x)
}

rng_canonical_fields <- function(fields) {
  if (!is.list(fields) || is.null(names(fields)) ||
    any(!nzchar(names(fields))) || anyDuplicated(names(fields))) {
    stop("RNG stream fields must be a uniquely named list", call. = FALSE)
  }
  rng_canonical_value(fields)
}

.rng_sha256 <- function(text) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("package 'digest' is required for named RNG streams", call. = FALSE)
  }
  digest::digest(charToRaw(enc2utf8(text)), algo = "sha256", serialize = FALSE)
}

.rng_hash_to_seed <- function(hash, seed_max) {
  digits <- strsplit(tolower(hash), "", fixed = TRUE)[[1]]
  values <- match(digits, c(as.character(0:9), letters[1:6])) - 1L
  if (anyNA(values)) stop("invalid hexadecimal stream hash", call. = FALSE)
  acc <- 0
  for (z in values) acc <- (acc * 16 + z) %% seed_max
  as.integer(acc + 1)
}

rng_stream_registry <- function(
  streams, namespace,
  key_columns = c("condition_id", "replication", "stream"),
  seed_max = 2147483646L,
  kind = DEFAULT_RNG_KIND[["kind"]],
  normal.kind = DEFAULT_RNG_KIND[["normal.kind"]],
  sample.kind = DEFAULT_RNG_KIND[["sample.kind"]]
) {
  if (!is.data.frame(streams) || !nrow(streams) || anyDuplicated(names(streams))) {
    stop("streams must be a non-empty data frame with unique columns",
      call. = FALSE
    )
  }
  if (!is.character(namespace) || length(namespace) != 1L ||
    is.na(namespace) || !nzchar(namespace)) {
    stop("RNG namespace must be one non-empty string", call. = FALSE)
  }
  missing_keys <- setdiff(key_columns, names(streams))
  if (length(missing_keys)) {
    stop("named streams are missing identity columns: ",
      paste(missing_keys, collapse = ", "),
      call. = FALSE
    )
  }
  for (key in key_columns) {
    z <- streams[[key]]
    if (anyNA(z) || any(!nzchar(trimws(as.character(z))))) {
      stop("named stream identity column has missing/blank values: ", key,
        call. = FALSE
      )
    }
  }
  reserved <- c(
    "rng_namespace", "stream_key", "stream_id", "raw_seed",
    "seed", "collision_probes", "collision_checked", "rng_kind",
    "rng_normal_kind", "rng_sample_kind"
  )
  conflict <- intersect(reserved, names(streams))
  if (length(conflict)) {
    stop("streams use reserved registry columns: ",
      paste(conflict, collapse = ", "),
      call. = FALSE
    )
  }
  seed_max <- as.numeric(seed_max)
  if (length(seed_max) != 1L || !is.finite(seed_max) || seed_max < 1 ||
    seed_max > 2147483646 || seed_max != floor(seed_max)) {
    stop("seed_max must be an integer in [1, 2147483646]", call. = FALSE)
  }
  if (nrow(streams) > seed_max) {
    stop("more named streams than available collision-free integer seeds",
      call. = FALSE
    )
  }

  stream_key <- vapply(seq_len(nrow(streams)), function(i) {
    rng_canonical_fields(c(
      list(namespace = namespace),
      as.list(streams[i, , drop = FALSE])
    ))
  }, character(1))
  if (anyDuplicated(stream_key)) {
    stop("duplicate named RNG stream identity", call. = FALSE)
  }
  stream_id <- vapply(stream_key, .rng_sha256, character(1))
  if (anyDuplicated(stream_id)) {
    stop("SHA-256 collision in named RNG stream identities", call. = FALSE)
  }
  raw_seed <- vapply(stream_id, .rng_hash_to_seed, integer(1),
    seed_max = seed_max
  )

  ord <- order(stream_id, stream_key, method = "radix")
  used <- new.env(hash = TRUE, parent = emptyenv())
  seed <- integer(nrow(streams))
  probes <- integer(nrow(streams))
  for (i in ord) {
    candidate <- raw_seed[i]
    probe <- 0L
    while (exists(as.character(candidate), envir = used, inherits = FALSE)) {
      candidate <- if (candidate >= seed_max) 1L else candidate + 1L
      probe <- probe + 1L
      if (probe >= seed_max) {
        stop("named RNG collision resolver exhausted the seed space",
          call. = FALSE
        )
      }
    }
    assign(as.character(candidate), TRUE, envir = used)
    seed[i] <- as.integer(candidate)
    probes[i] <- probe
  }

  out <- streams
  out$rng_namespace <- namespace
  out$stream_key <- stream_key
  out$stream_id <- stream_id
  out$raw_seed <- raw_seed
  out$seed <- seed
  out$collision_probes <- probes
  out$collision_checked <- TRUE
  out$rng_kind <- kind
  out$rng_normal_kind <- normal.kind
  out$rng_sample_kind <- sample.kind
  assert_rng_stream_registry(out,
    expected_rows = nrow(streams),
    seed_max = seed_max
  )
  out
}

assert_rng_stream_registry <- function(registry, expected_rows = NULL,
                                       seed_max = 2147483646L) {
  required <- c(
    "stream_key", "stream_id", "seed", "collision_probes",
    "collision_checked", "rng_kind", "rng_normal_kind",
    "rng_sample_kind"
  )
  if (!is.data.frame(registry) || !all(required %in% names(registry))) {
    stop("invalid named RNG stream registry schema", call. = FALSE)
  }
  if (!is.null(expected_rows) && nrow(registry) != as.integer(expected_rows)) {
    stop("named RNG stream registry row count mismatch", call. = FALSE)
  }
  if (anyNA(registry[required]) || anyDuplicated(registry$stream_key) ||
    anyDuplicated(registry$stream_id) || anyDuplicated(registry$seed) ||
    any(registry$seed < 1 | registry$seed > seed_max) ||
    any(registry$collision_probes < 0) ||
    !all(registry$collision_checked)) {
    stop("named RNG stream registry failed uniqueness/completeness checks",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

with_named_rng <- function(stream, code) {
  if (!is.data.frame(stream) || nrow(stream) != 1L) {
    stop("with_named_rng requires exactly one registry row", call. = FALSE)
  }
  assert_rng_stream_registry(stream, expected_rows = 1L)
  previous <- capture_rng_state()
  on.exit(restore_rng_state(previous), add = TRUE)
  pin_rng(
    as.character(stream$rng_kind),
    as.character(stream$rng_normal_kind),
    as.character(stream$rng_sample_kind)
  )
  set.seed(as.integer(stream$seed))
  force(code)
}
