# =====================================================================
# analysis_utils.R
# Utility functions and checkpoint system
# =====================================================================

# -------------------------------
# Checkpoint system
# -------------------------------

CHECKPOINT_DIR <- file.path(
  Sys.getenv("PROJECT_ROOT", "/gpfs/data1/vclgp/lmaden/chpt1"),
  "checkpoints"
)
dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)

checkpoint_exists <- function(name) {
  file.exists(file.path(CHECKPOINT_DIR, paste0(name, ".rds")))
}

save_checkpoint <- function(name, data) {
  ts_start <- Sys.time()
  file_path <- file.path(CHECKPOINT_DIR, paste0(name, ".rds"))
  
  # Add metadata
  checkpoint_data <- list(
    data = data,
    timestamp = Sys.time(),
    session_info = sessionInfo()$R.version$version.string
  )
  
  saveRDS(checkpoint_data, file_path, compress = "xz")
  ts_end <- Sys.time()
  elapsed <- as.numeric(difftime(ts_end, ts_start, units = "secs"))
  
  file_size <- file.info(file_path)$size
  size_mb <- round(file_size / 1024^2, 2)
  
  log_progress(sprintf("✓ Checkpoint saved: %s (%.2f MB, %.1f sec)", 
                       name, size_mb, elapsed))
  invisible(file_path)
}

load_checkpoint <- function(name) {
  file_path <- file.path(CHECKPOINT_DIR, paste0(name, ".rds"))
  if (!file.exists(file_path)) {
    stop(sprintf("Checkpoint not found: %s", name))
  }
  
  ts_start <- Sys.time()
  checkpoint_data <- readRDS(file_path)
  ts_end <- Sys.time()
  elapsed <- as.numeric(difftime(ts_end, ts_start, units = "secs"))
  
  file_size <- file.info(file_path)$size
  size_mb <- round(file_size / 1024^2, 2)
  
  log_progress(sprintf("✓ Checkpoint loaded: %s (%.2f MB, %.1f sec)", 
                       name, size_mb, elapsed))
  
  checkpoint_data$data
}

clear_checkpoint <- function(name) {
  file_path <- file.path(CHECKPOINT_DIR, paste0(name, ".rds"))
  if (file.exists(file_path)) {
    unlink(file_path)
    log_progress(sprintf("✗ Checkpoint cleared: %s", name))
  }
}

clear_all_checkpoints <- function() {
  files <- list.files(CHECKPOINT_DIR, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) > 0) {
    unlink(files)
    log_progress(sprintf("✗ Cleared %d checkpoints", length(files)))
  }
}

list_checkpoints <- function() {
  files <- list.files(CHECKPOINT_DIR, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    cat("No checkpoints found.\n")
    return(invisible(NULL))
  }
  
  info <- data.frame(
    name = basename(files),
    size_mb = round(file.info(files)$size / 1024^2, 2),
    modified = file.info(files)$mtime,
    stringsAsFactors = FALSE
  )
  info <- info[order(info$modified), ]
  
  cat("\nAvailable checkpoints:\n")
  cat("─────────────────────────────────────────────────────────────\n")
  print(info, row.names = FALSE)
  cat("─────────────────────────────────────────────────────────────\n")
  invisible(info)
}

# -------------------------------
# Logging functions
# -------------------------------

log_progress <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, msg))
  flush.console()
}

log_section <- function(title) {
  cat("\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat(sprintf("  %s\n", title))
  cat("════════════════════════════════════════════════════════════════\n")
  flush.console()
}

log_subsection <- function(title) {
  cat(sprintf("\n── %s ──\n", title))
  flush.console()
}

# -------------------------------
# Statistical utility functions
# -------------------------------

nm_ad <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  1.4826 * median(abs(x - median(x)), na.rm = TRUE)
}

q95_abs <- function(x) {
  x <- abs(x[is.finite(x)])
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, 0.95, na.rm = TRUE))
}

zscale <- function(x) as.numeric(scale(x))

central_limits <- function(x, probs = c(0.01, 0.99), symmetric_center = NULL) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  if (!is.null(symmetric_center)) {
    q <- stats::quantile(abs(x - symmetric_center), probs = max(probs), na.rm = TRUE)
    c(symmetric_center - q, symmetric_center + q)
  } else {
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
    c(max(q[1], min(x)), min(q[2], max(x)))
  }
}

wq <- function(y, w, probs = c(0.10, 0.50, 0.90)) {
  i <- is.finite(y) & is.finite(w) & (w > 0)
  if (!any(i)) return(rep(NA_real_, length(probs)))
  y <- y[i]; w <- w[i]
  o <- order(y); y <- y[o]; w <- w[o]
  cw <- cumsum(w); cw <- cw / cw[length(cw)]
  d  <- !duplicated(cw)
  approx(x = cw[d], y = y[d], xout = probs, method = "linear", ties = "ordered", rule = 2)$y
}

choose_hex_bins <- function(n, min_bins = 40, max_bins = 120) {
  k <- floor(sqrt(max(n, 1)) / 4)
  max(min_bins, min(max_bins, k))
}

top_levels <- function(f, k = 6L) {
  lv <- sort(table(f), decreasing = TRUE)
  names(lv)[seq_len(min(k, length(lv)))]
}

maybe_sample_df <- function(df) {
  if (!isTRUE(EDA_ENABLE)) return(df)
  n <- nrow(df); if (!n) return(df)
  group_col <- NULL
  for (g in EDA_GROUP_HINTS) if (g %in% names(df)) { group_col <- g; break }
  if (!is.null(group_col)) {
    df <- df %>%
      dplyr::group_by(.data[[group_col]]) %>%
      dplyr::group_modify(~{
        .x <- .x
        n_cap <- min(EDA_MAX_PER_GROUP, nrow(.x))
        if (EDA_SAMPLE_FRAC < 1) .x <- dplyr::slice_sample(.x, prop = EDA_SAMPLE_FRAC, replace = FALSE)
        if (nrow(.x) > n_cap) .x <- dplyr::slice_sample(.x, n = n_cap, replace = FALSE)
        .x
      }) %>%
      dplyr::ungroup()
  } else {
    target <- min(EDA_MAX_N, ceiling(n * EDA_SAMPLE_FRAC))
    if (n > target) df <- dplyr::slice_sample(df, n = target, replace = FALSE)
  }
  if (nrow(df) > EDA_MAX_N) df <- dplyr::slice_sample(df, n = EDA_MAX_N, replace = FALSE)
  df
}

maybe_sample_vec <- function(x) {
  if (!isTRUE(EDA_ENABLE)) return(x)
  n <- length(x)
  target <- min(EDA_MAX_N, ceiling(n * EDA_SAMPLE_FRAC))
  if (n > target) x <- x[sample.int(n, target)]
  x
}

# -------------------------------
# Plotting helper functions
# -------------------------------

ribbon_summary <- function(df, x_col, y_col, w_col,
                           xlim = NULL, ylim = NULL,
                           nbins = 120, min_n_per_bin = 25) {
  df <- maybe_sample_df(df)
  dt <- data.table::as.data.table(df[, c(x_col, y_col, w_col)])
  setnames(dt, c("x","y","w"))
  dt <- dt[is.finite(x) & is.finite(y) & is.finite(w)]
  if (is.null(xlim)) xlim <- central_limits(dt$x)
  if (is.null(ylim)) ylim <- central_limits(dt$y, symmetric_center = 0)
  dt <- dt[x >= xlim[1] & x <= xlim[2] & y >= ylim[1] & y <= ylim[2]]
  breaks <- seq(xlim[1], xlim[2], length.out = nbins + 1L)
  mids   <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
  dt[, bin := cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)]
  rib <- dt[, {
    qs  <- wq(y, w, probs = c(0.05, 0.10, 0.50, 0.90, 0.95))
    list(n = .N, q05 = qs[1], q10 = qs[2], q50 = qs[3], q90 = qs[4], q95 = qs[5])
  }, by = bin][!is.na(bin)]
  rib <- rib[n >= min_n_per_bin]
  rib[, x := mids[bin]]
  list(summary = rib, xlim = xlim, ylim = ylim)
}

fit_bam_smooth_weighted <- function(df, x_col, y_col, w_col,
                                    xlim = NULL, ylim = NULL,
                                    k = 80, nthreads = MGCV_THREADS, ndraw = 600) {
  df <- maybe_sample_df(df)
  dat <- df[, c(x_col, y_col, w_col)]
  names(dat) <- c("x","y","w")
  dat <- dat[is.finite(dat$x) & is.finite(dat$y) & is.finite(dat$w), , drop = FALSE]
  if (is.null(xlim)) xlim <- central_limits(dat$x)
  if (is.null(ylim)) ylim <- central_limits(dat$y, symmetric_center = 0)
  dat <- dat[dat$x >= xlim[1] & dat$x <= xlim[2] & dat$y >= ylim[1] & dat$y <= ylim[2], ]
  if (nrow(dat) < 200) {
    lo <- stats::lowess(dat$x, dat$y, f = 0.6)
    pred <- data.frame(x = lo$x, fit = lo$y, lo = NA_real_, hi = NA_real_)
    return(list(mod = NULL, pred = pred, xlim = xlim, ylim = ylim))
  }
  ku <- max(10L, min(k, floor(length(unique(dat$x)) / 3L)))
  mod <- mgcv::bam(y ~ s(x, k = ku), data = dat, weights = dat$w,
                   family = gaussian(), method = "fREML",
                   discrete = TRUE, nthreads = nthreads)
  newx <- seq(xlim[1], xlim[2], length.out = ndraw)
  pr   <- predict(mod, newdata = data.frame(x = newx), type = "response", se.fit = TRUE)
  pred <- data.frame(x = newx,
                     fit = as.numeric(pr$fit),
                     lo  = as.numeric(pr$fit - 1.96 * pr$se.fit),
                     hi  = as.numeric(pr$fit + 1.96 * pr$se.fit))
  list(mod = mod, pred = pred, xlim = xlim, ylim = ylim)
}

make_hex_trend_plot <- function(df, x_col, y_col, w_col,
                                x_lab, y_lab, title,
                                n_bins_x = NULL, n_bins_y = NULL,
                                min_n_per_bin = 25, nthreads = MGCV_THREADS,
                                show_ribbon95 = TRUE) {
  ribs  <- ribbon_summary(df, x_col, y_col, w_col, nbins = 140,
                          min_n_per_bin = min_n_per_bin)
  trend <- fit_bam_smooth_weighted(df, x_col, y_col, w_col,
                                   xlim = ribs$xlim, ylim = ribs$ylim,
                                   k = 80, nthreads = nthreads, ndraw = 600)
  dfp <- maybe_sample_df(df)
  n <- sum(is.finite(dfp[[x_col]]) & is.finite(dfp[[y_col]]))
  if (is.null(n_bins_x)) n_bins_x <- choose_hex_bins(n)
  if (is.null(n_bins_y)) n_bins_y <- n_bins_x
  p <- ggplot(dfp, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    stat_binhex(bins = c(n_bins_x, n_bins_y), aes(fill = after_stat(count))) +
    scale_fill_viridis_c(trans = "log10",
                         breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::label_number(accuracy = 1),
                         name   = "count (log10)") +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2, alpha = 0.8) +
    coord_cartesian(xlim = trend$xlim, ylim = trend$ylim, expand = FALSE) +
    geom_ribbon(data = ribs$summary, aes(x = x, ymin = q10, ymax = q90),
                inherit.aes = FALSE, alpha = 0.20) +
    { if (show_ribbon95)
      geom_ribbon(data = ribs$summary, aes(x = x, ymin = q05, ymax = q95),
                  inherit.aes = FALSE, alpha = 0.10) else NULL } +
    geom_line(data = trend$pred, aes(x = x, y = fit),
              inherit.aes = FALSE, linewidth = 1.1) +
    labs(title = title, x = x_lab, y = y_lab) +
    guides(fill = guide_colorbar(barheight = grid::unit(60, "pt")))
  bld <- ggplot_build(p)
  hex_max <- tryCatch(max(bld$data[[1]]$count, na.rm = TRUE), error = function(e) NA_real_)
  attr(p, "hex_max") <- hex_max
  p
}

bivar_mean_map <- function(df, x_col, y_col, w_col, z_col,
                           nx = NULL, ny = NULL, min_n = 50,
                           xlab = NULL, ylab = NULL, title = NULL,
                           lim_q = 0.98, fixed_limits = NULL) {
  df <- maybe_sample_df(df)
  dt <- data.table::as.data.table(df[, c(x_col, y_col, w_col, z_col)])
  setnames(dt, c("x","y","w","z"))
  dt <- dt[is.finite(x) & is.finite(y) & is.finite(w) & is.finite(z)]
  xlim <- central_limits(dt$x); ylim <- central_limits(dt$y)
  dt <- dt[x >= xlim[1] & x <= xlim[2] & y >= ylim[1] & y <= ylim[2]]
  n <- nrow(dt)
  if (is.null(nx)) nx <- choose_hex_bins(n)
  if (is.null(ny)) ny <- nx
  bx <- seq(xlim[1], xlim[2], length.out = nx + 1L)
  by <- seq(ylim[1], ylim[2], length.out = ny + 1L)
  dt[, gx := cut(x, breaks = bx, include.lowest = TRUE, labels = FALSE)]
  dt[, gy := cut(y, breaks = by, include.lowest = TRUE, labels = FALSE)]
  grd <- dt[, .(
    n  = .N,
    mu = if (.N > 0) matrixStats::weightedMean(z, w, na.rm = TRUE) else NA_real_
  ), by = .(gx, gy)][!is.na(gx) & !is.na(gy)]
  grd[, x := 0.5 * (bx[-1] + bx[-length(bx)])[gx]]
  grd[, y := 0.5 * (by[-1] + by[-length(by)])[gy]]
  grd <- grd[n >= min_n]
  L <- quantile(abs(grd$mu), probs = lim_q, na.rm = TRUE)
  fill_limits <- if (!is.null(fixed_limits)) fixed_limits else c(-L, L)
  ggplot(grd, aes(x, y, fill = mu, alpha = log10(pmax(n, 1)))) +
    geom_tile() +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_fill_distiller(type = "div", palette = "RdBu",
                         limits = fill_limits, oob = scales::squish,
                         name = "Mean error [m]") +
    scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    labs(title = if (is.null(title)) "" else title,
         x = if (is.null(xlab)) x_col else xlab,
         y = if (is.null(ylab)) y_col else ylab)
}

plot_hist_qq <- function(vec, title_prefix, bins = 60) {
  vec <- maybe_sample_vec(vec)
  xlim <- central_limits(vec)
  
  # Histogram (unchanged)
  p1 <- ggplot(data.frame(x = vec), aes(x = x)) +
    geom_histogram(bins = bins, alpha = 0.6) +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    labs(x = "Error (m)", y = "Count", title = paste0(title_prefix, " — histogram"))
  
  # Improved Q-Q plot with standardized axes
  # Remove NA values and sort
  vec_clean <- vec[!is.na(vec)]
  n <- length(vec_clean)
  sorted_vec <- sort(vec_clean)
  
  # Standardize the sample values
  sample_mean <- mean(sorted_vec)
  sample_sd <- sd(sorted_vec)
  sample_quantiles <- (sorted_vec - sample_mean) / sample_sd
  
  # Calculate theoretical quantiles for normal distribution
  probs <- ppoints(n)  # (1:n - 0.5)/n
  theoretical_quantiles <- qnorm(probs)
  
  # Create Q-Q data frame with standardized values
  qq_data <- data.frame(
    theoretical = theoretical_quantiles,
    sample = sample_quantiles
  )
  
  # Determine axis limits for equal aspect ratio
  lim_range <- max(abs(range(c(theoretical_quantiles, sample_quantiles)))) * c(-1, 1)
  
  # Create Q-Q plot with 45-degree reference line
  p2 <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.5) +  # 45-degree line
    geom_point(alpha = 0.5, size = 1) +
    coord_equal(xlim = lim_range, ylim = lim_range) +  # Equal axes for 45-degree line
    labs(
      x = "Theoretical Quantiles (Normal)",
      y = "Sample Quantiles (Standardized)",
      title = paste0(title_prefix, " — Q-Q (Normal)")
    ) +
    theme(aspect.ratio = 1)  # Ensure square plot
  
  p1 | p2
}


plot_ecdf_by_lc <- function(df, err_col, lc_col, top_classes, title) {
  dd <- df %>% dplyr::filter(.data[[lc_col]] %in% top_classes, is.finite(.data[[err_col]]))
  dd <- maybe_sample_df(dd)
  dd <- dd %>% dplyr::mutate(abs_err = abs(.data[[err_col]]))
  ggplot(dd, aes(x = abs_err, colour = .data[[lc_col]])) +
    stat_ecdf(geom = "step") +
    coord_cartesian(xlim = c(0, quantile(dd$abs_err, 0.99, na.rm = TRUE))) +
    labs(x = "|error| (m)", y = "ECDF", colour = "LC L1", title = title)
}

corr_plot <- function(df, cols, title, out_png) {
  dd <- df %>% dplyr::select(all_of(cols)) %>% na.omit()
  dd <- maybe_sample_df(dd)
  if (nrow(dd) > 2) {
    cm <- cor(dd, method = "pearson")
    pdf(sub("\\.png$", ".pdf", out_png), width = 10, height = 7.5)
    corrplot(cm, method = "color", tl.cex = 0.8, mar = c(0,0,1,0), title = title)
    dev.off()
  }
}

make_pairs <- function(df, cols, out_png, nmax = 20000) {
  set.seed(42)
  dd <- df %>% dplyr::select(all_of(cols)) %>% na.omit()
  if (nrow(dd) > nmax) dd <- dd %>% dplyr::slice_sample(n = nmax)
  p <- GGally::ggpairs(dd, progress = FALSE)
  ggsave(sub("\\.png$", ".pdf", out_png), p, width = 14, height = 14, bg = "white")
}

do_pca_plot <- function(df, cols, title, out_png) {
  dd <- df %>% dplyr::select(all_of(cols)) %>% na.omit()
  if (!nrow(dd)) return(invisible(NULL))
  dd_sample <- maybe_sample_df(dd)
  pr <- prcomp(dd_sample, center = TRUE, scale. = TRUE)
  var_expl <- pr$sdev^2 / sum(pr$sdev^2)
  pc <- as.data.frame(pr$x[, 1:2])
  p1 <- ggplot(pc, aes(PC1, PC2)) + 
    geom_point(size = 0.8, alpha = 0.6) +
    labs(title = sprintf("%s — PC1=%.1f%%, PC2=%.1f%%", 
                         title, 100*var_expl[1], 100*var_expl[2]))
  ggsave(sub("\\.png$", ".pdf", out_png), p1, width = 7, height = 5, bg = "white")
  invisible(pr)
}

rose_error_plot <- function(df, err_col, title) {
  dd <- df %>% 
    dplyr::select(meta_target_azimuth_avg, !!rlang::sym(err_col)) %>%
    dplyr::filter(is.finite(meta_target_azimuth_avg), is.finite(.data[[err_col]])) %>%
    dplyr::mutate(abs_err = abs(.data[[err_col]]),
                  az_bin  = cut(meta_target_azimuth_avg, 
                                breaks = seq(0,360,by=10), include.lowest=TRUE))
  dd <- maybe_sample_df(dd)
  dd2 <- dd %>% 
    dplyr::group_by(az_bin) %>%
    dplyr::summarize(mu = mean(abs_err, na.rm=TRUE), n = dplyr::n(), .groups="drop") %>%
    tidyr::separate_wider_delim(az_bin, delim=",", names=c("lo","hi"), too_many="merge") %>%
    dplyr::mutate(
      lo = as.numeric(gsub("\\(|\\[","",lo)),
      hi = as.numeric(gsub("\\]","",hi)),
      mid = 0.5*(lo+hi)
    )
  ggplot(dd2, aes(x = mid, y = mu)) +
    geom_col(width = 10, alpha=0.8) +
    coord_polar(start = 0) +
    scale_x_continuous(limits=c(0,360), breaks=seq(0,330,30)) +
    labs(x="Azimuth [deg]", y="Mean |error| (m)", title=title)
}

cyclic_trend_plot <- function(df, az_col, err_col, w_col, title, 
                              k = 24, nthreads = MGCV_THREADS) {
  dd <- df %>% 
    dplyr::select(any_of(c(az_col, err_col, w_col))) %>% 
    dplyr::filter(is.finite(.data[[az_col]]), is.finite(.data[[err_col]]))
  dd <- maybe_sample_df(dd)
  rng <- c(0, 360)
  dd$az <- (dd[[az_col]] %% 360)
  dd$y  <- dd[[err_col]]
  dd$w  <- if (w_col %in% names(dd)) dd[[w_col]] else 1.0
  if (nrow(dd) < 200) {
    return(ggplot(dd, aes(az, y)) + geom_point(alpha = 0.2) + 
             labs(title = title, x = "Azimuth [deg]", y = "Error (m)"))
  }
  mod <- mgcv::bam(y ~ s(az, bs="cc", k = k), data = dd, weights = dd$w,
                   family = gaussian(), method = "fREML",
                   discrete = TRUE, nthreads = nthreads,
                   knots = list(az = rng))
  newx <- data.frame(az = seq(0, 360, by = 1))
  pr <- predict(mod, newdata = newx, type = "response", se.fit = TRUE)
  ggplot() +
    geom_point(data = dd, aes(az, y), alpha = 0.05) +
    geom_line(aes(newx$az, pr$fit), linewidth = 1.1) +
    geom_ribbon(aes(newx$az, ymin = pr$fit - 1.96*pr$se.fit, 
                    ymax = pr$fit + 1.96*pr$se.fit), alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = 2) +
    coord_cartesian(xlim = c(0,360)) +
    labs(title = title, x = "Azimuth [deg]", y = "Error (m)")
}

vendor_long <- function(df) {
  keep <- c("v_GE01","v_WV01","v_WV02","v_WV03")
  dd <- df %>% 
    dplyr::select(any_of(keep)) %>% 
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::pivot_longer(-row, names_to = "vendor", values_to = "ratio") %>%
    dplyr::filter(is.finite(ratio)) %>% 
    dplyr::group_by(vendor) %>%
    dplyr::summarize(mean_ratio = mean(ratio, na.rm=TRUE), .groups="drop")
  dd$vendor <- factor(dd$vendor, levels = keep)
  dd
}

p_vendors <- function(vd, title) {
  ggplot(vd, aes(vendor, mean_ratio)) +
    geom_col() + ylim(0,1) + 
    labs(x="Vendor", y="Mean ratio", title=title)
}

missing_heat <- function(df, vars, title, outfile) {
  dd <- df %>% dplyr::select(all_of(vars))
  miss <- sapply(dd, function(x) mean(!is.finite(x)))
  plot_df <- data.frame(var = names(miss), frac_missing = as.numeric(miss))
  p <- ggplot(plot_df, aes(x=reorder(var, frac_missing), y=frac_missing)) +
    geom_col() + coord_flip() + scale_y_continuous(labels=scales::percent) +
    labs(x="", y="% missing (non-finite)", title=title)
  ggsave(sub("\\.png$", ".pdf", outfile), p, width=7, height=6, bg="white")
}