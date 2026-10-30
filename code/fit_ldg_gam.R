#!/usr/bin/env Rscript
# Fit per-habitat LDG models with mgcv on grid-level aggregated data.
# Usage: Rscript fit_ldg_gam.R <grid_csv> <out_prefix> [habitat_filter]
# Input columns expected: habitat, grid_id, glat, glon, lat_abs, richness (grid median
#   coverage-standardized), n_samples, depth_m (median, may be NA), tech_frac_wgs
# Outputs: <out_prefix>_curves.csv (predicted GAM curve per habitat),
#          <out_prefix>_strength.csv (per-habitat gradient strength metrics),
#          <out_prefix>_moransI.csv (residual spatial autocorrelation)
suppressMessages(library(mgcv))

args <- commandArgs(trailingOnly = TRUE)
in_csv <- args[1]
out_prefix <- args[2]
hab_filter <- if (length(args) >= 3) args[3] else NA
if (!is.na(hab_filter) && (hab_filter == "NA" || hab_filter == "")) hab_filter <- NA
force_no_depth <- if (length(args) >= 4) args[4] == "no_depth" else FALSE

df <- read.csv(in_csv)
if (!is.na(hab_filter)) df <- df[df$habitat == hab_filter, ]
df <- df[!is.na(df$richness) & df$richness > 0, ]
df$log_rich <- log(df$richness)
df$log_effort <- log(df$n_samples)
df$lat_z <- as.numeric(scale(df$lat_abs))

fit_one <- function(d, hab) {
  d <- d[order(d$lat_abs), ]
  n <- nrow(d)
  k_lat <- max(4, min(10, floor(n / 30)))
  has_depth <- !force_no_depth && sum(!is.na(d$depth_m)) >= 30 && length(unique(d$depth_m[!is.na(d$depth_m)])) > 5
  # GAM with smooth latitude + sampling effort (+ depth where available) + tech
  f <- if (has_depth) {
    log_rich ~ s(lat_abs, bs = "tp", k = k_lat) + s(log_effort, bs = "tp", k = 4) +
      s(log1p(depth_m), bs = "tp", k = 4) + tech_frac_wgs
  } else {
    log_rich ~ s(lat_abs, bs = "tp", k = k_lat) + s(log_effort, bs = "tp", k = 4) + tech_frac_wgs
  }
  g <- tryCatch(gam(f, data = d, method = "REML"), error = function(e) NULL)
  # quadratic parametric model for comparable gradient-strength coefficient
  q <- lm(log_rich ~ lat_z + I(lat_z^2) + log_effort + tech_frac_wgs, data = d)
  list(gam = g, quad = q, n = n, has_depth = has_depth)
}

habitats <- unique(df$habitat)
curves_all <- list()
strength_all <- list()
moran_all <- list()

for (hab in habitats) {
  d <- df[df$habitat == hab, ]
  if (nrow(d) < 30) next
  # drop NA-depth rows up front when depth will be used, so residuals align
  if (!force_no_depth && sum(!is.na(d$depth_m)) >= 30 && length(unique(d$depth_m[!is.na(d$depth_m)])) > 5) {
    d <- d[!is.na(d$depth_m), ]
  }
  fit <- fit_one(d, hab)
  # GAM curve prediction at fixed effort/tech (median values)
  if (!is.null(fit$gam)) {
    newd <- data.frame(
      lat_abs = seq(0, max(d$lat_abs), length.out = 120),
      log_effort = median(d$log_effort),
      tech_frac_wgs = median(d$tech_frac_wgs)
    )
    if (fit$has_depth) newd$depth_m <- median(d$depth_m, na.rm = TRUE)
    pr <- predict(fit$gam, newdata = newd, type = "response", se.fit = TRUE)
    curves_all[[hab]] <- data.frame(habitat = hab, lat_abs = newd$lat_abs,
                                    fit = as.numeric(pr$fit), se = as.numeric(pr$se.fit))
    # Moran's I on residuals (distance-based weights, k=8 nearest)
    res <- residuals(fit$gam)
    coords <- cbind(d$glon, d$glat)
    dd <- as.matrix(dist(coords))
    diag(dd) <- Inf
    k <- min(8, nrow(d) - 1)
    W <- apply(dd, 1, function(r) { o <- order(r)[1:k]; w <- numeric(nrow(d)); w[o] <- 1; w })
    W <- (W + t(W)) > 0
    W <- W / rowSums(W)
    rbar <- mean(res)
    moran <- (nrow(d) / sum(W)) * sum(W * outer(res - rbar, res - rbar)) / sum((res - rbar)^2)
    moran_all[[hab]] <- data.frame(habitat = hab, morans_I = moran, n = nrow(d))
  }
  # gradient strength from quadratic model
  cq <- summary(fit$quad)$coefficients
  lat2_est <- if ("I(lat_z^2)" %in% rownames(cq)) cq["I(lat_z^2)", "Estimate"] else NA
  lat2_se  <- if ("I(lat_z^2)" %in% rownames(cq)) cq["I(lat_z^2)", "Std. Error"] else NA
  # equator (|lat|=5) vs high-latitude (|lat|=55) predicted richness ratio from GAM curve
  eq_ratio <- NA; eq_ratio_lo <- NA; eq_ratio_hi <- NA
  if (!is.null(fit$gam) && exists("curves_all") && !is.null(curves_all[[hab]])) {
    cv <- curves_all[[hab]]
    if (max(cv$lat_abs) >= 55) {
      i5 <- which.min(abs(cv$lat_abs - 5))
      i55 <- which.min(abs(cv$lat_abs - 55))
      dfit <- cv$fit[i5] - cv$fit[i55]
      dse <- sqrt(cv$se[i5]^2 + cv$se[i55]^2)  # conservative: ignores covariance
      eq_ratio <- exp(dfit)
      eq_ratio_lo <- exp(dfit - 1.96 * dse)
      eq_ratio_hi <- exp(dfit + 1.96 * dse)
    }
  }
  strength_all[[hab]] <- data.frame(
    habitat = hab, n_grids = nrow(d),
    lat2_coef = lat2_est, lat2_se = lat2_se,
    lat2_t = lat2_est / lat2_se,
    eq_vs_highlat_ratio = eq_ratio, ratio_lo = eq_ratio_lo, ratio_hi = eq_ratio_hi,
    quad_r2 = summary(fit$quad)$adj.r.squared,
    gam_dev_expl = if (!is.null(fit$gam)) summary(fit$gam)$dev.expl else NA,
    has_depth_cov = fit$has_depth
  )
}

curves <- do.call(rbind, curves_all)
strength <- do.call(rbind, strength_all)
moran <- if (length(moran_all)) do.call(rbind, moran_all) else data.frame()
write.csv(curves, paste0(out_prefix, "_curves.csv"), row.names = FALSE)
write.csv(strength, paste0(out_prefix, "_strength.csv"), row.names = FALSE)
write.csv(moran, paste0(out_prefix, "_moransI.csv"), row.names = FALSE)
cat("GAM_FIT_DONE\n")
