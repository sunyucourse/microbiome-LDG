#!/usr/bin/env Rscript
# P0-2: direction + strength metrics per habitat, no extrapolation.
# Direction: predicted ratio 10deg/50deg (posterior-simulation CI) + linear slope within 10-50deg.
# Strength: amplitude (max/min over robust observed range) + deviance explained (latitude-only GAM).
# Latitude centered at 30deg for the linear slope; slope units = log-richness per 10 degrees.
suppressMessages(library(mgcv)); suppressMessages(library(MASS))

df <- read.csv("/mnt/results/tables/grid_table.csv")
df <- df[!is.na(df$richness) & df$richness > 0, ]
df$log_rich <- log(df$richness)
df$log_effort <- log(df$n_samples)

set.seed(1)
out <- list()
for (hab in sort(unique(df$habitat))) {
  d <- df[df$habitat == hab, ]
  if (nrow(d) < 30) next
  rng <- as.numeric(quantile(d$lat_abs, c(0.025, 0.975)))
  k_lat <- max(4, min(10, floor(nrow(d) / 30)))
  # full model (as in original analysis, no depth covariate)
  g <- gam(log_rich ~ s(lat_abs, bs = "tp", k = k_lat) + s(log_effort, bs = "tp", k = 4) +
             tech_frac_wgs, data = d, method = "REML")
  # latitude-only model for deviance explained
  g0 <- gam(log_rich ~ s(lat_abs, bs = "tp", k = k_lat), data = d, method = "REML")
  dev_lat <- summary(g0)$dev.expl
  dev_full <- summary(g)$dev.expl
  # posterior simulation on prediction grid within robust observed range
  xs <- seq(rng[1], rng[2], length.out = 200)
  nd <- data.frame(lat_abs = xs, log_effort = median(d$log_effort),
                   tech_frac_wgs = median(d$tech_frac_wgs))
  Xp <- predict(g, newdata = nd, type = "lpmatrix")
  sims <- mvrnorm(2000, coef(g), vcov(g))
  pr <- Xp %*% t(sims)   # 200 x 2000 log-scale predictions
  med <- apply(pr, 1, median)
  ampl <- exp(apply(pr, 2, function(p) max(p) - min(p)))
  i_peak <- which.max(med)
  # direction ratio at 10 vs 50 deg (only if robust range covers both; P1-7 no extrapolation)
  r10_50 <- r_lo <- r_hi <- NA
  if (rng[1] <= 10 && rng[2] >= 50) {
    i10 <- which.min(abs(xs - 10)); i50 <- which.min(abs(xs - 50))
    rr <- exp(pr[i10, ] - pr[i50, ])
    r10_50 <- median(rr); r_lo <- quantile(rr, 0.025); r_hi <- quantile(rr, 0.975)
  }
  # linear slope within common window 10-50 deg (centered at 30)
  dw <- d[d$lat_abs >= 10 & d$lat_abs <= 50, ]
  sl <- sl_se <- sl_p <- NA
  if (nrow(dw) >= 30) {
    m <- lm(log_rich ~ I(lat_abs - 30), data = dw)
    cf <- summary(m)$coefficients
    sl <- 10 * cf[2, 1]; sl_se <- 10 * cf[2, 2]; sl_p <- cf[2, 4]
  }
  out[[hab]] <- data.frame(
    habitat = hab, n_grids = nrow(d), obs_lo = round(rng[1], 1), obs_hi = round(rng[2], 1),
    lat_of_peak = round(xs[i_peak], 1),
    amplitude = round(median(ampl), 2), ampl_lo = round(quantile(ampl, 0.025), 2),
    ampl_hi = round(quantile(ampl, 0.975), 2),
    ratio_10v50 = round(r10_50, 2), ratio_lo = round(r_lo, 2), ratio_hi = round(r_hi, 2),
    slope_per10 = round(sl, 4), slope_se = round(sl_se, 4), slope_p = signif(sl_p, 3),
    dev_lat_only = round(dev_lat, 4), dev_full = round(dev_full, 4))
}
res <- do.call(rbind, out)
write.csv(res, "/workspace/data/processed/p02_gradient_metrics.csv", row.names = FALSE)
print(res)
