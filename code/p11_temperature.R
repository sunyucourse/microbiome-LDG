#!/usr/bin/env Rscript
# P1-1: per-habitat temperature sensitivity.
# log_rich ~ temp (per 10 degC) + effort + tech; 1/kT activation-energy form;
# GAM(temp) deviance vs GAM(lat) deviance for comparison.
suppressMessages(library(mgcv))

gt <- read.csv("/mnt/results/tables/grid_table.csv")
cl <- read.csv("/workspace/data/processed/p11_grid_climate.csv")
d <- merge(gt, cl[, c("habitat","grid_id","bio1","bio4","bio12","sst","temp")],
           by = c("habitat","grid_id"))
d <- d[!is.na(d$richness) & d$richness > 0 & !is.na(d$temp), ]
d$log_rich <- log(d$richness); d$log_effort <- log(d$n_samples)
k_B <- 8.617e-5  # eV/K
d$inv_kT <- 1 / (k_B * (d$temp + 273.15))

out <- list()
for (hab in sort(unique(d$habitat))) {
  x <- d[d$habitat == hab, ]
  if (nrow(x) < 30 || length(unique(round(x$temp, 1))) < 10) next
  x$temp10 <- (x$temp - mean(x$temp)) / 10   # centered, per 10 degC
  m <- lm(log_rich ~ temp10 + log_effort + tech_frac_wgs, data = x)
  cf <- summary(m)$coefficients
  # 1/kT form: ln(rich) = c - E * (1/kT); E = -slope
  x$inv_kT_c <- x$inv_kT - mean(x$inv_kT)
  mk <- lm(log_rich ~ inv_kT_c + log_effort + tech_frac_wgs, data = x)
  cfk <- summary(mk)$coefficients
  # GAM deviance: temperature vs latitude
  k_s <- max(4, min(10, floor(nrow(x) / 30)))
  gt_ <- gam(log_rich ~ s(temp, bs = "tp", k = k_s) + s(log_effort, bs = "tp", k = 4) + tech_frac_wgs,
             data = x, method = "REML")
  gl_ <- gam(log_rich ~ s(lat_abs, bs = "tp", k = k_s) + s(log_effort, bs = "tp", k = 4) + tech_frac_wgs,
             data = x, method = "REML")
  out[[hab]] <- data.frame(
    habitat = hab, n_grids = nrow(x),
    temp_range = paste0(round(min(x$temp),1), " to ", round(max(x$temp),1)),
    slope_per10C = round(cf["temp10","Estimate"], 4),
    slope_se = round(cf["temp10","Std. Error"], 4),
    slope_p = signif(cf["temp10","Pr(>|t|)"], 3),
    E_eV = round(-cfk["inv_kT_c","Estimate"], 3),
    E_se = round(cfk["inv_kT_c","Std. Error"], 3),
    E_p = signif(cfk["inv_kT_c","Pr(>|t|)"], 3),
    dev_gam_temp = round(summary(gt_)$dev.expl, 4),
    dev_gam_lat = round(summary(gl_)$dev.expl, 4))
}
res <- do.call(rbind, out)
write.csv(res, "/workspace/data/processed/p11_temperature_sensitivity.csv", row.names = FALSE)
print(res)
