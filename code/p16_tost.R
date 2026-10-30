#!/usr/bin/env Rscript
# P1-6: bootstrap lat2 for endo/ecto gut + TOST equivalence test.
# Equivalence margin: |delta lat2| < 0.10 declared equivalent (pre-specified:
# corresponds to a ratio difference of ~e^0.10 = 1.11-fold across 55 degrees).
suppressMessages(library(mgcv)); suppressMessages(library(MASS))
set.seed(42)

d <- read.csv("/workspace/data/processed/grid_h2_gut.csv")
d$log_rich <- log(d$richness); d$log_effort <- log(d$n_samples)

lat2_of <- function(x) {
  x$lat_z <- scale(x$lat_abs)[, 1]
  m <- lm(log_rich ~ lat_z + I(lat_z^2) + log_effort + tech_frac_wgs, data = x)
  coef(m)["I(lat_z^2)"]
}

e <- d[d$habitat == "gut_ectotherm", ]; en <- d[d$habitat == "gut_endotherm", ]
B <- 2000
diffs <- numeric(B)
for (b in 1:B) {
  eb <- e[sample(nrow(e), replace = TRUE), ]
  nb <- en[sample(nrow(en), replace = TRUE), ]
  diffs[b] <- tryCatch(lat2_of(nb) - lat2_of(eb), error = function(e) NA)
}
diffs <- diffs[!is.na(diffs)]
obs <- lat2_of(en) - lat2_of(e)
ci <- quantile(diffs, c(0.025, 0.975))
# TOST: two one-sided tests against margins +/-0.10 using bootstrap SE
se <- sd(diffs)
t1 <- (obs - (-0.10)) / se; t2 <- (obs - 0.10) / se
p_eq <- 1 - pnorm(t1, lower.tail = FALSE)  # placeholder, computed properly below
p_low <- pnorm((obs + 0.10) / se)          # H0: delta <= -0.10
p_high <- 1 - pnorm((obs - 0.10) / se)     # H0: delta >= +0.10
p_tost <- max(p_low, p_high)
# minimal detectable effect (MDE) at 80% power, alpha=0.05 two-sided
mde <- 2.8 * se
cat(sprintf("obs_delta=%.4f  boot_SE=%.4f  95CI=[%.4f, %.4f]\n", obs, se, ci[1], ci[2]))
cat(sprintf("TOST p=%.4f (margin 0.10) -> %s\n", p_tost, ifelse(p_tost < 0.05, "equivalent", "NOT equivalent / indeterminate")))
cat(sprintf("MDE (80%% power)=%.4f\n", mde))
write.csv(data.frame(obs_delta=obs, boot_se=se, ci_lo=ci[1], ci_hi=ci[2],
                     tost_p=p_tost, mde80=mde, n_ecto=nrow(e), n_endo=nrow(en)),
          "/workspace/data/processed/p16_tost.csv", row.names = FALSE)
