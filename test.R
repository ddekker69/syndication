library(relevent)
library(dplyr)
library(lubridate)

## Paths: set FUNDING_INVESTOR_EL in .Renviron or edit the default below.
DATA_CSV <- Sys.getenv(
  "FUNDING_INVESTOR_EL",
  unset = "C:/Users/david/Dropbox/EBS-HW/Research/venture capital/syndication/funding_investor_edgelist.csv"
)

if (interactive()) {
  dn <- gsub("/", "\\\\", tcltk::tclvalue(tcltk2::tkchooseDirectory()))
  if (nzchar(dn)) setwd(dn)
}

funding_investor_el <- read.csv(DATA_CSV)
names(funding_investor_el)
el <- funding_investor_el[, 2:4]
nodes <- unique(c(el$Sender, el$Receiver))
el[["Sender"]] <- match(el$Sender, nodes)
el[["Receiver"]] <- match(el$Receiver, nodes)
names(el) <- c("time", "src", "dest")
dte <- as.Date(el$time)
base_time <- as.numeric(as.POSIXct(dte, tz = "CET"))
if (anyNA(base_time)) {
  stop("Funding.Date contains missing or invalid values; cannot build event times.")
}
# Within each day, enforce strict temporal ordering with millisecond offsets.
within_day_idx <- ave(seq_along(base_time), dte, FUN = seq_along) - 1
el$time <- base_time + within_day_idx * 1e-3
el <- el[order(el$time), ]
# Guard against any residual ties from numeric precision.
if (any(diff(el$time) <= 0)) {
  el$time <- cummax(el$time + seq_along(el$time) * 1e-6)
}
attr(el, "n") <- length(nodes)
fiel <- as.sociomatrix.eventlist(el, length(nodes))

## ---------------------------------------------------------------------------
## Model terms aligned with slide 9 ("Relational Event Model") in
## docs/every penny counts.pptx — four endogenous statistics (4 df in that
## table). Mapping to relevent::rem.dyad effect names:
##   Degree (# syndications) -> NTDegSnd
##   Recency                 -> RRecSnd
##   Referral                -> OSPSnd  (outbound shared partners)
##   Clustering              -> ISPSnd  (inbound shared partners)
## (OSPSnd/ISPSnd labels follow common dyadic closure readings; adjust if you
##  prefer two-path statistics OTPSnd/ITPSnd for either substantive row.)
## ---------------------------------------------------------------------------
effects_slide9 <- c("NTDegSnd", "RRecSnd", "OSPSnd", "ISPSnd")
# covar <- list(...)  # e.g. CovEvent with funding_investor_el$Is.Deadpooled

fit.ord <- rem.dyad(el, length(nodes), effects = effects_slide9, hessian = TRUE)
summary(fit.ord)

fit.time <- rem.dyad(el, length(nodes), effects = effects_slide9, ordinal = FALSE, hessian = TRUE)
summary(fit.time)
par(mfrow = c(1, 1))
plot(fit.ord$coef, fit.time$coef, asp = 1, xlab = "Ordinal BPM coef", ylab = "Temporal BPM coef")
abline(0, 1)

## BSIR: same data and effects as empirical fits (manual recommends larger sir.expand in practice).
fit.bsir <- rem.dyad(
  el,
  length(nodes),
  effects = effects_slide9,
  fit.method = "BSIR",
  sir.draws = 100,
  sir.expand = 5
)
summary(fit.bsir)

p <- length(fit.bsir$coef)
nr <- ceiling(sqrt(p))
nc <- ceiling(p / nr)
par(mfrow = c(nr, nc)) # e.g. four parameters -> 2x2 marginal histograms
for (i in seq_len(p)) {
  hist(fit.bsir$post[, i], main = names(fit.bsir$coef)[i], prob = TRUE, xlab = "Draw")
}

## Simulate from the fitted temporal model (structural effects only — no FESnd/FERec balance plots).
sim <- simulate(fit.time, nsim = 50000)
head(sim)

sim.pre <- sim[1:10, ]
sim2 <- simulate(fit.time, nsim = 20, edgelist = sim.pre)
all(sim2[1:10, ] == sim.pre)

sim2.t <- simulate(fit.time, nsim = 20, edgelist = sim.pre, redraw.timing = TRUE)
sim2.e <- simulate(fit.time, nsim = 20, edgelist = sim.pre, redraw.events = TRUE)

## ---------------------------------------------------------------------------
## Package manual example (relevent.pdf): FESnd+FERec on simulated 10-node data.
## Coefficient recovery plots apply only to that DGP. Set TRUE to run.
## ---------------------------------------------------------------------------
RUN_PACKAGE_DEMO <- FALSE

if (RUN_PACKAGE_DEMO) {
  roweff <- rnorm(10)
  roweff <- roweff - roweff[1]
  coleff <- rnorm(10)
  coleff <- coleff - coleff[1]
  lambda <- exp(outer(roweff, coleff, "+"))
  diag(lambda) <- 0
  ratesum <- sum(lambda)
  esnd <- as.vector(row(lambda))
  erec <- as.vector(col(lambda))
  time <- 0
  edgelist <- vector()
  while (time < 15) {
    drawsr <- sample(length(esnd), 1, prob = as.vector(lambda))
    time <- time + rexp(1, ratesum)
    if (time <= 15) {
      edgelist <- rbind(edgelist, c(time, esnd[drawsr], erec[drawsr]))
    } else {
      edgelist <- rbind(edgelist, c(15, NA, NA))
    }
  }

  effects_fe <- c("FESnd", "FERec")
  fit.ord.demo <- rem.dyad(edgelist, 10, effects = effects_fe, hessian = TRUE)
  summary(fit.ord.demo)
  par(mfrow = c(1, 2))
  plot(roweff[-1], fit.ord.demo$coef[1:9], asp = 1, xlab = "true", ylab = "est")
  abline(0, 1)
  plot(coleff[-1], fit.ord.demo$coef[10:18], asp = 1, xlab = "true", ylab = "est")
  abline(0, 1)

  fit.time.demo <- rem.dyad(edgelist, 10, effects = effects_fe, ordinal = FALSE, hessian = TRUE)
  fit.bsir.demo <- rem.dyad(
    edgelist,
    10,
    effects = effects_fe,
    fit.method = "BSIR",
    sir.draws = 100,
    sir.expand = 5
  )
  par(mfrow = c(3, 3))
  for (i in 1:9) {
    hist(fit.bsir.demo$post[, i], main = names(fit.bsir.demo$coef)[i], prob = TRUE)
    abline(v = roweff[i + 1], col = 2, lwd = 3)
  }
  for (i in 10:18) {
    hist(fit.bsir.demo$post[, i], main = names(fit.bsir.demo$coef)[i], prob = TRUE)
    abline(v = coleff[i - 8], col = 2, lwd = 3)
  }

  sim.demo <- simulate(fit.time.demo, nsim = 50000)
  par(mfrow = c(1, 2))
  esnd_demo <- exp(c(0, fit.time.demo$coef[1:9]))
  esnd_demo <- esnd_demo / sum(esnd_demo) * 5e4
  erec_demo <- exp(c(0, fit.time.demo$coef[10:18]))
  erec_demo <- erec_demo / sum(erec_demo) * 5e4
  plot(esnd_demo, tabulate(sim.demo[, 2]), xlab = "Expected Out-events", ylab = "Out-events")
  abline(0, 1, col = 2)
  plot(erec_demo, tabulate(sim.demo[, 3]), xlab = "Expected In-events", ylab = "In-events")
  abline(0, 1, col = 2)
}
