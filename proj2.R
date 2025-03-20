rm(list = ls(all = TRUE))

source('funcs.R')

# load data
unis <- read.csv("data/finalbis.csv")
unis$continuation <- as.numeric(unis$continuation)

# calculate ethnic diversity index
races <- c('Other', 'Mixed', 'Asian', 'Black', 'White')
unis$ethnic_diversity_index <- 1 - rowSums(unis[, paste0(races, '.ethnic.group')] ^ 2)

# ensure career after 15 is within 0-1
unis$career_after_15_month <- unis$career_after_15_month / 100

# add Russell group indicator
russell <- rep(FALSE, NROW(unis))
russell_codes <- c(
  "B32", "B78", "C05", "C15", 
  "D86", "E56", "E84", "G28", 
  "I50", "K60", "L23", "L41", 
  "L72", "M20", "N21", "N84", 
  "O33", "Q50", "Q75", "S18", 
  "S27", "U80", "W20","Y50"
)
russell[unis$INSTITUTION_CODE %in% russell_codes] <- TRUE
unis$russell <- russell

# estimate the gini coefficients for each uni
unis$estimated_gini <- sapply(unis$INSTITUTION_CODE, \(code) {
  cleaned <- clean_uni(code)
  q <- estimate_q(cleaned$value, q = 10000)
  gini(q[-1])
})

# scale all the covariates
scale_cols <- c(
  'Men',
  'estimated_gini',
  'ethnic_diversity_index'
)
for(col in scale_cols) {
  unis[, paste0(col, '_ctr')] <- unis[, col] - mean(unis[, col], na.rm = TRUE)
}

# subset data and create linear predictor
unis.mod <- unis[, c('career_after_15_month', paste0(scale_cols, '_ctr'), 'russell')]
formula <- career_after_15_month ~ estimated_gini_ctr + ethnic_diversity_index_ctr + Men_ctr + russell

# define prior distributions
prec.prior <- list(phi = list(prior = "loggamma", param = c(1, 0.01)))
beta.prior <- list(
  mean.intercept = 0, 
  prec.intercept = 0.01,
  mean = rep(0, 4), 
  prec = rep(0.01, 4)
)

# fit model in INLA
inlamod <- INLA::inla(
  formula = formula,
  data = unis.mod,
  family = 'Beta',
  control.family = list(hyper = prec.prior),
  control.fixed = beta.prior,
  control.predictor = list(compute = TRUE), 
  control.compute   = list(dic = TRUE, cpo = TRUE)
)

summary(inlamod)

coeffs <- inlamod$summary.fixed

ilogit(coeffs[rownames(coeffs) == '(Intercept)', 'mean'])
ilogit(
  coeffs[rownames(coeffs) == 'russellTRUE', 'mean'] + coeffs[rownames(coeffs) == '(Intercept)', 'mean']
) /
  ilogit(coeffs[rownames(coeffs) == '(Intercept)', 'mean'])

png(filename="plots/beta_margs.png", width = 8, height = 3, units = 'in', res = 1200)
par(mfrow = c(1, 3), mar = c(2.5, 2.5, 5, 2.5))
plot(
  INLA::inla.tmarginal(\(x) exp(0.1 * x), inlamod$marginals.fixed[['Men_ctr']]),
  type = 'l',
  main = expression(exp(0.1 * beta[1])),
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
plot(
  INLA::inla.tmarginal(\(x) exp(0.1 * x), inlamod$marginals.fixed[['estimated_gini_ctr']]),
  type = 'l',
  main = expression(exp(0.1 * beta[2])),
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
plot(
  INLA::inla.tmarginal(\(x) exp(0.1 * x), inlamod$marginals.fixed[['ethnic_diversity_index_ctr']]),
  type = 'l',
  main = expression(exp(0.1 * beta[3])),
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
dev.off()
par(mfrow = c(1, 1))

calc_metrics(career_after_15_month ~ russell)

# fitting with prediction
png(
  filename = "plots/predictions.png",
  width = 12,
  height = 4,
  units = 'in',
  res = 1200
)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
plot_predictions(
  inlamod,
  unis,
  formula = formula,
  variable_name = 'Men',
  xlab = '% Males',
  ylab = '% Positive Outcomes',
  cex.lab = 2.5,
  cex.axis = 1.5
)
plot_predictions(
  inlamod,
  unis,
  formula = formula,
  variable_name = 'estimated_gini',
  xlab = 'Gini',
  ylab = '% Positive Outcomes',
  cex.lab = 2.5,
  cex.axis = 1.5
)
plot_predictions(
  inlamod,
  unis,
  formula = formula,
  variable_name = 'ethnic_diversity_index',
  xlab = 'ELF',
  ylab = '% Positive Outcomes',
  cex.lab = 2.5,
  cex.axis = 1.5
)
dev.off()
par(mfrow = c(1, 1))
