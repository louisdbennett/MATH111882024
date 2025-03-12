rm(list = ls(all = TRUE))

source('funcs.R')

# Load data
unis <- read.csv("data/finalbis.csv")
unis$continuation <- as.numeric(unis$continuation)
# covert POLAR4 %ages to numbers
unis[, paste0('POLAR4.Q', 1:5)] <- unis[, paste0('POLAR4.Q', 1:5)] * unis$Total

# Add Ethnic Diversity Index
c1 <- 100
c2 <- -100 * sqrt(5 * (5 - 1)) / (5 - 1)
unis$ethnic_diversity_index <- c1 +
  c2 *
    sqrt(
      (unis$White.ethnic.group - 0.2)^2 +
        (unis$Black.ethnic.group - 0.2)^2 +
        (unis$Asian.ethnic.group - 0.2)^2 +
        (unis$Mixed.ethnic.group - 0.2)^2 +
        (unis$Other.ethnic.group - 0.2)^2
    )

# ensure career after 15 is within 0-1
unis$career_after_15_month <- unis$career_after_15_month / 100

# Add Russell Group Uni factor
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

# estimate the gini coefficient for each uni
unis$estimated_gini <- sapply(unis$INSTITUTION_CODE, \(code) {
  cleaned <- clean_uni(code)
  q <- estimate_q(cleaned$value, q = 10000)
  gini(q[-1])
})

# scale all the covariates as this makes the interpration of coefficients later better
scale_cols <- c(
  'Men',
  # 'satisfied_teaching',
  # 'satisfied_feedback',
  # 'avg_entry_tariff',
  # 'spent_per_student',
  # 'added_value',
  # 'continuation',
  'estimated_gini',
  'ethnic_diversity_index'
  # 'students_staff_ratio',
  # 'Total'
)

for(col in scale_cols) {
  unis[, paste0(col, '_ctr')] <- unis[, col] - mean(unis[, col], na.rm = TRUE)
}

unis.mod <- unis[, c('career_after_15_month', paste0(scale_cols, '_ctr'), 'russell')]
# Task: Predicting graduate employability and the most predictive features
# Part 1: Basic Linear Regression Model

# define a formula here
formula <- career_after_15_month ~ estimated_gini_ctr + ethnic_diversity_index_ctr + Men_ctr + russell

# define priors here
prec.prior <- list(phi = list(prior = "loggamma", param = c(1, 0.01)))

beta.prior <- list(
  mean.intercept = 0, 
  prec.intercept = 0.01,
  mean = 0, 
  prec = 0.01
)

inlamod <- INLA::inla(
  formula = formula,
  data = unis.mod,
  family = 'Beta',
  control.family = list(hyper = prec.prior),
  control.fixed = beta.prior,
  control.predictor = list(compute = TRUE), 
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

summary(inlamod)

plot(inlamod$marginals.fixed[['(Intercept)']], type = 'l')
plot(inlamod$marginals.fixed[['estimated_gini_ctr']], type = 'l')
plot(inlamod$marginals.fixed[['Men_ctr']], type = 'l')
plot(inlamod$marginals.fixed[['ethnic_diversity_index_ctr']], type = 'l')

# fitting with predictions
ethnic_diversity_index <- seq(0, 100, 2.5)
ethnic_diversity_index_ctr <- ethnic_diversity_index - mean(unis$ethnic_diversity_index)

gini <- seq(0, 0.5, 0.01)
gini_ctr <- gini - mean(unis$gini)

Men <- seq(0.2, 0.8, 0.05)
Men_ctr <- Men - mean(unis$Men)

newdata <- rbind(
  data.frame(
    career_after_15_month = NA,
    russell = TRUE,
    estimated_gini = 0,
    estimated_gini_ctr = 0,
    ethnic_diversity_index = 0,
    ethnic_diversity_index_ctr  = 0,
    Men = Men,
    Men_ctr = Men_ctr 
  ),
  data.frame(
    career_after_15_month = NA,
    russell = FALSE,
    estimated_gini = 0,
    estimated_gini_ctr = 0,
    ethnic_diversity_index = 0,
    ethnic_diversity_index_ctr = 0,
    Men = Men,
    Men_ctr = Men_ctr
  )
)

inlamod.pred <- INLA::inla(
  formula = formula,
  data = dplyr::bind_rows(newdata, unis.mod),
  family = 'Beta',
  control.family = list(hyper = prec.prior),
  control.fixed = beta.prior,
  control.predictor = list(compute = TRUE, link = 1),
  control.compute = list(config = TRUE)
)

pred <- inlamod.pred$summary.fitted.values$mean

plot(x = unis$Men, y = unis$career_after_15_month)
lines(Men, pred[1:length(Men)], col="red",lwd=3)
lines(Men, pred[(length(Men) + 1):(2 * length(Men))], col="blue",lwd=3)
