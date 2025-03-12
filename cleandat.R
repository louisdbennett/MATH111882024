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
  'satisfied_teaching',
  'satisfied_feedback',
  'avg_entry_tariff',
  'spent_per_student',
  'added_value',
  'continuation',
  'estimated_gini',
  'ethnic_diversity_index',
  'students_staff_ratio',
  'Total'
)

for(col in scale_cols) {
  unis[, col] <- unis[, col] - mean(unis[, col], na.rm = TRUE)
}

unis.mod <- unis[, c('career_after_15_month', scale_cols, 'russell')]
# Task: Predicting graduate employability and the most predictive features
# Part 1: Basic Linear Regression Model

# define a formula here
formula <- career_after_15_month ~
  estimated_gini +
  ethnic_diversity_index +
  russell

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
  data = unis,
  family = 'Beta',
  control.family = list(hyper = prec.prior),
  control.fixed = beta.prior
)

summary(inlamod)
