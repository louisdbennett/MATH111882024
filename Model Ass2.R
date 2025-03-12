
  rm(list = ls(all = TRUE))
# Load data
unis <- read.csv("C:/Users/aparn/Downloads/finalbis.csv")


# Ethnic Diversity Index Computation: https://www.ed-data.org/article/Ethnic-Diversity-Index

# Add Ethnic Diversity Index
c1 = 100 
c2 = -100*sqrt(5*(5-1))/(5-1)
unis$ethnic_diversity_index = c1 + c2*sqrt((unis$White.ethnic.group - 0.2)^2 + (unis$Black.ethnic.group - 0.2)^2
                                           + (unis$Asian.ethnic.group - 0.2)^2 + (unis$Mixed.ethnic.group - 0.2)^2
                                           + (unis$Other.ethnic.group - 0.2)^2)

# Add Russell Group Uni factor
unis$russell = FALSE
russell_codes = c("B32", "B78", "C05", "C15", "D86", "E56", "E84", "G28", "I50", 
                  "K60", "L23", "L41", "L72", "M20", "N21", "N84", "O33", "Q50", "Q75",
                  "S18", "S27", "U80", "W20", "Y50")
unis[unis$INSTITUTION_CODE %in% russell_codes, ]$russell = TRUE
unis$continuation <- as.numeric(unis$continuation)
# covert POLAR4 %ages to numbers
unis[, paste0('POLAR4.Q', 1:5)] <- unis[, paste0('POLAR4.Q', 1:5)] * unis$Total

clean_uni <- function(code) {
  select <- data.table::setDT(unis[unis$INSTITUTION_CODE == code, ])
  
  select.m <- data.table::melt(
    select, 
    id.vars = c('INSTITUTION_CODE'),
    measure.vars = paste0('POLAR4.Q', 1:5)
  )
  
  start <- data.table::data.table(INSTITUTION_CODE = code, variable = 'BASE', value = 0)
  
  select.m <- rbind(start, select.m)
  
  select.m <- select.m[order(value)]
  
  select.m[
    , `:=`(cum.value = cumsum(value) / sum(value), cum.prop = seq(0, 1, by = 0.2))
  ]
}

plot_lorentz_curve <- function(c.prop, c.value) {
  ggplot2::ggplot(mapping = ggplot2::aes(x = c.prop, y = c.value)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymax = c.prop, ymin = c.value), fill = "#007a3e", alpha = 0.2
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_abline(
      slope = 1, intercept = 0
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(limits = c(0, NA)) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank()
    )
}

gini <- function(dist) {
  dist <- sort(dist)
  n <- length(dist)
  w <- rep(1 / n, n) 
  w.c <- cumsum(w) / sum(w)
  dist.c <- cumsum(dist) / sum(dist)
  
  g <- t(dist.c[-1]) %*% w.c[-n] - t(dist.c[-n]) %*% w.c[-1]
  as.numeric(g)
}

estimate_q <- function(dist, q = 10) {
  o <- dist / sum(dist)
  
  fn <- function(params) {
    x <- pbeta(q = seq(0.2, 1, length.out = 5), params[1], params[2])
    t <- x - c(0, x[-length(x)])
    
    sum((o[-1] - t) ^ 2)
  }
  
  fit <- optim(c(1, 1), fn)
  
  x <- pbeta(q = seq(0, 1, length.out = q + 1), fit$par[1], fit$par[2])
  x - c(0, x[-length(x)])
}

unis$RAW_GINI <- sapply(unis$INSTITUTION_CODE, \(code) {
  cleaned <- clean_uni(code)
  gini(cleaned$value[-1])
})

unis$ESTIMATED_GINI <- sapply(unis$INSTITUTION_CODE, \(code) {
  cleaned <- clean_uni(code)
  q <- estimate_q(cleaned$value)
  gini(q[-1])
})

code <- 'Y75'
cleaned <- clean_uni(code)

plot_lorentz_curve(cleaned$cum.prop, cleaned$cum.value) +
  ggplot2::ggtitle(glue::glue("Lorentz curve of student proportion by POLAR4 quintile in {code}"))

est.c.val <- cumsum(estimate_q(cleaned$value))
c.prop <- seq(0, 1, length.out = length(est.c.val))

plot_lorentz_curve(c.prop, est.c.val) +
  ggplot2::ggtitle(glue::glue("Estimated Lorentz curve of student proportion by POLAR4 quintile in {code}"))

unis

# Task: Predicting graduate employability and the most predictive features
# Part 1: Basic Linear Regression Model
unis$career_prop = unis$career_after_15_month / 100
unis$satisfied_teaching_scaled=scale(unis$satisfied_teaching) 
unis$spent_per_student_scaled=scale(unis$spent_per_student)
unis$Men_scaled=scale(unis$Men)
unis$continuation_scaled=scale(unis$continuation)
unis$ESTIMATED_GINI_scaled=scale(unis$ESTIMATED_GINI)
unis$russell_scaled=scale(unis$russell)
unis$ethnic_diversity_index_scaled=scale(unis$ethnic_diversity_index)



prec.prior <- list(
  phi= list(prior = "loggamma", param = c(1, 0.01))
)


prior.beta <- list(mean = 0, prec = 0.0001, 
                   mean = 0, prec = 0.0001) 


library(INLA)
# 1. Define your formula
formula_inla <- career_prop ~ Men_scaled+satisfied_teaching_scaled + spent_per_student_scaled +
  continuation_scaled + ESTIMATED_GINI_scaled + russell_scaled + ethnic_diversity_index_scaled 


# 2. Fit the Beta model
res_inla_beta <- inla(
  formula_inla,
  family = "beta",             # Use the Beta likelihood
  data = unis,
  control.family = list(hyper = prec.prior),  # Beta precision prior
  control.fixed  = prior.beta,  
  control.predictor = list(compute = TRUE), 
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# 3. Summarize results
summary(res_inla_beta)


# 2. Fit the model
res_inla_gaussian <- inla(
  formula_inla,
  family = "gaussian",      # linear (normal) model
  data = unis,
  control.family=list(hyper=prec.prior),control.fixed=prior.beta,
  control.predictor = list(compute = TRUE),       # to get fitted values
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)  # model fit metrics
)

# 3. Summarize results
summary(res_inla_gaussian)

plot(res_inla_gaussian$marginals.fixed$`(Intercept)`, type ="l",xlab="x",ylab="Density", 
     main='Posterior density of beta0') 


library(INLA)


unis$career_prop = unis$career_after_15_month / 100
prec.prior <- list(
  phi= list(prior = "loggamma", param = c(1, 0.01))
)


prior.beta <- list(mean = 0, prec = 0.0001, 
                   mean = 0, prec = 0.0001) 



plot(res_inla_beta$marginals.fixed$`(Intercept)`, type ="l",xlab="x",ylab="Density", 
     main='Posterior density of beta0') 







