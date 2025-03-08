unis <- read.csv("data/finalbis.csv")
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

estimate_q <- function(dist, q = 10, nbsamp = 10000) {
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
  q <- estimate_q(cleaned$value, q = 100)
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
