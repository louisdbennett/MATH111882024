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

plot_predictions <- function(model, unis, formula, variable_name, n.grid = 50, ...) {
  var <- unis[[variable_name]]
  range <- seq(min(var), max(var), length.out = n.grid)

  # extract variables from the model object
  poss.vars <- gsub('_ctr', '', names(model$marginals.fixed))
  poss.vars <- poss.vars[poss.vars != "(Intercept)"]
  categorical.vars <- gsub('TRUE', '', poss.vars[endsWith(poss.vars, 'TRUE')])
  poss.vars <- poss.vars[!endsWith(poss.vars, 'TRUE')]

  newdata <- lapply(categorical.vars, \(x) {
    temp <- data.frame(
      x = c(rep(TRUE, length(range)), rep(FALSE, length(range)))
    )
    colnames(temp) <- x

    for(i in poss.vars) {
      temp[paste0(poss.vars, '_ctr')] <- 0
    }

    temp[variable_name] <- range
    temp[paste0(variable_name, '_ctr')] <- range - mean(var)

    temp
  }) |> 
    dplyr::bind_rows()

  newdata[is.na(newdata)] <- FALSE

  inlamod.pred <- INLA::inla(
    formula = formula,
    data = dplyr::bind_rows(newdata, unis),
    family = 'Beta',
    control.family = list(hyper = prec.prior),
    control.fixed = beta.prior,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(config = TRUE)
  )
  
  pred <- inlamod.pred$summary.fitted.values$mean

  plot(x = var, y = unis$career_after_15_month, ...)
  lines(range, pred[1:length(range)], col = "red", lwd = 3)
  lines(
    range,
    pred[(length(range) + 1):(2 * length(range))],
    col = "blue",
    lwd = 3
  )
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

ilogit <- function(x) 1 / (1 + exp(-x))

calc_metrics <- function(formula) {
  inlamod.form <- INLA::inla(
    formula = formula,
    data = unis.mod,
    family = 'Beta',
    control.family = list(hyper = prec.prior),
    control.fixed = beta.prior,
    control.predictor = list(compute = TRUE), 
    control.compute   = list(dic = TRUE, cpo = TRUE)
  )
  
  list(nslcpo = -sum(log(inlamod.form$cpo$cpo)), dic = inlamod.form$dic$dic)
}