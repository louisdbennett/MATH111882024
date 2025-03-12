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