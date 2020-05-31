iter <- 100
dose <- cbindrnorm(100), 2=rnorm(100, 1, 3))
risk <- data.frame(happy=rnorm(100), sad=rnorm(100, 1, 3))
illness <- data.frame(happy=rnorm(100), sad=rnorm(100, 1, 3))

melt(dose)
plots <- list(
'1' = qmra.treatmentPlots(dose, risk, illness),
'2' = qmra.treatmentPlots(dose, risk, illness)
)
grid.draw(plots[[2]])
grid.arrange(grobs=plots, nrows=2, ncol=1, heights=unit(rep(4,5), rep('in', 5)))
