#' this script generates example AUPR plots to explain 
#' the metric computation framework in the manuscript.

library(PRROC)

set.seed(101)
plot_fn = "results/aupr_example_plots_for_manuscript.pdf"

### open plot file
pdf(plot_fn)

for(tmp in 1:3){
  ### generate data
  n = 1000
  x = abs(rnorm(n=n, mean = 0, sd = 0.5))
  y = x * abs(rnorm(n=n))
  y = y / max(y)
  
  # binary label
  yb = y
  yb[yb<quantile(yb, 0.8)] = 0
  yb[yb>0] = 1
  
  ### visualize data
  # hist(x, breaks = 30)
  # hist(y)
  # hist(yb)
  # plot(x,y)
  
  ### precision-recall curve
  probj = pr.curve(
    scores.class0 = x,
    weights.class0 = y,
    curve = T,
    max.compute = T,
    min.compute = T,
    rand.compute = T
  )
  plot(
    probj,
    max.plot = TRUE,
    min.plot = T,
    fill.area = T,
    color = "red",
    main = "",
    xlab = "",
    ylab = "",
    ann = F,
    col.lab = "white",
    col.axis = "white",
    bg = "white"
  )
  
  ### precision-recall curve with binary classification
  probj = pr.curve(
    scores.class0 = x,
    weights.class0 = yb,
    curve = T,
    max.compute = T,
    min.compute = T,
    rand.compute = T
  )
  plot(
    probj,
    max.plot = TRUE,
    min.plot = T,
    fill.area = T,
    color = "red",
    main = "",
    xlab = "",
    ylab = "",
    ann = F,
    col.lab = "white",
    col.axis = "white",
    bg = "white"
  )
}


### close plot file
dev.off()
