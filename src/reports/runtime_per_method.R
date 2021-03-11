methods = c(
  "random",
  "clr",
  "wgcna",
  "aracne",
  "pcorr",
  "spice",
  "mrnetb",
  "mrnet",
  "glasso_likelihood",
  "et_genie3"
)

dir = "/work-zfs/abattle4/ashis/progres/spice_anlysis/gtex_v8/results/Brain_Cerebellum/corrected/AXPVAR/5000"
runtimes = sapply(methods, function(method){
  time_fn = sprintf("%s/%s_time.rds", dir, method)
  t = readRDS(time_fn)
  return(t[3])
})
x = data.frame(method = methods, runtime = round(as.numeric(runtimes)))
x$sec = sprintf("%s sec", round(x$runtime))
x$min = sprintf("%s min", round(x$runtime/60))
x$hour = sprintf("%0.1f hr", x$runtime/3600) 
x


