#rlda.multinomial<-function(data,){
  #Create a stop point
#  stopifnot(inherits(data, "matrix"))
  #Use a function not exported
#  Rlda:::lda_multinomial()
  #Create the class
#  class(res) <- c("list", "rlda")
#}

#plot.rlda <- function(d){
#  op <- par(mar = c(4, 4, 1, 1)) on.exit(par(op))
#  plot.new()
#  plot.window(xlim = range(d$x, na.rm = TRUE), ylim = range(d$y, na.rm = TRUE))
#  text(d$x, d$y, labels = d$labels)
#  axis(side = 1, range(d$x, na.rm = TRUE))
#  axis(side = 2, range(d$y, na.rm = TRUE))
#  invisible(d)
#}

#summary.rlda <-function(res){

#}
