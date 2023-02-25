fit.expr = function(DE, log = FALSE){
  expr = read.csv("../Data/fitted_DEmiRNA.csv", check.names = FALSE, row.names = 1)
  expr.DE = expr[,match(DE, colnames(expr))]
  
  expr.DE.df = cbind.data.frame(expr[,c(1:2)], expr.DE)
  expr.DE.df = cbind.data.frame(expr.DE.df[,1],
                                as.factor(expr.DE.df$Subject_Group),
                                lapply(expr.DE.df[,c(3:dim(expr.DE.df)[2])], as.numeric))
  colnames(expr.DE.df)[1:2] = c("Time", "Subject_Group")
  if(log == FALSE){
    expr.DE = expr.DE
    expr.DE.df = expr.DE.df
  }else{
    expr.DE = log2(expr.DE + 1)
    expr.DE.df = cbind.data.frame(expr.DE.df[,c(1:2)], expr.DE)
  }
  expr.DE.df$Bin = rep(seq(1, 20), 2) # creating bins for the fitted spline
  expr.DE.df = expr.DE.df[,c(1:2, dim(expr.DE.df)[2],
                             3:(dim(expr.DE.df)[2]-1))]
  expr.DE.df$Bin = paste0(expr.DE.df$Subject_Group, expr.DE.df$Bin)
  rownames(expr.DE.df) = expr.DE.df$Bin
  expr.DE.df = expr.DE.df[, -c(1:3)]
  expr.DE.df = as.data.frame(t(as.matrix(expr.DE.df)))
  expr.DE.m = as.matrix(expr.DE.df)
  return(list(expr.DE.df, expr.DE.m))
}
