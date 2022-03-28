##==============================##
##    Functions used in p.Rmd   ##
## =============================##


#marginal prior
makeMargPrior1 <- function(I = 5, M = 10000, b){
  a1 <- matrix(rnorm((I - 1) * M, 0, b), ncol = (I - 1))
  a2 <- pnorm(t(apply(a1, 1, sort)))
  a3 <- cbind(rep(0, M), a2, rep(1, M))
  prob <- t(apply(a3, 1, diff))
  p <- as.data.frame.table(prob)
  colnames(p) <- c('b','rating','prob')
  p$b <- b
  return(p)
}

#plot differences
diff_plot <- function(data, cats = NULL, ylim = NULL,
                      xlab = NULL, ylab = NULL, title = NULL,
                      mid_col = "white", end_col = "black", ...){
  
  x <- 0:nrow(data)
  y <- c(0, data[,1] - data[,2])
  col <- rep(mid_col, length(x))
  col[c(1, length(col))] <- end_col
  
  p <- ggplot(data.frame(x, y), aes(x, y)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line() + 
    geom_point(shape = 21, size = 2, fill = col) +
    theme_classic() +
    coord_capped_cart(bottom='both', left = "both") +
    theme(axis.ticks.length = unit(.2, "cm"),
          axis.text = element_text(color = "black"),
          text = element_text(size = 10))
  
  if(!is.null(cats)){
    if(length(cats) != nrow(data))
      stop("Number of category labels must match number of categories")
    else{
      p <- p + 
        scale_x_continuous(labels = c("N/A", cats)) +
        theme(axis.text.x = element_text(angle = 45, 
                                         vjust = 1,
                                         hjust = 1))
    }
  }
  
  if(!is.null(xlab))
    p <- p + labs(x = xlab)
  
  if(!is.null(ylab))
    p <- p + labs(y = ylab)
  
  if(!is.null(title))
    p <- p + ggtitle(title)
  
  if(!is.null(ylim))
    p <- p + ylim(ylim)
  
  p <- p + theme(...)
  
  return(p)
}

# plot roc
roc_plot <- function(data, cats = NULL, 
                     xlab = NULL, ylab = NULL, title = NULL,
                     mid_col = "white", end_col = "black", ...){
  
  x <- c(0, data[,1])
  y <- c(0, data[,2])
  col <- rep(mid_col, length(x))
  col[c(1, length(col))] <- end_col
  
  p <- ggplot(data.frame(x, y), aes(x, y)) +
    geom_abline(slope = 1, linetype = "dashed") +
    geom_line() + 
    geom_point(shape = 21, size = 2, fill = col) +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    scale_y_continuous(breaks = seq(0, 1, .2)) +   
    theme_classic() +
    coord_capped_cart(bottom="both", left = "both") +
    theme(axis.ticks.length = unit(.2, "cm"),
          axis.text = element_text(color = "black"),
          text = element_text(size = 10))
  
  if(!is.null(xlab))
    p <- p + labs(x = xlab)
  
  if(!is.null(ylab))
    p <- p + labs(y = ylab)
  
  if(!is.null(title))
    p <- p + ggtitle(title)
  
  p <- p + theme(...)
  
  return(p)
}


pipeline <- function(y1, y2, b = c(1, 1/3), M = 2e5, small = rep(5e-4, 2)){
  I <- length(y1)
  val <- rbind(y1, y2)
  tot <- rowSums(val)
  prop <- val / tot
  m <- rowSums(prop * rbind(1:I, 1:I))
  cum <- t(apply(prop, 1, cumsum))
  bf <- likertBF(cbind(y1, y2), b = b, M = M, small = small)
  return(list(val = val,
              tot = tot,
              prop = prop,
              m = m,
              cum = cum,
              bf = bf,
              b = b))
}