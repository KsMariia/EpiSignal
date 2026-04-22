
########## RUNNING INSTRUCTIONS:

###### load this code via:
#source('.../EpiSignalEM.R', chdir = TRUE)
#load required packages: BAMMtools, modi, bayestestR, ggplot2, gridExtra, dplyr, tidyr, ggrepel, RColorBrewer

#### let X be a dataframe, genes in rows and samples in columns. Then run:
#output_EpiSig <- EpiSignalEM(X)
#summary(output_EpiSig) 
#plot(output_EpiSig)

############################

##########################################################
##   BELOW: HELPER FUNCTION 

## CV: cap_floor_ was missing!!!
cap_floor_ <- function(x, cap = Inf, floor = -Inf) {
  x[x > cap] <- cap
  x[x < floor] <- floor
  x
}

kullback_leibler_cont_appr <- function(p, q) {
  
  if (!is.numeric(p)) {
    stop("kullback_leibler_cont_appr(): `p` must be numeric", call. = FALSE)
  }
  if (!is.numeric(q)) {
    stop("kullback_leibler_cont_appr(): `p` must be numeric", call. = FALSE)
  }
  if (length(p) != length(q)) {
    stop("kullback_leibler_cont_appr(): `p` must be the same length as `q`", call. = FALSE)
  }
  cap_floor_((cap_floor_(sum(log(p) * p)) - cap_floor_(sum(log(q) * p))) / length(p))
}


##########################################################
##   BELOW: MAIN FUNCTION
##########################################################

EpiSignalEM <- function(data,xl=NULL,mu_0m=NULL,sigma_0m=NULL,mu_1m=NULL,sigma_1=NULL,conv_thresh=10^{-8}){

 
  ############# GENERAL PACKAGE AND DATA CHECK  
  if (!requireNamespace("BAMMtools")) {
    stop("Please install package BAMMtools")
  }
  if (!requireNamespace("modi")) {
    stop("Please install package modi")
  }
  if (!requireNamespace("bayestestR")) {
    stop("Please install package bayestestR")
  }
  if (!requireNamespace("ggplot2")) {
    stop("Please install package ggplot2")
  }
  if (!requireNamespace("gridExtra")) {
    stop("Please install package gridExtra")
  }
  if (!requireNamespace("dplyr")) {
    stop("Please install package dplyr")
  }
  if (!requireNamespace("tidyr")) {
    stop("Please install package tidyr")
  }
  if (!requireNamespace("ggrepel")) {
    stop("Please install package ggrepel")
  }
  if (!requireNamespace("RColorBrewer")) {
    stop("Please install package RColorBrewer")
  }
  if (!is.numeric(conv_thresh) | (conv_thresh <= 0) | (conv_thresh > 0.01))  {
    stop("EpiSignalEM(): `conv_thresh`, the convergence threshold, must be between 0 and 0.01", call. = FALSE)
  }
  

    data <- as.data.frame(data)
    
  
    allzero_data <- subset(data, rowSums(data, na.rm = TRUE) == 0)
    if(!is.null(xl)){
      xl <- xl[which(rowSums(data, na.rm = TRUE) != 0)]
    }
    data <- subset(data, rowSums(data, na.rm = TRUE) != 0)
  
  
  ############# AUTOMATED INITIALISATION 
  if(sum(sapply(list(xl,mu_0m,sigma_0m,mu_1m,sigma_1),is.null))>0){
  
    print("Determining Starting Values...")
    
    dim_dat_2 <- dim(data)[2]
    if(is.null(dim_dat_2)){
      dim_dat_2 <- 1
    }
    dim_dat_1 <- dim(data)[1]
    if(is.null(dim_dat_1)){
      dim_dat_1 <- 1
    }
    
    # Jenk's breaks (1-dim k-means clustering)
    br_data <- rep(0,dim_dat_2)
    for(i in 1:(dim_dat_2)){
      br_data[i]<-getJenksBreaks(data[,i],k=3)[2]
    }
    # matrix index for genes on/off based on the breaks
    gind_data <- matrix(0,nrow=dim_dat_1,ncol=dim_dat_2)
    for(i in 1:(dim_dat_2)){
      for(j in 1:(dim_dat_1)){
        gind_data[j,i]<-as.numeric(data[j,i]>br_data[i])
      }
    }
    # set starting values from clusters based on the breaks
    gxl <- gind_data
    xl <- rowMeans(gxl)
    mu_0m <- as.numeric(colSums((1 - xl)*data)/colSums(1-gxl)) 
    sigma_0m <- sqrt(apply(data,2,w=1-xl,weighted.var))
    mu_1m <- colSums(xl*data)/sum(xl)
    vVars <- apply(data,2,w=xl,weighted.var)
    sigma_1 <- sqrt(mean(vVars))
  }
  

   if(sum(is.numeric(c(xl,mu_0m,sigma_0m,mu_1m,sigma_1))==FALSE)>0){
     stop("EpiSignalEM(): all input variables must be numeric", call. = FALSE)
   }
    if(any(sigma_0m<0.01)){
      sigma_0m <- sigma_0m + 0.01
    } 
    if(any(mu_0m<0.05)){
      mu_0m <- mu_0m + 0.05
    } 
    
   if(any(c(mu_0m,sigma_0m,mu_1m,sigma_1)<0.0001)){
     stop("EpiSignalEM(): all distribution parameters must be positive", call. = FALSE)
   }
   if(any(xl<0)|any(xl>1)){
     stop("EpiSignalEM(): entries in 'xl' must be probabilities", call. = FALSE)
   }
   if (length(sigma_1) != 1) {
     stop("EpiSignalEM(): standard deviation `sigma_1` must be of length 1", call. = FALSE)
   }
   
   if (all(sapply(list(length(mu_0m),length(sigma_0m),length(mu_1m)), function(x) x == length(mu_1m)) )!=TRUE) {
     stop("EpiSignalEM(): means `mu_0m`, `mu_1m` and standard deviation `sigma_0m` must be vectors of equal length", call. = FALSE)
   }
   
   if ( length(xl)!=nrow(data) ) {
     stop("EpiSignalEM(): probability vector 'xl' must have the same number of entries as the data has rows", call. = FALSE)
   }
   
  
  print("Initialising EM...")
  Q <- 0
  Q[1] <- 0
  LH <- 0
  M <- dim(data)[2]
  
  mu_0m_track <- matrix(NA, ncol=M)
  sigma_0m_track <- matrix(NA, ncol=M)
  mu_1m_track <- matrix(NA, ncol=M)
  sigma_1_track <- vector()
  
  mu_0m_track[1,] <- mu_0m
  sigma_0m_track[1,] <-sigma_0m
  mu_1m_track[1,] <- mu_1m
  sigma_1_track[1] <-  sigma_1
  
  ####################################################################
  for (i in 1:M) {
    LH_a <- (1-xl)*dnorm(data[,i], mu_0m[i], sigma_0m[i],log=TRUE)
    LH_b <- (xl)*dnorm(dnorm(data[,i], mu_1m[i], sigma_1),log=TRUE)
    LH <- LH + sum(LH_a+LH_b)
  }
  
  Q[1] <- LH 
  ##print(LH)
  
  ####################################################################
  ################ E step 
  k <- 2
      repeat {
        xl<-rep(1, nrow(data))                          
        onemxl<-rep(1, nrow(data))
        for (i in 1:M) {
        
        xl <- xl * dnorm(data[,i], mu_1m[i], sigma_1)
        onemxl <- onemxl * dnorm(data[,i], mu_0m[i], sigma_0m[i]) 
        
        }
       
        de <- xl+onemxl
        de[de==0] <- 4.940656e-324 
        xl <- xl/(de)
      
      ################ M step            
        mu_0m <- (colSums((1 - xl)*data)/sum(1-xl)) # 
        vVar0=apply(data,2,w=1-xl,weighted.var)
        ##vVar0=vVar0*(sum(1-xl)-1)/sum(1-xl)
        sigma_0m <- sqrt(vVar0)
        mu_1m <- colSums(xl*data)/sum(xl)
        vVars <- apply(data,2,w=xl,weighted.var)
        ##vVars=vVars*(sum(xl)-1)/sum(xl)
        sigma_1 <- sqrt(mean(vVars))
        
        if (sum(mu_0m<mu_1m)!=length(mu_0m) ){
          print("EpiSignalEM(): mean signal for at least one gene is lower than the average noise", call. = FALSE)
        }
        
        mu_0m_track <- rbind(mu_0m_track,mu_0m)
        sigma_0m_track <- rbind(sigma_0m_track,sigma_0m)
        mu_1m_track <- rbind(mu_1m_track,mu_1m)
        sigma_1_track <- c(sigma_1_track,sigma_1)
        
   ####################################################################
        LH <- 0
        for (i in 1:M) {
          LH_a <- (1-xl)*dnorm(data[,i], mu_0m[i], sigma_0m[i],log=TRUE)
          LH_b <- (xl)*dnorm(dnorm(data[,i], mu_1m[i], sigma_1),log=TRUE)
          LH <- LH + sum(LH_a+LH_b)
          ##print(c(sum(LH_a),sum(LH_b)))
          ##LH_max <- pmax(LH_a, LH_b)
          ##log_sum_exp <- LH_max + log(exp(LH_a - LH_b) + exp(LH_b - LH_max))
          ##log_sum_exp <- LH_max + log(exp(LH_a - LH_max) + exp(LH_b - LH_max))
          ##LH <- LH + sum(log_sum_exp)
        }
        
        Q[k] <- LH
        ##print(LH)
        
        ################
        if (abs(Q[k]-Q[k-1])<conv_thresh){break}
        
        if(k%%10==0){
          print(c("Iteration:",k))
        }
        
        k <- k + 1
      }
  
     hard_thresh <- xl>=0.5
     data_hardsel <- as.matrix(data[hard_thresh, ])
     norm_hard_thresh <-  data_hardsel- mu_1m[col(data_hardsel)]
     
     samples <- list(likelihood=Q,mu_0m=mu_0m_track, sigma_0m=sigma_0m_track, mu_1m=mu_1m_track, sigma_1=sigma_1_track)
     estimates <- list(mu_0m=as.numeric(mu_0m_track[k,]), sigma_0m=as.numeric(sigma_0m_track[k,]), mu_1m=as.numeric(mu_1m_track[k,]), sigma_1=as.numeric(sigma_1_track[k]))
     
     res <- list(likelihood=Q[k], k=k, xl=xl, estimates=estimates, samples=samples, data=data, filtered_out=rownames(allzero_data),
                norm_hard_thresh=norm_hard_thresh,hard_thresh_index=hard_thresh)
     
     
     info <- list()
     info$model <- list(model_name = "episignal_em",
                   iterations = k,
                   number_genes = nrow(data),
                   number_replicates = ncol(data))
     info$object <- list(size=format(utils::object.size(res), "Mb"), date= as.character(Sys.time()))
     
     
     print(info)
     class(res) <- "episignal_em"
     invisible(res)
     
}





##########################################################
##   BELOW: CUSTOM SUMMARY FUNCTION

summary.episignal_em <- function(object, ...) {
  
  data <- object$data
  assignment_summary <- table(object$hard_thresh_index, dnn = "")
  
  n_assignment <- table(object$hard_thresh_index,dnn="")
  names(n_assignment) <- c("Nr. Noise Genes", "Nr. Signal Genes")
  
 
  kl_output <- matrix(NA,nrow=ncol(object$data),ncol=2)
  names(kl_output)
  overlap_coef <- rep(NA,ncol(object$data))
  
  kl_list <- rep(NA,1000)
  for(i in 1:ncol(object$data)){
    dens_data1 <- stats::density(object$data[!(object$hard_thresh_index),i])
    dens_data2 <- stats::density(object$data[(object$hard_thresh_index),i])
    overlap_coef[i] <- as.numeric(bayestestR::overlap(dens_data1$y,dens_data2$y))
    
   for (j in 1:500) {
    dens_sim1 <- stats::density(stats::rnorm(n=as.numeric(n_assignment[1]), object$estimates$mu_0m[i], object$estimates$sigma_0m[i]))
    dens_sim2 <- stats::density(stats::rnorm(n=as.numeric(n_assignment[2]), object$estimates$mu_1m[i], object$estimates$sigma_1))
    kl_list[j] <- kullback_leibler_cont_appr(dens_data1$y, dens_sim1$y)
    kl_list[j+500] <- kullback_leibler_cont_appr(dens_data2$y, dens_sim2$y)
   } 
    kl_output[i,1] <- mean(kl_list[1:500])
    kl_output[i,2] <- mean(kl_list[501:1000])
  }
  data_distr <- cbind(kl_output,overlap_coef)
  colnames(data_distr) <- c("KL-Divergence_noise_distributions_(input:estimated)","KL_Divergence_signal_distributions_(input:estimated)","overlap_percentage_(noise:signal)")
  rownames(data_distr) <- paste0(rep("sample ",nrow(data_distr)),1:nrow(data_distr))
  
  
  kl_output <- rep(NA,ncol(object$norm_hard_thresh))
  sw_test_stat <- rep(NA,ncol(object$norm_hard_thresh))
  sw_test_p <- rep(NA,ncol(object$norm_hard_thresh))
  
  kl_list <- rep(NA,1500)
    
  for(i in 1:ncol(object$norm_hard_thresh)){
    dens_data3 <- stats::density(object$norm_hard_thresh[,i])
    for (j in 1:1500) {
      dens_sim3 <- stats::density(stats::rnorm(n=as.numeric(n_assignment[2]), 0 , object$estimates$sigma_1))
      kl_list[j] <- kullback_leibler_cont_appr(dens_data3$y, dens_sim3$y)
  
    } 
    kl_output[i] <- mean(kl_list)
    sw <- stats::shapiro.test(dens_data3$y)
    sw_test_stat[i] <-sw$statistic
    sw_test_p[i] <- sw$p.value
  }
  norm_distr <- cbind(kl_output,sw_test_stat,sw_test_p)
  colnames(norm_distr) <- c("KL-Divergence_(normalised:zero-centred)","Shapiro-Wilks_Statistic_(normalised)","Shapiro-Wilks_pval_(normalised)")
  rownames(norm_distr) <- paste0(rep("sample ",nrow(norm_distr)),1:nrow(norm_distr))

  
  ll_info <- c(mean(object$samples$likelihood[2:length(object$samples$likelihood)]),
               stats::sd(object$samples$likelihood[2:length(object$samples$likelihood)]),
               stats::median(object$samples$likelihood[2:length(object$samples$likelihood)]))
  names(ll_info) <- c("mean", "sd", "median")
  
  
  
  summary_res <- list("filtered_out"= as.numeric(length(object$filtered_out)),
                       "noise_vs_signal"= as.numeric(assignment_summary),
                      "estimated_means_noise" = as.numeric(object$estimates$mu_0m),
                      "estimated_sds_noise" = as.numeric(object$estimates$sigma_0m),
                      "estimated_means_signal" = as.numeric(object$estimates$mu_1m),
                      "estimated_sd_signal" = as.numeric(object$estimates$sigma_1),
                      "input_data_summary" = data_distr,
                      "output_data_summary" = norm_distr,
                      "log_likelihood" = ll_info)
  
  cat("Number of Filtered Out Genes (no signal in any sample):\n")
  print(summary_res$filtered_out)
  cat("\n")
  
   
  cat("Number of Noise and Signal Genes:\n")
  print(summary_res$noise_vs_signal)
  cat("\n")
  
  cat("Estimated Means for Noise Genes per Sample:\n")
  print(summary_res$estimated_means_noise)
  cat("\n")
  
  cat("Estimated Sds for Noise Genes per Sample:\n")
  print(summary_res$estimated_sds_noise)
  cat("\n")
  
  cat("Estimated Means for Signal Genes per Sample:\n")
  print(summary_res$estimated_means_signal)
  cat("\n")
  
  cat("Estimated Sd for Signal Genes across Samples:\n")
  print(summary_res$estimated_sd_signal)
  cat("\n")
  
  cat("Input Data Characteristics per Sample:\n")
  print(summary_res$input_data_summary)
  cat("\n")
  
  cat("Normalised Data Characteristics per Sample:\n")
  print(summary_res$output_data_summary)
  cat("\n")
  
  cat("Log Likelihood:\n")
  print(summary_res$log_likelihood)
  cat("\n")
  
  invisible(summary_res)
}


##########################################################
##   BELOW: CUSTOM PLOTTING FUNCTION

plot.episignal_em <- function(x,
                                 ...) {
  
  data <- x$data
  norm_data <- x$norm_hard_thresh
  idx <- x$k
  
  c24 <- brewer.pal(ncol(x$data), "Paired")
  
  if(ncol(data)>12){print("Warning: Plots may become messy because there are too many samples (for e.g. colours no longer unique, space issues, etc.)")}
  
  Conv.plots <- function(){
    
    par(mar = c(5,5,4,10))
    matplot(x$samples$mu_0m[2:length(x$samples$likelihood),],type="l",col=c24,lwd=2,xlab="Iteration", main="Convergence Check",ylab="Estimates of Mean Noise Genes", panel.first = grid())
    par(xpd=TRUE)
    legend(par('usr')[2], par('usr')[4], legend = colnames(x$samples$mu_0m),fill=c24,cex=0.42,bty="n")

    matplot(x$samples$mu_1m[2:length(x$samples$likelihood),],type="l",col=c24,lwd=2,xlab="Iteration",main="Convergence Check",ylab="Estimates of Mean Signal Genes",  panel.first = grid())
    par(xpd=TRUE)
    legend(par('usr')[2], par('usr')[4], legend = colnames(x$samples$mu_1m),fill=c24,cex=0.42,bty="n")


    matplot(x$samples$sigma_0m[2:length(x$samples$likelihood),],type="l",col=c24,lwd=2,xlab="Iteration",ylab="Estimates of SD of Noise Genes", main="Convergence Check", panel.first = grid())
    par(xpd=TRUE)
    legend(par('usr')[2], par('usr')[4], legend = colnames(x$samples$sigma_0m),fill=c24,cex=0.42,bty="n")

    par(mar = c(5.1, 4.1, 4.1, 2.1))
    par(mfrow=c(1,1))
    par(xpd = FALSE)
    matplot(x$samples$sigma_1[2:length(x$samples$likelihood)],type="l",col="black",lwd=2,xlab="Iteration",ylab="Estimates of SD of Signal Genes", main="Convergence Check", panel.first = grid())
  
    plot(x$samples$likelihood[2:length(x$samples$likelihood)],type="l",col="black",lwd=2,xlab="Iteration",ylab="Log-Likelihood", main="Convergence Check", panel.first = grid())
    
    
    }
  Conv.plots()
  
  #
  
  fitplotlist <- list()
  ll <- length(x$estimates$mu_0m)
  for(j in 1:ll){
    dist <- rep("norm", 2)    
    params <- list()
    params[[1]] <- c(as.numeric(x$estimates$mu_0m[j]), x$estimates$sigma_0m[j])
    params[[2]] <- c(as.numeric(x$estimates$mu_1m[j]), x$estimates$sigma_1)
    weight <- c(1-sum(x$hard_thresh_index)/length(x$hard_thresh_index),sum(x$hard_thresh_index)/length(x$hard_thresh_index))
    resulting_mixture_of_normals <- mistr::mixdist(dist, params, weights = weight)
    
    nn <- colnames(x$data)[j]
    
    sample_data_long <- as.data.frame(x$data[,j])
    colnames(sample_data_long) <-"value"
    raw_hist <-ggplot(sample_data_long,aes(x = value))+geom_histogram(bins=50)+
      ggplot2::labs(x = "Observed Expression Level", 
                    y = "Counts") +
      ggtitle(paste0(nn))
    
   
    emp_data <- density(x$data[,j]) 
    pdfplot <- mistr::autoplot.dist(data=emp_data,resulting_mixture_of_normals,which="pdf")+
      ggplot2::labs(x = "Expression", 
                    y = "Density") +
      theme_bw()+
      theme(legend.position="bottom", legend.box = "horizontal")+
      ggtitle("Observed vs Fitted Densities (zoomed in)")+
      guides(colour=FALSE)+
      scale_fill_discrete(labels=c("fitted distribution","signal","noise"))
    
    
    qqplot <- mistr::QQplotgg(x$data[,j], resulting_mixture_of_normals, col = "black", line_col = "blue") +
      ggplot2::labs(x = "Fitted Mixture of Normals", 
                    y = "Data") +
      ggtitle("QQ-Plot")+
      ggplot2::theme_get()
    
    grid.arrange(raw_hist,pdfplot,qqplot)
    
  }
  
  #
  
  data_long <- as.data.frame(data) %>%
    mutate(ID = rownames(data)) %>% 
    pivot_longer(cols = -ID, names_to = "Sample", values_to = "Expression")
  
   original_density <- ggplot(data_long, aes(x = Expression))+
     geom_density(aes(fill=Sample),alpha = 0.23) +
     scale_fill_brewer(palette="Paired")+ 
    labs(title = "Original Data - Density Plots", x = "Expression Level", y = "Density")
   print(original_density)
   
  
  #
   
  eval_df_mean_noise <- data.frame(task_number = as.factor(c(replicate(ncol(data), "mean_noise"))),
                             Sample= paste0(colnames(x$data)),      
                             eval1 = c(apply(data[!(x$hard_thresh_index),],2,mean)-x$estimates$mu_0m),
                             eval2 = c((apply(data[!(x$hard_thresh_index),],2,mean)+x$estimates$mu_0m))/2)
  
  eval_df_mean_sig <- data.frame(task_number = as.factor(c(replicate(ncol(data), "mean_signal"))),
                                   Sample= paste0(colnames(x$data)),      
                                   eval1 = c(apply(data[(x$hard_thresh_index),],2,mean)-x$estimates$mu_1m),
                                   eval2 = c((apply(data[(x$hard_thresh_index),],2,mean)+x$estimates$mu_1m))/2)
  
  eval_df_sd_noise <- data.frame(task_number = as.factor(c(replicate(ncol(data), "sd_noise"))),
                                   Sample= paste0(colnames(x$data)),      
                                   eval1 = c(apply(data[!(x$hard_thresh_index),],2,sd)/x$estimates$sigma_0m))
  
  eval_df_sd_sig <- data.frame(task_number = as.factor(c(replicate(ncol(data), "sd_signal"))),
                                 Sample= paste0(colnames(x$data)),      
                                 eval1 = c(apply(data[(x$hard_thresh_index),],2,sd)/x$estimates$sigma_1))
  
  

  
  strip_mean_noise <- ggplot(eval_df_mean_noise, aes(x = eval2, y = eval1,color=Sample,label=Sample)) + geom_point(shape=4,size=3)+
    scale_color_brewer(palette="Paired")+  
    geom_hline(yintercept = 0,col="grey")+ 
     theme(legend.position="none")+
    labs(title = "Observed Data - Mean Noise", x = "Difference from Estimate", y = "Average with Estimate")+  
    geom_label_repel(aes(label = Sample),
                   box.padding   = 1, 
                   point.padding = 1,
                   size=1.2,
                   segment.color = 'grey50')
  
  strip_mean_signal <- ggplot(eval_df_mean_sig, aes(x = eval2, y = eval1,color=Sample)) + geom_point(shape=4,size=3)+
    scale_color_brewer(palette="Paired")+  
    geom_hline(yintercept = 0,col="grey")+ 
    theme(legend.position="none")+
    labs(title = "Observed Data - Mean Signal", x = "Difference from Estimate", y = "Average with Estimate")+
    geom_label_repel(aes(label = Sample),
                     box.padding   = 1, 
                     point.padding = 1,
                     size=1.2,
                     segment.color = 'grey50')
  
  strip_sd_noise <- ggplot(eval_df_sd_noise, aes(x = Sample, y = eval1,color=Sample)) + geom_point(shape=4,size=3)+
    scale_color_brewer(palette="Paired")+ 
    geom_hline(yintercept = 1,col="grey")+ 
    theme(legend.position="none")+
    labs(title = "Observed Data - Sd Noise", x = "", y = "Ratio vs Estimate")+
    geom_label_repel(aes(label = Sample),
                     box.padding   = 1, 
                     point.padding = 1,
                     size=1.2,
                     segment.color = 'grey50')+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  
  strip_sd_signal <- ggplot(eval_df_sd_sig, aes(x = Sample, y = eval1,color=Sample)) + geom_point(shape=4,size=3)+
    scale_color_brewer(palette="Paired")+ 
    geom_hline(yintercept = 1,col="grey")+ 
    theme(legend.position="none")+
    labs(title = "Observed Data - Sd Signal", x = "", y = "Ratio vs Estimate")+
    geom_label_repel(aes(label = Sample),
                     size=1.5,
                     box.padding   = 1, 
                     point.padding = 1,
                     segment.color = 'grey50')+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  

 grid.arrange(strip_mean_noise,strip_mean_signal,strip_sd_noise,strip_sd_signal)

  
 #
  
  Threshold.plots <- function(){
    par(mfrow=c(2,1))
    plot(sort(x$xl),xlab="Genes Sorted by Signal Probability",ylab="Probability Assignment to Signal",pch="x",
         panel.first= grid())
    abline(h=0.5)
    
    plot(x$xl,xlab="Genes In Input Order",ylab="Probability Assignment to Signal",pch="x",
         panel.first = grid())
    abline(h=0.5)
    par(mfrow=c(1,1))
  }
  Threshold.plots()

  # 
  thresh <- x$xl
  thresh <- replace(thresh,thresh<0.125,125)
  thresh <-replace(thresh,(thresh>=0.125)&(thresh<0.25),250)
  thresh <-replace(thresh,(thresh>=0.25)&(thresh<0.375),375)
  thresh <-replace(thresh,(thresh>=0.375)&(thresh<0.5),500)
  thresh <-replace(thresh,(thresh>=0.5)&(thresh<0.625),625)
  thresh <-replace(thresh,(thresh>=0.625)&(thresh<0.75),750)
  thresh <-replace(thresh,(thresh>=0.75)&(thresh<0.875),825)
  thresh <-replace(thresh,(thresh>=0.875)&(thresh<=1) ,1000)
  thresh <- thresh/1000
  cate <- table(thresh)
  
  data2 <- data
  data2$thresh <- as.character(thresh)
  thresh_mean <- aggregate(. ~ thresh, data2, FUN=mean)
  thresh_var <- aggregate(. ~ thresh, data2, FUN=var)
  
  data3 <- as.data.frame(cbind(rowSums(thresh_mean[,2:ncol(thresh_mean)])/(ncol(thresh_mean)-1),
                 sqrt(rowSums(thresh_var[,2:ncol(thresh_mean)])/(ncol(thresh_var)-1))))
  data3$thresh <- thresh_mean$thresh
  data3$cate <- cate
  
  gene_summary_plot <- ggplot(data3, aes(x = thresh, y = V1)) +
    geom_errorbar(aes(ymin = V1 - V2, ymax = V1 + V2), 
                  width = 0.2) +
    geom_point(size = 3) + 
    labs(title = "Stripchart of Average Expression (+/-SD) Across All Genes in Each Category", y = "Expression",x="Probability of Assignment to Signal (Rounded)")+
    geom_text(aes(label = paste0("n=",cate)),position=position_nudge(x=0.38,y=0.2))
   print(gene_summary_plot)
  
   #
   norm_data <- x$norm_hard_thresh
   sim_ref_signal_vals <- data.frame(rnorm(nrow(norm_data), mean=0, sd=x$estimates$sigma_1))
   colnames(sim_ref_signal_vals) <- c("Reference")
   norm_data2 <- cbind(norm_data,sim_ref_signal_vals)
   
   norm_long <- as.data.frame(norm_data2) %>%
     mutate(ID = rownames(norm_data)) %>%
     pivot_longer(cols = -ID, names_to = "Sample", values_to = "Expression")
  

  
   norm_sig_density <- ggplot(data=norm_long,aes(x = Expression, colour=Sample)) +
    stat_density(data=function(x){x[(x$Sample!="Reference"),]},geom="line",position="jitter") +
    scale_color_brewer(palette="Paired") +
    labs(title = "Normalised Signal - Density Plots", x = "Expression Level", y = "Density")+
    stat_density(data=function(x){x[x$Sample=="Reference",]},geom="line",col="black",linewidth=1.3)
  
  
  norm_sig_qq <- ggplot(norm_long, aes(sample = Expression, colour = Sample)) +
    stat_qq(data=function(x){x[(x$Sample!="Reference"),]}) +
    stat_qq_line(data=function(x){x[(x$Sample!="Reference"),]})+
    #scale_fill_manual(values = c24[1:ncol(x$data)])+
    scale_colour_brewer(palette="Paired")+ 
    labs(title = "Normalised Signal - QQ-Plot", x = "Theoretical Quantiles", y = "Sample Quantiles")
  
  
  eval_normsig <- data.frame(task_number = as.factor(c(replicate(ncol(data), "mean_difference"), 
                                              replicate(ncol(data), "sd_ratio_-1"))), 
                    Sample= rep(colnames(x$data),2),      
                    S = c(apply(norm_data,2,mean), (apply(norm_data,2,sd)/x$estimates$sigma_1)-1))
  
  eval_normsig_plot <- ggplot(eval_normsig, aes(x = task_number, y = S,color=Sample)) + geom_point(position="jitter",shape=4,size=3)+
    #scale_fill_manual(values = c24[1:ncol(x$data)])+ 
    scale_colour_brewer(palette="Paired") +
    labs(title = "Normalised Data Characteristics", x = "Summary Statistics", y = "Value")+  
    geom_hline(yintercept=0)
  

  grid.arrange(norm_sig_density,norm_sig_qq,eval_normsig_plot)
 

}



