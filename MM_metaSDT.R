# INDIVIDUAL SIGNAL DETECTION ANALYSIS FOR TYPE II DATA.
# Maxime Maheu, 2014.

# TEST:
  # setwd("~/Documents/R/metaSDT")
  # source("MM_metaSDT.r")
  # example_data <- read.csv("example_data.csv", header = TRUE, sep = ";")
  # results <- MM_metaSDT(input_data = data, min_conf = 1, max_conf = 6, design = 2, coefficient = 1, typeI = 1, output = 0)
  # results

MM_metaSDT <- function(input_data, min_conf, max_conf, design = 2, coefficient = 1, typeI = 0, output = 0) {
  
  # Print title
  cat("\nType I and type II signal detection analysis.\n")
  cat("Copyright Maxime Maheu (2014).\n\n")
  
  # Check if the data input is a data frame
  if (is.data.frame(input_data) == FALSE) {cat("Warning! The input table is not a data frame. Please convert it.\n")}
  
  # Complete the data frame
  input_data <- na.omit(input_data)
  if (ncol(input_data) == 2) {colnames(input_data) <- c("Accuracy", "Confidence")}
  if (ncol(input_data) == 3) {colnames(input_data) <- c("Accuracy", "Confidence", "Label")}
  
  decimals <- 4

  # DESCRIPTIVE ANALYSIS
  
  ntrials <- nrow(input_data)
  conf_levels <- seq(min_conf, max_conf)
  nconf <- length(conf_levels)
  mean_conf <- mean(as.matrix(input_data[, 2]))
  sd_conf <- sd(as.matrix(input_data[, 2]))
  med_conf <- round(median(conf_levels))
  
  # Display some warnings
  if (ntrials < 100) {
    cat("Warning! The SDT analysis will only be made on less than 100 trials which can led to poor estimation of subject's actual metacognitive accuracy and, in some cases, to completly wrong computations.\n")}
  
  # CORRELATION ANALYSIS
  
  phi <- as.numeric(cor.test(as.matrix(input_data[, 2]), as.matrix(input_data[, 1]))[4])
  p_cor <- as.numeric(cor.test(as.matrix(input_data[, 2]), as.matrix(input_data[, 1]))[3])
  
  # PROBABILITY ANALYSIS
  
  trials_per_conf <- table(factor(input_data[, 1], levels = 0 : 1), factor(input_data[, 2], levels = min_conf : max_conf))
# You should probably make bins of your criterion (confidence ratings, control levels, ...) scale
  proportion_correct <- trials_per_conf[2, ] / table(factor(input_data[,2], levels = min_conf : max_conf)) # c
  subjective_probability <- seq(0.5, 1, length.out = nconf) # f
  trials_num_per_conf <- table(factor(input_data[, 2], levels = min_conf : max_conf))

  outcome <- (sum((proportion_correct * trials_num_per_conf), na.rm = TRUE) / ntrials) * (1 - (sum((proportion_correct * trials_num_per_conf), na.rm = TRUE) / ntrials))
  calibration <- sum((trials_num_per_conf * ((subjective_probability - proportion_correct) ^ 2)), na.rm = TRUE) / ntrials
  resolution <- sum((trials_num_per_conf * (proportion_correct - (sum(proportion_correct * trials_num_per_conf) / ntrials)) ^ 2), na.rm = TRUE) / ntrials
  
  PS <- outcome + calibration - resolution
  
  bias_coordinates <- matrix(c(subjective_probability, proportion_correct), nrow = 2, byrow = TRUE)
  bias_coordinates <- bias_coordinates[, complete.cases(t(bias_coordinates))]
  #if(nconf >= 10) {
  #  bias_coordinates <- t(subset(t(bias_coordinates), t(bias_coordinates)[,2] > 0))
  #  bias_coordinates <- t(subset(t(bias_coordinates), t(bias_coordinates)[,2] < 1))
  #  bias_coordinates <- t(subset(t(bias_coordinates), t(bias_coordinates)[,2] != 0.5))}
  
  # SIGNAL DETECTION ANALYSIS

  # Define the cut-off for first versus second part of trials
  cut_off <- round(ntrials/2)
  
  # Find a usefull index
  equation <- function(x) {min_conf + x - 1}
  index <- as.numeric(uniroot(equation, lower = -100000, upper = 100000)[1])
  
  # Launch the signal detection theory analysis  
  for (procedure in c(1, 2, 3)) {
      
      # Decide on which data subset working on
      if (procedure == 1) {data <- input_data[1 : cut_off,]}
      if (procedure == 2) {data <- input_data[cut_off + 1 : nrow(data),]}
      if (procedure == 3) {data <- input_data}      
      
      # Initialize some variables
      
      T1_beta0 <- NA
      T1_beta1 <- NA
      T1_r2 <- NA
      T1_p <- NA
      c <- NA
      d_prime <- NA
      
      T2_beta0 <- NA
      T2_beta1 <- NA
      T2_r2 <- NA
      T2_p <- NA
      meta_d_prime <- 0
      big_d_prime <- NA
      Aroc <- 0
    
      # If the type I analysis is required, do it. But this needs a SDT classification column.
      if (typeI == 1) {
        
        # Create a confidence ratings by answer types table
        table <- table(factor(data[, 3], levels = 1:4), factor(data[, 2], levels = min_conf : max_conf))
        typeI_table <- table + 0.5 
        
        # Initialize the first order probabilities variables
        p_hit <- vector()
        p_false_alarm <- vector()
        p_hit[1] <- 0
        p_false_alarm[1] <- 0
        
        # Compute the first order cumulative probabilities of hits and false alarms (used to plot the curve)
        i <- 2
        for (confidence in seq(max_conf, min_conf, by = -1)) {
          p_hit[i] <- (typeI_table[1, confidence + index] / (sum(typeI_table[1,]) + sum(typeI_table[3,]))) + p_hit[i - 1]
          p_false_alarm[i] <- (typeI_table[2, confidence + index] / (sum(typeI_table[2,]) + sum(typeI_table[4,]))) + p_false_alarm[i - 1]
          i <- i + 1}
        i = length(conf_levels) + 2
        for (confidence in conf_levels) {
          p_hit[i] <- (typeI_table[3, confidence + index] / (sum(typeI_table[1,]) + sum(typeI_table[3,]))) + p_hit[i - 1]
          p_false_alarm[i] <- (typeI_table[4, confidence + index] / (sum(typeI_table[2,]) + sum(typeI_table[4,]))) + p_false_alarm[i - 1]
          i <- i + 1}
        
        # Get categorical answers type
        p_hits <- sum(typeI_table[1,]) / (sum(typeI_table[1,]) + sum(typeI_table[3,]))
        p_false_alarms <- sum(typeI_table[2,]) / (sum(typeI_table[2,]) + sum(typeI_table[4,]))
        z_hits <- qnorm(p_hits, 0, 1)
        z_false_alarms <- qnorm(p_false_alarms, 0, 1)
        
        # Compute the type I d' and criterion
        if (design == 1) {
          d_prime <- z_hits - z_false_alarms}
        if (design == 2) {
          d_prime <- (z_hits - z_false_alarms) * (1 / sqrt(2))}
        c <- -(1 / 2) * (z_hits + z_false_alarms)
        
        # Fit linear regression to z transformed type I ROC (least squares fit)
        fit_typeI_cor <- cor.test(qnorm(p_hit[2 : (length(p_hit) - 1)], 0, 1), qnorm(p_false_alarm[2 : (length(p_false_alarm) - 1)], 0, 1))
        fit_typeI_poly <- lm(qnorm(p_hit[2 : (length(p_hit) - 1)], 0, 1) ~ qnorm(p_false_alarm[2 : (length(p_false_alarm) - 1)], 0, 1))
        T1_beta0 <- as.numeric(fit_typeI_poly$coefficients[2])
        T1_beta1 <- as.numeric(fit_typeI_poly$coefficients[1])
        T1_r2 <- as.numeric(fit_typeI_cor[4])^2
        T1_p <- as.numeric(fit_typeI_cor[3])}
      
      # Split into p(confidence|correct) and p(confidence|incorrect)
      typeII_table <- table(factor(data[, 1], levels = c(0, 1)), factor(data[, 2], levels = min_conf : max_conf)) + 0.5
      H <- typeII_table[2, (min_conf + index) : (max_conf + index)] + 0.5
      FA <- typeII_table[1, (min_conf + index) : (max_conf + index)] + 0.5

      # Get categorical answers type
      pH <- vector()
      pFA <- vector()
      for (confidence in conf_levels) {
        pH[confidence + index] <- H[confidence + index] / sum(H)
        pFA[confidence + index] <- FA[confidence + index] / sum(FA)}
      
      # Compute z score of type II hits and false alarms
      zH_cum <- qnorm(sum(pH[med_conf : max_conf]), 0, 1)
      zFA_cum <- qnorm(sum(pFA[med_conf : max_conf]), 0, 1)
      
      # Compute the type II d'
      if (design == 1) {
        meta_d_prime <- zH_cum - zFA_cum}
      if (design == 2) {
        meta_d_prime <- (zH_cum - zFA_cum) * (1 / sqrt(2))}
            
      # Compute the second order cumulative probabilities of hits and false alarms (used to plot the curve)
      pH_cum <- vector()
      pFA_cum <- vector()
      pH_cum[1] <- 0
      pFA_cum[1] <- 0
      for (confidence in seq((min_conf + index + 1), (max_conf + index + 1))) {
        pH_cum[confidence] <- pH[(max_conf + index + 2) - confidence] + pH_cum[confidence - 1]
        pFA_cum[confidence] <- pFA[(max_conf + index + 2) - confidence] + pFA_cum[confidence - 1]}

      # Fit linear regression to z transformed type II ROC (least squares fit)
      fit_typeII_cor <- cor.test(qnorm(pH_cum[2 : (length(pH_cum) - 1)], 0, 1), qnorm(pFA_cum[2 : (length(pFA_cum) - 1)], 0, 1))
      fit_typeII_poly <- lm(qnorm(pH_cum[2 : (length(pH_cum) - 1)], 0, 1) ~ qnorm(pFA_cum[2 : (length(pFA_cum) - 1)], 0, 1))
      T2_beta0 <- as.numeric(fit_typeII_poly$coefficients[2])
      T2_beta1 <- as.numeric(fit_typeII_poly$coefficients[1])
      T2_r2 <- as.numeric(fit_typeII_cor[4])^2
      T2_p <- as.numeric(fit_typeII_cor[3])
      
      # Calculate the ka and kb values
      ka <- vector()
      kb <- vector()
      for (confidence in seq((min_conf + index + 1), med_conf)) {
        ka[confidence - 1] <- (pH_cum[confidence] - pFA_cum[confidence - 1])^2 - (pH_cum[confidence - 1] - pFA_cum[confidence])^2}
      for (confidence in seq((med_conf + 1), (max_conf + index + 1))) {
        kb[confidence - med_conf] <- (pH_cum[confidence] - pFA_cum[confidence - 1])^2 - (pH_cum[confidence - 1] - pFA_cum[confidence])^2}
      ka <- sum(ka)
      kb <- sum(kb)
      
      # Compute the area under the type II ROC curve
      if (coefficient == 1) {
        Aroc <- (1 / 2) + ((1 / 4) * ka) + ((1 / 4) * kb)
        if (ka > 0 && kb > 0) {Bk <- log(((1/4) * ka) / ((1/4) * kb))}
        if (ka <= 0 | kb <= 0) {Bk <- -Inf}}
      if (coefficient == 2) {
        Aroc <- (1 / 2) + ((med_conf / (2 * length(conf_levels))) * ka) + (((med_conf + 1) / (2 * length(conf_levels))) * kb)
        if (ka > 0 && kb > 0) {Bk <- log(((med_conf / (2 * length(conf_levels))) * ka) / (((med_conf + 1) / (2 * length(conf_levels))) * kb))}
        if (ka <= 0 | kb <= 0) {Bk <- -Inf}}
      
      # Save splitted Aroc values
      if (procedure == 1) {first_half <- Aroc}
      if (procedure == 2) {second_half <- Aroc}
  }
  
  # FIRST PART OF THE RESULTS
  
  # Save the results in a matrix
  big_d_prime <- meta_d_prime / d_prime
  error <- first_half - second_half
  summary_table <- matrix(c(ntrials, mean_conf, sd_conf, phi, p_cor, outcome, calibration, resolution, PS, d_prime, c, meta_d_prime, big_d_prime, Bk, Aroc, error, T1_beta0, T1_beta1, T1_r2, T1_p, T2_beta0, T2_beta1, T2_r2, T2_p), nrow = 1, ncol = 24)
  colnames(summary_table) <- c("N", "Mean rate", "SD rate", "phi", "p", "O", "C", "R", "PS", "d'", "c", "meta-d'", "D'", "Bk", "Aroc", "Error", "T1 beta0", "T1 beta1", "T1 r2", "T1 p", "T2 beta0", "T2 beta1", "T2 r2", "T2 p")
  rownames(summary_table) <- ""
  
  if (output == 0) {
    
    # Print them
    cat("DESCRIPTIVE ANALYSIS\n N =", ntrials, "\n Mean rating =", round(mean_conf, decimals), "\n SD rating =", round(sd_conf, decimals), "\n")
    cat("\nCORRELATION ANALYSIS\n Phi =", round(phi, decimals), "\n p =", round(p_cor, decimals), "\n")
    cat("\nPROBABILITIES ANALYSIS\n O =", round(outcome, decimals), "\n C =", round(calibration, decimals), "\n R =", round(resolution, decimals), "\n PS =", round(PS, decimals), "\n")
    
    # FIRST PLOT
    
    plot.new()
    plot.window(xlim  = c(0.5, 1), ylim = c(0.5, 1), xaxs = "i", yaxs = "i")
    
    for (i in seq(0, 0.5, by = 0.001)) {
      segments(x0 = 0.5, y0 = (i + 0.5), x1 = (1 - i), y1 = 1, col = rgb(160/255, 032/255, 240/255, (i * 2) ^ 2))
      segments(x0 = (i + 0.5), y0 = 0.5, x1 = 1, y1 = (1 - i), col = rgb(255/255, 165/255, 000/255, (i * 2) ^ 2))}
    
    lines(x = bias_coordinates[1, ], y = bias_coordinates[2, ], col = "black", lty = 1, lwd = 5)
    lines(x = c(0.5, 1), y = c(0.5, 1), col = "black", lty = 1, lwd = 1)
    
    text(x = 0.725, y = 0.75, "No bias", srt = 45, col = "black")
    text(x = 0.6, y = 0.9, "Undercriterion", srt = 45, col = "purple")
    text(x = 0.9, y = 0.6, "Overcriterion", srt = 45, col = "orange")
    
    axis(1, xaxp = c(0.5, 1, 5))
    axis(2, yaxp = c(0.5, 1, 5))
    title(xlab = "Subjective probability", ylab = "Proportion correct", main = "Type I accuracy according to type II score")
    box(lwd = 2)
    
    # Wait before displaying the next plot
    cat ("\nPress [ENTER] to pursue the analysis.")
    line <- readline()
    
    # SECOND PART OF THE RESULTS
    
    cat("FIRST ORDER SDT\n d' =", round(d_prime, decimals), "\n c =", round(c, decimals), "\n beta0 = ", round(T1_beta0, decimals), "\n beta1 =", round(T1_beta1, decimals), "\n r2 =", round(T1_r2, decimals), "\n p =", round(T1_p, decimals), "\n")
    cat("\nSECOND ORDER SDT\n meta-d' =", round(meta_d_prime, decimals), "\n D' =", round(big_d_prime, decimals), "\n beta0 =", round(T2_beta0, decimals), "\n beta1 =", round(T2_beta1, decimals), "\n r2 =", round(T2_r2, decimals), "\n p =", round(T2_p, decimals), "\n")
    cat("\nROC ANALYSIS\n Bk =", round(Bk, decimals), "\n Aroc =", round(Aroc, decimals), "\n Error =", round(error, decimals), "\n\n")    
    
    # SECOND PLOT
    
    plot.new()
    plot.window(xlim  = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    
    if (typeI == 1) {
      # Plot the type I ROC curve
      polygon(x = p_false_alarm, y = p_hit, col = "grey")
      lines(x = p_false_alarm, y = p_hit, col = "black", lty = 1, lwd = 2)}
    
    # Plot the type II ROC curve
    polygon(x = pFA_cum, y = pH_cum, col = "royalblue3")
    lines(x = pFA_cum, y = pH_cum, col = "black", lty = 1, lwd = 2)
    
    # Plot the diagonals
    lines(x = c(0, 1), y = c(0, 1), col = "black", lty = 1, lwd = 1)
    lines(x = c(0, 0.5), y = c(1, 0.5), col = "black", lty = 2, lwd = 1)
    
    # Define some other graph parameters
    axis(1, xaxp = c(0, 1, 10))
    axis(2, yaxp = c(0, 1, 10))
    title(xlab = "P (rating | incorrect) or FA rate", ylab = "P (rating | correct) or H rate", main = "ROC curves")
    if (typeI == 1) {legend("bottomright", c("Type I ROC curve", "Type II ROC curve"), col = c("grey", "royalblue3"), pch = c(15, 15), bg = "white")}
    box(lwd = 2)
  }
  
  # Return the results table
  if (output == 0) {} # Display the results and plot the curves
  if (output == 1) {return(summary_table)} # Summary table
  if (output == 2) {return(bias_coordinates)} # Bias curve coordinates
  if (output == 3) {return(matrix(c(pFA_cum, pH_cum), nrow = 2, byrow = TRUE))} # ROC curve coordinates
}