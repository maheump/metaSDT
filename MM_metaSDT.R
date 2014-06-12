# INDIVIDUAL SIGNAL DETECTION ANALYSIS FOR TYPE II DATA.
# Maxime Maheu, 2014.

# setwd("~/Documents/R/metaSDT")
# source("MM_metaSDT.r")
# example_data <- read.csv("example_data.csv", header = TRUE, sep = ";")
# results <- MM_metaSDT(example_data, 1, 6)
# results

MM_metaSDT <- function(input_data, min_conf, max_conf, design = 2, coefficient = 1) {

  # Print title
  cat("\nType I and type II signal detection analysis\n")
  cat("Copyright Maxime Maheu (2014)\n")
  
  # SUPPRIMER LES input_data$* ET LES REMPLACER PAR LE NUMÉRO DE LA COLONNE CORRESPONDANTE
  
  # Check if the data input is a data frame
  if (is.data.frame(input_data) == FALSE) {cat("Warning! The input table is not a data frame. Please convert it.\n")}
  
  # Complete the data frame
  input_data <- na.omit(input_data)
  for (trial in seq(1, nrow(input_data))) {
    if (input_data[trial, 1] == 1) {input_data[trial, 3] = 1}
    if (input_data[trial, 1] == 2) {input_data[trial, 3] = 0}
    if (input_data[trial, 1] == 3) {input_data[trial, 3] = 0}
    if (input_data[trial, 1] == 4) {input_data[trial, 3] = 1}}
  colnames(input_data) <- c("Label", "Confidence", "Outcome")
  
  decimals <- 4

  # DESCRIPTIVE ANALYSIS
  
  ntrials <- nrow(input_data)
  conf_levels <- seq(min_conf, max_conf)
  nconf <- length(conf_levels)
  mean_conf <- mean(as.matrix(input_data[2]))
  sd_conf <- sd(as.matrix(input_data[2]))
  
  # Display some warnings
  if (ntrials < 100) {
    cat("Warning! The SDT analysis will only be made on less than 100 trials which can led to poor estimation of subject's actual metacognitive accuracy and, in some cases, to completly wrong computations.\n")}
  
  # CORRELATION ANALYSIS
  
  phi <- as.numeric(cor.test(as.matrix(input_data[2]), as.matrix(input_data[3]))[4])
  p_cor <- as.numeric(cor.test(as.matrix(input_data[2]), as.matrix(input_data[3]))[3])
  
  # PROBABILITY ANALYSIS
  
  trials_per_conf <- table(factor(input_data$Outcome, levels = 0 : 1), factor(input_data$Confidence, levels = min_conf : max_conf))
  proportion_correct <- trials_per_conf[2,] / table(factor(input_data$Confidence, levels = min_conf : max_conf)) # c
  subjective_probability <- ((conf_levels * 0.5) / max_conf) + 0.5 # f
  trials_num_per_conf <- table(factor(input_data$Confidence, levels = min_conf : max_conf))
  
  outcome <- (sum(proportion_correct * trials_num_per_conf) / ntrials) * (1 - (sum(proportion_correct * trials_num_per_conf) / ntrials))
  calibration <- sum(trials_num_per_conf * ((subjective_probability - proportion_correct) ^ 2)) / ntrials
  resolution <- sum(trials_num_per_conf * (proportion_correct - (sum(proportion_correct * trials_num_per_conf) / ntrials)) ^ 2) / ntrials
  
  PS <- outcome + calibration - resolution

  # SIGNAL DETECTION ANALYSIS

  # Define the cut-off for first versus second part of trials
  cut_off <- round(ntrials/2)
  
  # Find a usefull index
  equation <- function(x) {min_conf + x - 1}
  index <- as.numeric(uniroot(equation, lower = -100000, upper = 100000)[1])

  # Launch the signal detection theory analysis  
  for (procedure in c(1, 2, 3)) {
      
      # Decide on which data subset work on
      if (procedure == 1) {data <- input_data[1 : cut_off,]}
      if (procedure == 2) {data <- input_data[cut_off + 1 : nrow(data),]}
      if (procedure == 3) {data <- input_data}
      
      # Initialize some variables
      d_prime <- 0
      meta_d_prime <- 0
      big_d_prime <- 0
      Aroc <- 0
    
      # Create a confidence ratings by answer types table
      table <- table(factor(data$Label, levels = 1:4), factor(data$Confidence, levels = min_conf : max_conf))
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
      i = length(seq(min_conf, max_conf)) + 2
      for (confidence in seq(min_conf, max_conf)) {
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
      
      # Split into p(confidence|correct) and p(confidence|incorrect)
      hit <- table[1, (min_conf + index) : (max_conf + index)]
      false_alarm <- table[2, (min_conf + index) : (max_conf + index)]
      miss <- table[3, (min_conf + index) : (max_conf + index)]
      correct_reject <- table[4, (min_conf + index) : (max_conf + index)]
      H <- hit + correct_reject + 0.5
      FA <- false_alarm + miss + 0.5
      
      # Get categorical answers type
      pH <- vector()
      pFA <- vector()
      for (confidence in seq(min_conf, max_conf)) {
        pH[confidence + index] <- H[confidence + index] / sum(H)
        pFA[confidence + index] <- FA[confidence + index] / sum(FA)}
      
      # Compute z score of type II hits and false alarms
      med_conf <- round(median(seq(min_conf, max_conf)))
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
        Bk <- log(((1/4) * ka) / ((1/4) * kb))}
      if (coefficient == 2) {
        Aroc <- (1 / 2) * ((med_conf / (2 * length(seq(min_conf, max_conf)))) * ka) * (((med_conf + 1) / (2 * length(seq(min_conf, max_conf)))) * kb)
        Bk <- log(((med_conf / (2 * length(seq(min_conf, max_conf)))) * ka) / (((med_conf + 1) / (2 * length(seq(min_conf, max_conf)))) * kb))}
      
      # Save splitted Aroc values
      if (procedure == 1) {first_half <- Aroc}
      if (procedure == 2) {second_half <- Aroc}
  }
  
  # FIRST PART OF THE RESULTS
  
  # Save the results in a matrix
  big_d_prime <- meta_d_prime / d_prime
  error <- first_half - second_half
  summary_table <- matrix(c(ntrials, mean_conf, sd_conf, phi, p_cor, outcome, calibration, resolution, PS, d_prime, c, meta_d_prime, big_d_prime, Bk, Aroc, error, as.numeric(fit_typeI_poly$coefficients[2]), as.numeric(fit_typeI_poly$coefficients[1]), as.numeric(fit_typeI_cor[4])^2, as.numeric(fit_typeI_cor[3]), as.numeric(fit_typeII_poly$coefficients[2]), as.numeric(fit_typeII_poly$coefficients[1]), as.numeric(fit_typeII_cor[4])^2, as.numeric(fit_typeII_cor[3])), nrow = 1, ncol = 24)
  colnames(summary_table) <- c("N", "Mean conf.", "SD conf.", "phi", "p", "O", "C", "R", "PS", "d'", "c", "meta-d'", "D'", "Bk", "Aroc", "Error", "T1 beta0", "T1 beta1", "T1 r2", "T1 p", "T2 beta0", "T2 beta1", "T2 r2", "T2 p")
  rownames(summary_table) <- ""
  
  # Print them
  cat("\nDESCRIPTIVE ANALYSIS\n N =", ntrials, "\n Mean confidence =", round(mean_conf, decimals), "\n SD confidence =", round(sd_conf, decimals), "\n")
  cat("\nCORRELATION ANALYSIS\n Phi =", round(phi, decimals), "\n p =", round(p_cor, decimals), "\n")
  cat("\nPROBABILITIES ANALYSIS\n O =", round(outcome, decimals), "\n C =", round(calibration, decimals), "\n R =", round(resolution, decimals), "\n PS =", round(PS, decimals), "\n")
  
  # FIRST PLOT
  
  plot.new()
  plot.window(xlim  = c(0.5, 1), ylim = c(0.5, 1), xaxs = "i", yaxs = "i")
  grid(5, 5, col = "lightgrey", lty = 2, lwd = 0.5)
  
  for (i in seq(0, 0.5, by = 0.001)) {
    segments(x0 = 0.5, y0 = (i + 0.5), x1 = (1 - i), y1 = 1, col = rgb(160/255, 32/255, 240/255, (i * 2)))
    segments(x0 = (i + 0.5), y0 = 0.5, x1 = 1, y1 = (1 - i), col = rgb(255/255, 165/255, 0/255, (i * 2)))}
  
  lines(x = subjective_probability, y = proportion_correct, col = "black", lty = 1, lwd = 5)
  lines(x = c(0.5, 1), y = c(0.5, 1), col = "black", lty = 1, lwd = 1)
  
  text(x = 0.725, y = 0.75, "No bias", srt = 45, col = "black")
  text(x = 0.6, y = 0.9, "Underconfidence", srt = 45, col = "purple")
  text(x = 0.9, y = 0.6, "Overconfidence", srt = 45, col = "orange")
  
  axis(1, xaxp = c(0.5, 1, 5))
  axis(2, yaxp = c(0.5, 1, 5))
  title(xlab = "Subjective probability", ylab = "Proportion correct", main = "Type I accuracy according to type II score")
  box(lwd = 2)
  
  # Wait before displaying the next plot
  cat ("\nPress [ENTER] to pursue the analysis.")
  line <- readline()
  
  # SECOND PART OF THE RESULTS
  
  cat("FIRST ORDER SDT\n d' =", round(d_prime, decimals), "\n c =", round(c, decimals), "\n beta0 = ", round(as.numeric(fit_typeI_poly$coefficients[2]), decimals), "\n beta1 =", round(as.numeric(fit_typeI_poly$coefficients[1]), decimals), "\n r2 =", round(as.numeric(fit_typeI_cor[decimals])^2, decimals), "\n p =", round(as.numeric(fit_typeI_cor[3]), decimals), "\n")
  cat("\nSECOND ORDER SDT\n meta-d' =", round(meta_d_prime, decimals), "\n D' =", round(big_d_prime, decimals), "\n beta0 =", round(fit_typeII_poly$coefficients[2], decimals), "\n beta1 =", round(as.numeric(fit_typeII_poly$coefficients[1]), decimals), "\n r2 =", round(as.numeric(fit_typeII_cor[decimals])^2, decimals), "\n p =", round(as.numeric(fit_typeII_cor[3]), decimals), "\n")
  cat("\nROC ANALYSIS\n Bk =", round(Bk, decimals), "\n Aroc =", round(Aroc, decimals), "\n Error =", round(error, decimals), "\n")

  # SECOND PLOT
  
  plot.new()
  plot.window(xlim  = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  grid(10, 10, col = "lightgrey", lty = 2, lwd = 0.5)
  
  # Plot the type I ROC curve
  polygon(x = p_false_alarm, y = p_hit, col = "grey")
  lines(x = p_false_alarm, y = p_hit, col = "black", lty = 1, lwd = 2)
  
  # Plot the type II ROC curve
  polygon(x = pFA_cum, y = pH_cum, col = "royalblue3")
  lines(x = pFA_cum, y = pH_cum, col = "black", lty = 1, lwd = 2)
  
  # Plot the diagonals
  lines(x = c(0, 1), y = c(0, 1), col = "black", lty = 1, lwd = 1)
  lines(x = c(0, 0.5), y = c(1, 0.5), col = "black", lty = 2, lwd = 1)
  
  # Define some other graph parameters
  axis(1, xaxp = c(0, 1, 10))
  axis(2, yaxp = c(0, 1, 10))
  title(xlab = "P (confidence | incorrect) or FA rate", ylab = "P (confidence | correct) or H rate", main = "ROC curves")
  legend("bottomright", c("Type I ROC curve", "Type II ROC curve"), col = c("grey", "royalblue3"), pch = c(15, 15), bg = "white")
  box(lwd = 2)
  
  # Return the results table
  return(summary_table)
}