# INDIVIDUAL SIGNAL DETECTION ANALYSIS FOR TYPE II DATA (Maxime Maheu)

# Only requires an "input_data" data.frame with answers labels in the first column and confidence ratings in the second one; and confidence borns !
# It can also take into account design (detecion or choice), ponderation type, and homogeneity analysis

# SDT classifications (i.e. answers labels):
  # 1: Hit
  # 2: False alarm
  # 3: Miss
  # 4: Correct rejection

# Usefull references to undestand what is going on here :
  # Fleming, S. M., Weil, R. S., Nagy, Z., Dolan, R. J., & Rees, G. (2010). Relating introspective accuracy to individual differences in brain structure. Science, 329(5998), 1541–1543.
  # Kornbrot, D. E. (2006). Signal detection theory, the approach of choice: Model-based and distribution-free measures and evaluation. Perception & Psychophysics, 68(3), 393–414.
  # Macmillan, N. A., & Creelman, C. D. (2004). Detection theory: A user's guide. Psychology press.
  # Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430.
  # Szczepanowski, R., & Pessoa, L. (2007). Fear perception: Can objective and subjective awareness measures be dissociated? Journal of Vision, 7(4), 1–17.

MM_metaSDT <- function(input_data, min_conf, max_conf, design = 2, coefficient = 1) {

  # Print 
  cat("\nType I and type II signal detection analysis\n")
  cat("Copyright Maxime Maheu (2014)\n\n")
  
  # Check if the data input is a data frame
  if (is.data.frame(input_data) == FALSE) {cat("Warning! The input table is not a data frame. Please convert it.\n")}
  
  # Delete NA trials
  input_data <- na.omit(input_data)
  
  # Display some warnings
  ntrials <- nrow(input_data)
  mean_conf <- sum(example_data[2]) / ntrials
  if (ntrials < 100) {
    cat("Warning! The SDT analysis will only be made on less than 100 trials which can led to poor estimation of subject's actual metacognitive accuracy and, in some cases, to completely wrong computations.\n")}
  
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
      
      #
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
      if (coefficient == 2) { # À MODIFIER !!!!!!
        Aroc <- (1/2) + ((51/(2*101))*ka) + ((52/(2*101+1))*kb)
        Bk <- log(((51/(2*101)) * ka) / ((51/(2*101)) * kb))}
      
      # Save splitted Aroc values
      if (procedure == 1) {first_half <- Aroc}
      if (procedure == 2) {second_half <- Aroc}}
  
  # Display the results tables
  big_d_prime <- meta_d_prime / d_prime
  error <- first_half - second_half
  summary_table <- matrix(c(ntrials, mean_conf, d_prime, c, meta_d_prime, big_d_prime, Bk, Aroc, error, fit_typeI_poly$coefficients[2], fit_typeI_poly$coefficients[1], fit_typeI_cor[4], fit_typeI_cor[3], fit_typeII_poly$coefficients[2], fit_typeII_poly$coefficients[1], fit_typeII_cor[4], fit_typeII_cor[3]), nrow = 1, ncol = 17)
  colnames(summary_table) <- c("N", "Mean conf.", "d'", "c", "meta-d'", "D'", "Bk", "Aroc", "Error", "T1 beta1", "T1 beta2", "T1 cor", "T1 p", "T2 beta1", "T2 beta2", "T2 cor", "T2 p")
  rownames(summary_table) <- ""
  
  # Create the graph window
  plot.new()
  plot.window(xlim  = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  grid(10, 10, col = "lightgrey", lty = 2, lwd = 0.5)
  
  # Plot the type I ROC curve
  polygon(x = p_false_alarm, y = p_hit, col = "grey")
  lines(x = p_false_alarm, y = p_hit, col = "black", lty = 1, lwd = 1)
  
  # Plot the type II ROC curve
  polygon(x = pFA_cum, y = pH_cum, col = "royalblue3")
  lines(x = pFA_cum, y = pH_cum, col = "black", lty = 1, lwd = 1)
  
  # Plot the diagonal
  lines(x = c(0,1), y = c(0,1), col = "black", lty = 1, lwd = 1)
  
  # Define some other graph parameters
  axis(1, xaxp = c(0, 1, 10))
  axis(2, yaxp = c(0, 1, 10))
  title(xlab = "P (confidence | incorrect) or FA rate", ylab = "P (confidence | correct) or H rate", main = "ROC curves")
  legend("bottomright", c("Type I ROC curve", "Type II ROC curve"), col = c("grey", "royalblue3"), pch = c(15, 15), bg = "white")
  box(lwd = 2)
  
  # Return the results table
  return(print(summary_table, digits = 4))
}