# Author: Mary Ann Binuya
# Last update: June 24, 2026

# References:
## Rentroia-Pacheco et al., 2024: https://github.com/emc-dermatology/ncc-evaluation/blob/main/Auxilliary%20Scripts/2023_03_03%20Weighted%20performance%20metrics%20functions.R
## Vickers and Elkin, 2006: https://pubmed.ncbi.nlm.nih.gov/17099194/

## Function to truncate until time "tp"
truncateFUP <- function(outcome_FUP, outcome_num, tp) {
  outcome_num_u <- outcome_num
  i.greaterFUP <- which(outcome_FUP > tp)
  if (length(i.greaterFUP) > 0) {
    outcome_num_u[outcome_FUP > tp] <- 0
  }
  outcome_FUP[outcome_FUP > tp] <- tp
  return(list("FUP" = outcome_FUP, "Status" = outcome_num_u))
}

## Function to calculate weighted C at tp. Note we choose C here as timeROC does not support weights.
computeCindex <- function(outcome_FUP, outcome_num, predictions, tp, weights = NULL) {
  # Truncate outcome at timepoint t:
  if (!is.null(tp)) {
    truncated_survival <- truncateFUP(outcome_FUP, outcome_num, tp)
    outcome_FUP <- truncated_survival[["FUP"]]
    outcome_num <- truncated_survival[["Status"]]
  }
  
  # Compute unweighted and weighted C-index:
  if (is.null(weights)) {
    cind <- 1 - rcorr.cens(predictions, Surv(outcome_FUP, outcome_num))[1]
  } else {
    cind <- intsurv::cIndex(outcome_FUP, event = outcome_num, predictions, weight = weights)[[1]]
  }
  names(cind) <- NULL
  return(cind)
}

## Function to compute C-related metrics
computeDiscMetrics <- function(outcome_FUP, outcome_num, predictions, cutoff, tp, weights = NULL) {
  
  # Obtain labels for the desired cutoff:
  labels_cutoff <- ifelse(predictions > cutoff, 1, 0)
  
  # Use equal weights if no weights are provided
  if (is.null(weights)) {
   weights <- rep(1, length(predictions))
  }
  
  # Compute survival probabilities for low-risk and high-risk groups:
  low_risk_surv <- summary(survfit(Surv(outcome_FUP[labels_cutoff == 0], outcome_num[labels_cutoff == 0]) ~ 1, weights = weights[labels_cutoff == 0]), times = tp, extend = TRUE)$surv
  high_risk_surv <- summary(survfit(Surv(outcome_FUP[labels_cutoff == 1], outcome_num[labels_cutoff == 1]) ~ 1, weights = weights[labels_cutoff == 1]), times = tp, extend = TRUE)$surv
  high_risk_n <- sum(ifelse(labels_cutoff == 1, weights, 0))
  low_risk_n <- sum(ifelse(labels_cutoff == 0, weights, 0))
  overall_surv <- summary(survfit(Surv(outcome_FUP, outcome_num) ~ 1, weights = weights), times = tp, extend = TRUE)$surv
  
  # Calculate True Negatives (TN), True Positives (TP), False Negatives (FN), and False Positives (FP):
  TN <- low_risk_surv * low_risk_n
  TP <- (1 - high_risk_surv) * high_risk_n
  FN <- (1 - low_risk_surv) * low_risk_n
  FP <- high_risk_surv * high_risk_n
  
  # Calculate overall Positive (Pos) and Negative (Neg) counts:
  Pos <- (1 - overall_surv) * sum(weights)
  Neg <- overall_surv * sum(weights)
  
  # Compute Sensitivity (SE):
  SE <- TP / Pos
  
  # Compute Specificity (SP):
  SP <- 1 - FP / Neg
  
  # Compute Positive Predictive Value (PPV):
  PPV <- TP / (TP + FP) #alternatively, PPV <- (1 - high_risk_surv)
  
  # Compute Negative Predictive Value (NPV):
  NPV <- TN / (TN + FN) #alternatively, NPV <- low_risk_surv
  
  # Compute Number Needed to Screen (NNS); NN to select to identify one TP event by tp
  Risk_unscreened <- (TP + FN) / sum(weights) # Weighted approximation of the baseline risk or prevalence.
  Risk_screened <- FN / sum(weights) # Weighted risk among those identified as low risk (missed cases).
  ARR <- Risk_unscreened - Risk_screened
  NNS <- 1 / ARR # Reduces to: sum(weights) / TP
  
  # Summarize all metrics:
  performance_metrics <- c("SE" = SE, "SP" = SP, "PPV" = PPV, "NPV" = NPV, "NNS" = NNS)
  
  return(performance_metrics)
}

## Function to calculate O/E ratio
computeOEratio <- function(outcome_FUP, outcome_num, predictions, tp, weights = NULL) {
  # Check if weights are provided
  if (is.null(weights)) {
    # Unweighted version
    obj <- summary(survfit(Surv(outcome_FUP, outcome_num) ~ 1), 
                   times = tp, extend = TRUE)
    OE_ratio <- (1 - obj$surv) / mean(predictions)
  } else {
    # Weighted version
    obj <- summary(survfit(Surv(outcome_FUP, outcome_num) ~ 1, weights = weights), 
                   times = tp, extend = TRUE)
    OE_ratio <- (1 - obj$surv) / weighted.mean(predictions, weights)
  }
  return(OE_ratio)
}

## Function to calculate calibration slope
computeSlope_loglog <- function(outcome_FUP, outcome_num, predictions, tp = NULL, weights = NULL) {
  # Truncate outcome at timepoint tp (if provided)
  if (!is.null(tp)) {
    truncated_survival <- truncateFUP(outcome_FUP, outcome_num, tp)
    outcome_FUP <- truncated_survival[["FUP"]]
    outcome_num <- truncated_survival[["Status"]]
  }
  
  # Compute log-log transformation of predictions
  log_log_pred <- log(-log(1 - predictions)) # Log-log of the survival predictions
  
  # Fit Cox proportional hazards model and extract calibration slope
  cal_slope <- coef(coxph(Surv(outcome_FUP, outcome_num) ~ log_log_pred, weights = weights))
  names(cal_slope) <- NULL
  
  return(cal_slope)
}
# For Cox-based absolute risk predictions: log(-log(1 - p_i(t))) = log(H0(t)) + LP_i.
# Fitting coxph(... ~ log(-log(1 - p_i(t)))) estimates a Cox calibration slope for the prognostic index.

## Alternative calibration slope using grouped observed risks (if absolute rather than log hazard scale is preferred)
computeSlope <- function(outcome_FUP, outcome_num, predictions, tp = NULL, weights = NULL, num_groups = 10) {
  # Truncate outcome at timepoint tp (if provided)
  if (!is.null(tp)) {
    truncated_survival <- truncateFUP(outcome_FUP, outcome_num, tp)
    outcome_FUP <- truncated_survival[["FUP"]]
    outcome_num <- truncated_survival[["Status"]]
  }
  
  # Create the data frame
  dat <- data.frame(outcome_FUP = outcome_FUP, outcome_num = outcome_num, predictions = predictions, weights = weights)
  
  # Add a small jitter to predictions if breaks are not unique
  jittered_predictions <- jitter(dat$predictions, factor = 1e-5)
  
  # Create groups (deciles by default)
  q <- cut(jittered_predictions, breaks = quantile(jittered_predictions, probs = seq(0, 1, length.out = num_groups + 1)), include.lowest = TRUE)
  dat$q_f <- factor(q, levels = levels(q), labels = paste0("q", 1:num_groups))
  
  # Fit a KM curve to estimate survival probabilities for each group
  obs <- survfit(Surv(outcome_FUP, outcome_num) ~ q_f, weights = weights, data = dat)
  obs_q <- summary(obs, times = tp, extend = TRUE)
  
  # Split data by q_f
  split_data <- split(dat, dat$q_f)
  
  # Calculate weighted means for each group
  weighted_means <- sapply(split_data, function(subset) {
    weighted.mean(subset$predictions, subset$weights)
  })
  
  # Create the calplot dataframe
  calplot <- data.frame(y = 1 - obs_q$surv,
                        x = as.numeric(weighted_means))
  
  # Fit a linear model to compute the calibration slope
  cal_model <- lm(y ~ x, data = calplot)
  cal_slope <- coef(cal_model)[2]
  
  return(cal_slope)
}

computeNRI <- function(outcome_FUP, outcome_num, predictions_old, predictions_new, tp, weights = NULL) {
  
  # Use equal weights if no weights are provided
  if (is.null(weights)) {
    weights <- rep(1, length(predictions_old))
  }
  
  # Define risk categories
  risk_cats <- c(0, 0.017, 0.03, Inf)
  labels <- c("Low", "Moderate", "High")
  
  # Assign risk categories for old and new predictions
  old_risk_cat <- cut(predictions_old, breaks = risk_cats, labels = labels, right = FALSE)
  new_risk_cat <- cut(predictions_new, breaks = risk_cats, labels = labels, right = FALSE)
  
  # Convert factor levels to numeric indices for comparisons
  old_risk_num <- as.numeric(old_risk_cat)
  new_risk_num <- as.numeric(new_risk_cat)
  
  # Define event and non-event status at tp using original follow-up
  events <- outcome_num == 1 & outcome_FUP <= tp
  non_events <- outcome_FUP >= tp & !events
  
  # Estimate censoring survival G(t) = P(not censored before t)
  censor_event <- as.numeric(outcome_num == 0)
  
  G_fit <- survival::survfit(
    survival::Surv(outcome_FUP, censor_event) ~ 1,
    weights = weights
  )
  
  get_G <- function(times) {
    if (length(times) == 0) return(numeric(0))
    pmax(summary(G_fit, times = times, extend = TRUE)$surv, 1e-6)
  }
  
  # IPCW weights
  ipcw <- rep(NA_real_, length(outcome_FUP))
  ipcw[events] <- 1 / get_G(outcome_FUP[events])
  ipcw[non_events] <- 1 / get_G(rep(tp, sum(non_events)))
  
  weights_ipcw <- weights * ipcw
  
  # Define upward and downward movements
  upward <- !is.na(old_risk_num) & !is.na(new_risk_num) & new_risk_num > old_risk_num
  downward <- !is.na(old_risk_num) & !is.na(new_risk_num) & new_risk_num < old_risk_num
  
  # Count upward and downward movements for events
  event_upward <- sum(weights_ipcw[events & upward], na.rm = TRUE)
  event_downward <- sum(weights_ipcw[events & downward], na.rm = TRUE)
  
  # Count upward and downward movements for non-events
  non_event_upward <- sum(weights_ipcw[non_events & upward], na.rm = TRUE)
  non_event_downward <- sum(weights_ipcw[non_events & downward], na.rm = TRUE)
  
  # Compute total number of events and non-events
  total_events <- sum(weights_ipcw[events], na.rm = TRUE)
  total_non_events <- sum(weights_ipcw[non_events], na.rm = TRUE)
  
  # Compute Net Reclassification Improvement components
  NRI_events <- ifelse(
    total_events > 0,
    (event_upward - event_downward) / total_events,
    NA_real_
  )
  
  NRI_non_events <- ifelse(
    total_non_events > 0,
    (non_event_downward - non_event_upward) / total_non_events,
    NA_real_
  )
  
  NRI <- NRI_events + NRI_non_events
  
  # Return results
  nri_results <- c(
    "NRI Events" = NRI_events,
    "NRI Non-Events" = NRI_non_events,
    "Overall NRI" = NRI
  )
  
  return(nri_results)
}

## Compute Net Benefit
computeNB <- function(outcome_FUP, outcome_num, predictions, weights, tp, thresholdmax){
  
  # Truncate outcome at timepoint t:
  if (!is.null(tp)){
    truncated_survival <- truncateFUP(outcome_FUP, outcome_num,tp)
    outcome_FUP <- truncated_survival[["FUP"]]
    outcome_num <- truncated_survival[["Status"]]
  }
  
  data <- data.frame(outcome_FUP, outcome_num, predictions, weights)
  
  # Apply NB to all thresholds until thresholdmax
  thresholds <- seq(0.01, thresholdmax, by = 0.001) 
  
  NB_f <- lapply(thresholds, function(pt) {
    
    # NB treat all
    m_all <- 1 - summary(survfit(Surv(outcome_FUP, outcome_num) ~ 1,
                         weights = weights, data = data),
                         times = tp, extend = TRUE)$surv
    NB_all <- m_all - (1-m_all) * (pt/(1-pt))
    
    # NB using predicted risk ("predictions")
      # Proportion of high risk subjects:
      prop_pred =  weighted.mean(ifelse(data$predictions >= pt, 1, 0), data$weights)
      
      surv <- try(
        summary(survfit(Surv(outcome_FUP, outcome_num) ~ 1,
                        weights = weights, data = data[data$predictions >= pt, ]),
                        times = tp, extend = TRUE), silent = TRUE)
      
      if (class(surv) == "try-error") {
        TP <- 0
        FP <- 0
        NB <- 0 #no observations above threshold
      } else {
        m_model <- 1 - surv$surv
        TP <- m_model * prop_pred
        FP <- (1 - m_model) * prop_pred
        NB <- TP - FP * (pt/(1-pt))
      }
    
      NBres <- data.frame("threshold" = pt,
                          "NB_all" = NB_all,
                          "NB" = NB)
    
  })
  
  NB_res <- do.call(rbind, NB_f)
  
  return(NB_res)
}

# Calibration plot function
calplot_smooth <- function(data, predicted_risk, tmax, main = '',
                       C, calslope, OEratio,
                       limit = 0.5, 
                       size_lab = 1, size_legend = 1,
                       markers = FALSE, g = 5) {
  df <- data
  
  # Predicted risks
  df$x <- df[, predicted_risk]
  
  df$x.ll <- log(-log(1 - df$x))
  
  model <- cph(Surv(futime, bc) ~ rcs(x.ll), data = df, x = TRUE, y = TRUE, weights = wt, surv = TRUE)
  
  # Observed risks
  xx <- seq(quantile(df$x, prob = 0.01), quantile(df$x, prob = 0.99), length = 100)
  xx.ll <- log(-log(1 - xx))
  xx.ll.df <- data.frame(x.ll = xx.ll)
  
  y <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$surv
  y.lower <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$upper
  y.upper <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$lower
  
  # Plot parameters
  xlim <- c(0, limit + 0.005)
  ylim <- c(0, limit + 0.005)
  xlab <- 'Predicted risk'
  ylab <- 'Observed proportion'
  
  # Plot
  plot(0, 0, type = 'n',
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab,
       main = main, cex.main = 1.2,
       cex.lab = size_lab, cex.axis = size_lab,
       las = 1)
  polygon(c(xx, rev(xx), xx[1]),
          c(y.lower, rev(y.upper), y.lower[1]),
          border = NA, density = 50, angle = -20, col = 'gray')
  abline(coef = c(0,1), lty = 2, col = 'gray') #diagonal line
  lines(xx, y, lwd = 2, col = 'black')
  
  # Markers
  if (markers) {
    q <- Hmisc::cut2(df$x, levels.mean = TRUE, g = g) # group predicted risks
    means <- as.double(levels(q))
    y1 <- 1 - survest(model, newdata = df, times = tmax)$surv
    prop <- tapply(y1, q, mean) # mean observed (fitted) risks
    points(means, prop, pch = 124, cex = 1, col = 'black') # add markers
  }
  
  # Texts
  limit_inc <- limit / 3 / 5
  text(x = 0, y = (limit), labels = paste('C-index: ', C, sep = ''), cex = size_legend, pos = 4)
  text(x = 0, y = (limit - limit_inc), labels = paste('O/E ratio: ', OEratio, sep = ''), cex = size_legend, pos = 4)
  text(x = 0, y = (limit - 2 * limit_inc), labels = paste('Slope: ', calslope, sep = ''), cex = size_legend, pos = 4)
}

calplot_quintiles <- function(predrisk_col, model_name, limit = 0.05) {
  dat$predrisk <- predrisk_col
  
  # Create groups using quintiles
  unique_breaks <- unique(quantile(dat$predrisk, probs = seq(0, 1, 0.2), na.rm = TRUE))
  q5 <- cut(dat$predrisk, breaks = unique_breaks, include.lowest = TRUE)
  dat$q5_f <- factor(q5, levels = levels(q5), labels = c("q1", "q2", "q3", "q4", "q5"))
  
  # Observed risk
  obs <- survfit(Surv(futime, bc) ~ q5_f, weights = wt, data = dat)
  obs_q5 <- summary(obs, times = tmax, extend = TRUE)
  
  # Split data by q5_f
  split_data <- split(dat, dat$q5_f)
  
  # Calculate weighted means for each group
  weighted_means <- sapply(split_data, function(subset) {
    weighted.mean(subset$predrisk, subset$wt)
  })
  
  # Create the calplot dataframe
  calplot <- data.frame(
    y = 1 - obs_q5$surv,
    y_l = 1 - obs_q5$upper,
    y_u = 1 - obs_q5$lower,
    x = as.numeric(weighted_means))
  
  # Plot
  par(las = 1, xaxs = "i", yaxs = "i")
  plot(calplot$x, calplot$y,
       type = "b", bty = "n", pch = 15, col = "black", lty = 1,
       xlim = c(0, limit), ylim = c(0, limit),
       xlab = paste("Predicted risk"),
       ylab = paste("Observed proportion"))
  
  # Add confidence intervals
  plotrix::plotCI(x = calplot$x, y = calplot$y,
                  li = calplot$y_l, ui = calplot$y_u,
                  xlim = c(0, limit), ylim = c(0, limit),
                  add = TRUE, col = "black", pch = 16)
  
  # Add a diagonal reference line
  abline(a = 0, b = 1, lwd = 1, lty = 2, col = "gray")
  
  # Set the title for the plot
  title(model_name, adj = 0.5)
  box()
}

calplot_deciles <- function(predrisk_col, model_name, limit = 0.05) {
  dat$predrisk <- predrisk_col
  
  # Create groups using deciles
  unique_breaks <- unique(quantile(dat$predrisk, probs = seq(0, 1, 0.1), na.rm = TRUE))
  q10 <- cut(dat$predrisk, breaks = unique_breaks, include.lowest = TRUE)
  dat$q10_f <- factor(q10, levels = levels(q10), labels = paste0("d", 1:length(levels(q10))))
  
  # Observed risk
  obs <- survfit(Surv(futime, bc) ~ q10_f, weights = wt, data = dat)
  obs_q10 <- summary(obs, times = tmax, extend = TRUE)
  
  # Split data by q10_f
  split_data <- split(dat, dat$q10_f)
  
  # Calculate weighted means for each group
  weighted_means <- sapply(split_data, function(subset) {
    weighted.mean(subset$predrisk, subset$wt)
  })
  
  # Create the calplot dataframe
  calplot <- data.frame(
    y = 1 - obs_q10$surv,
    y_l = 1 - obs_q10$upper,
    y_u = 1 - obs_q10$lower,
    x = as.numeric(weighted_means))
  
  # Plot
  par(las = 1, xaxs = "i", yaxs = "i")
  plot(calplot$x, calplot$y,
       type = "b", bty = "n", pch = 15, col = "black", lty = 1,
       xlim = c(0, limit), ylim = c(0, limit),
       xlab = paste("Predicted risk"),
       ylab = paste("Observed proportion"))
  
  # Add confidence intervals
  plotrix::plotCI(x = calplot$x, y = calplot$y,
                  li = calplot$y_l, ui = calplot$y_u,
                  xlim = c(0, limit), ylim = c(0, limit),
                  add = TRUE, col = "black", pch = 16)
  
  # Add a diagonal reference line
  abline(a = 0, b = 1, lwd = 1, lty = 2, col = "gray")
  
  # Set the title for the plot
  title(model_name, adj = 0.5)
  box()
}
