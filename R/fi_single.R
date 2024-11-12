#' Fragility Index in Single Arm Design with Time to Event Endpoints
#'
#' @param data a data matrix with at least first two columns -- time and 
#'             censor indicator takes two values -- 1 if an event and
#'             0 otherwise. The name of the columns must be "time" and
#'             "status" respectively
#' @param median_benchmark the desired median survival beyond which the outcome 
#'                         is desired.
#' @param target_probability the probability threshold which is desired. 
#'                           generally very high value, default is 0.7.
#'
#' @details In an Oncology single arm trial often the decision is to invest 
#'          more on the compound if posterior probability that the median
#'          survival (or progression free survival) is more than median_probability
#'          is higher than the target probability. Here the posterior probability
#'          of the median survival is computed assuming the data follows an
#'          Exponential distribution. 
#'          
#' @references Maity, A. K. and Basu, C. (2024) Fragility Index in Single Arm Time to Event Trial
#' 
#' @return 
#' \item{fi}{Fragility Index}
#' \item{out}{An out put with time, indicator that the censored individual had been 
#'            converted to an event, and the corresponding posterior probability} 
#' 
#' 
#' @examples 
#' require(survival)
#' data(lung)
#' n <- 30 # sample size for hypothtical expansion cohort
#' set.seed(100) # fix the seed to sample from the full data
#' sample_id <- sample(x = 1:nrow(lung), size = n, replace = FALSE, prob = NULL)
#' lung_sample <- lung[sample_id, ]
#' lung_sample$timem <- lung_sample$time/30. # time in months
#' lung_sample$time <- lung_sample$timem
#' lung_sample$status <- lung_sample$status - 1 # the censored indicator is 1/2
#' 
#' # Kaplan Meier Plot
#' require(tidyverse)
#' ggsurvfit::survfit2(Surv(timem, status) ~ 1, data = lung_sample) %>% 
#' ggsurvfit::ggsurvfit() +
#' labs(
#' x = "Months",
#' y = "Overall survival probability"
#' ) + 
#' ggsurvfit::add_confidence_interval() +
#' ggsurvfit::add_risktable()
#' 
#' survfit(Surv(timem, status) ~ 1, data = lung_sample) 
#' # median of data
#' 
#' fi.single(data = lung_sample, median_benchmark = 7)
#' # Fragility index is 5 because when the 6th censored observation
#' # is switched to an event, the posterior probability becomes lower
#' # than the thersold 0.7.
#' 
#' @export


fi.single <- function(data, median_benchmark, target_probability = 0.7)
  # data will have
{
  n <- nrow(data) # sample size for expansion cohort
  
  if(ncol(data) == 2)
  {
    data$temp <- rep(99, n)
  }
  
  # prior parameters
  alpha <- 0.5
  beta <- 0.5
  
  nevent <- length(which(data$status == 1))
  tevent <- sum(data$time)
  nevent
  tevent
  # Posterior probability that the median survival is more than 8 months
  post_prob <- invgamma::pinvgamma(q = median_benchmark, shape = alpha + nevent, 
                                   rate = log (2) * (beta + tevent), lower.tail = FALSE)
  post_prob
  
  
  data_order <- data[order(data$status, data$time),] 
  # dataset ordered by censoring status; 0 = censored, 1 = dead
  
  i <- 0
  output <- NULL
  
  while((post_prob > target_probability) && (i < n))
    # the loop will stop as soon as the posyterior probability becomes under 0.7
    # 0.7 is a desirable high probability but it can be updated when needed
  {
    i <- i + 1
    data_order$status[i] <- 1 # replace censor by event
    
    nevent <- length(which(data_order$status == 1))
    tevent <- sum(data_order$time)
    
    post_prob <- invgamma::pinvgamma(q = median_benchmark, shape = alpha + nevent, 
                                     rate = log (2) * (beta + tevent), lower.tail = FALSE)
    
    out <- c(data_order$time[i], data_order$status[i], post_prob)
    output <- rbind(output, out)
    
  }
  
  output
  fi <- nrow(output) - 1 # Fragility Index
  fi
  
  return(list(fi = fi, out = output))
}