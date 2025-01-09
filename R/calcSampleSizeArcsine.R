#' Calculate Sample Size Using Arcsine Transformation
#'
#' This function calculates the required sample size for single-arm survival studies
#' based on the arcsine transformation method. It accounts for uniform accrual
#' and exponential survival assumptions, including numeric integration for time
#' points that exceed the follow-up period.
#'
#' @param S0 Numeric. Survival probability under the null hypothesis (must be strictly between 0 and 1).
#' @param S1 Numeric. Survival probability under the alternative hypothesis (must be strictly between 0 and 1).
#' @param alpha Numeric. The one-sided Type I error rate. Default is 0.05.
#' @param power Numeric. Desired statistical power of the test (1 - beta). Default is 0.80.
#' @param accrual Numeric. Duration of the accrual period in months. Default is 24.
#' @param followup Numeric. Additional follow-up duration in months after accrual. Default is 24.
#' @param timePoint Numeric. Time of interest in months for evaluating survival probabilities. Default is 18.
#' @param steps Integer. Number of steps for numeric integration if \code{timePoint} exceeds follow-up duration. Default is 10,000.
#'
#' @return Integer. The required sample size, rounded up to the nearest whole number.
#'
#' @examples
#' # Calculate sample size for typical survival probabilities
#' calcSampleSizeArcsine(S0 = 0.90, S1 = 0.96)
#'
#' # Adjusting for lower survival probabilities and extended accrual
#' calcSampleSizeArcsine(
#'   S0 = 0.80,
#'   S1 = 0.85,
#'   accrual = 36,
#'   followup = 12,
#'   timePoint = 24
#' )
#'
#' @importFrom stats qnorm uniroot
#' @export

###############################################################################
# REVISED FUNCTION FOR ARCSINE-BASED SAMPLE SIZE UNDER
# UNIFORM ACCRUAL + EXPONENTIAL SURVIVAL
#
# Now properly handles t > followup by numeric integration of the
# Greenwood variance (Section 2.4 of Nagashima et al. 2021).
###############################################################################

calcSampleSizeArcsine <- function(S0, S1,
                                  alpha     = 0.05,   # one-sided alpha
                                  power     = 0.80,   # desired 1 - beta
                                  accrual   = 24,     # accrual duration (months)
                                  followup  = 24,     # additional follow-up (months)
                                  timePoint = 18,      # time of interest (months)
                                  steps     = 10000    # for numeric integration if t>followup
) {


  #-------------------------------------------------------------------------
  # 0) Basic checks
  #-------------------------------------------------------------------------
  if(S0 <= 0 || S0 >= 1) stop("S0 must be strictly between 0 and 1.")
  if(S1 <= 0 || S1 >= 1) stop("S1 must be strictly between 0 and 1.")
  if(S1 <= S0) {
    warning("S1 <= S0: effect size is non-positive or zero. Returning NA.")
    return(NA)
  }
  if(timePoint <= 0) stop("timePoint must be > 0.")
  if(accrual < 0 || followup < 0) stop("accrual and followup must be >= 0.")
  if((accrual + followup) < timePoint) {
    message("WARNING: timePoint is beyond (accrual + followup). Interpretation requires caution.")
  }



  #-------------------------------------------------------------------------
  # 1) Convert S(t) -> hazard lambda. We assume exponential: S(t)= exp(-lambda*t)
  #-------------------------------------------------------------------------
  lambda0 <- -log(S0) / timePoint
  lambda1 <- -log(S1) / timePoint

  #-------------------------------------------------------------------------
  # 2) Function to compute Var{ KM(t) } under:
  #     - Exponential survival with rate 'lambda'
  #     - Uniform accrual in [0, accrual], plus followup
  #
  #   Implements numeric Greenwood integral if t > followup
  #   per Section 2.4 of Nagashima et al. (2021).
  #-------------------------------------------------------------------------
  kmVarExpoUniform <- function(lambda, t, n) {
    # total study length: 0..accrual for enrollment + followup for observation
    totalTime <- accrual + followup
    # Survival at t
    St <- exp(-lambda * t)

    #---------------------------
    # If t <= followup:
    #    => no forced censoring by time t. (Everyone who enrolled before t can
    #       be observed at t). Then the usual binomial-like approx is valid:
    #---------------------------
    if(t <= followup) {
      return( St*(1 - St) / n )
    }

    #---------------------------
    # If t > followup:
    #   => partial censoring. Use Greenwood integral for expo + uniform accrual
    #---------------------------
    # The standard continuous Greenwood formula can be approximated by:
    #
    #    Var{S(t)}  ~  S(t)^2 * ∫[0..t] [ dH(u) / S(u)^2 ],
    #
    # where dH(u)= hazard(u)*du (under exponential => lambda du),
    # but we must scale by the fraction at risk.  In practice:
    #
    #    Var{S(t)}  ~  S(t)^2 * ∫[0..t] [ (#events in du ) / (Y(u)*S(u))^2 ],
    #
    # For no other random censoring, # events in du ~ lambda * Y(u)*du,
    # so the integrand ~ lambda*du / [Y(u)*S(u)]^2 * S(t)^2, etc.
    #
    # A simpler route (common in textbooks) is:
    #    Var(KM(t)) ~ S(t)^2 * ∫_0^t [ lambda / Y(u) ] du
    #
    # where Y(u) = n * fractionAtRisk(u), fractionAtRisk(u) from uniform accrual
    # and the fraction still alive by time u.
    #
    # We'll do a numeric approximation from 0.. min(t, totalTime).
    # Also note that once we exceed totalTime, there's nobody at risk.
    # So effectively integrate from 0.. min(t, totalTime).
    #

    # fractionAtRisk(u): fraction of the cohort still in the study at time u
    #    = 1/accrual * ∫[x=0..min(u,accrual)] exp{ -lambda*(u - x) } dx
    # because those who entered at time x (0<=x<=accrual) must survive (u-x).
    #
    # We'll implement fractionAtRisk(u) piecewise:

    fractionAtRisk <- function(u) {
      if(u <= 0) return(0)
      # If u <= accrual:
      #   fraction = (1/accrual)* ∫[0..u] e^{-lambda(u - x)} dx
      # If u > accrual:
      #   fraction = (1/accrual)* ∫[0..accrual] e^{-lambda(u - x)} dx
      if(u <= accrual) {
        # integral( e^{-lambda(u-x)} dx, x=0..u ) = (1 - exp(-lambda*u)) / lambda * exp(-lambda*(u-u))?
        # More simply:
        #   ∫[0..u] e^{-lambda(u-x)} dx = e^{-lambda u} ∫[0..u] e^{lambda x} dx
        #   = e^{-lambda u} * [ (e^{lambda u} - 1)/lambda ]
        #   = (1/lambda)(1 - e^{-lambda u}).
        return( (1/accrual) * (1/lambda) * (1 - exp(-lambda*u)) )
      } else {
        # integral( e^{-lambda(u-x)} dx, x=0..accrual )
        #   = e^{-lambda u} * ∫[0..accrual] e^{lambda x} dx
        #   = e^{-lambda u} * [ (e^{lambda*accrual} - 1)/lambda ]
        return( (1/accrual) * (1/lambda) * exp(-lambda*u) *
                  (exp(lambda*accrual) - 1) )
      }
    }

    # numeric integration of integrand(u)= lambda / [ n * fractionAtRisk(u) ]
    # from 0.. min(t, totalTime)
    upperLim <- min(t, totalTime)
    nSteps   <- steps
    du       <- upperLim / nSteps
    accum    <- 0

    for(i in seq_len(nSteps)) {
      # midpoint in each sub-interval for a simple Riemann/rectangle
      uMid   <- (i - 0.5)*du
      if(uMid >= 0 && uMid <= upperLim) {
        fr <- fractionAtRisk(uMid)
        if(fr > 1e-15) {
          accum <- accum + (lambda / (n*fr))
        }
      }
    }
    # multiply by step
    greenwoodIntegral <- accum * du

    # Then Var(KM(t)) ~ S(t)^2 * greenwoodIntegral
    varKM <- (St^2)* greenwoodIntegral
    return(varKM)
  } # end function kmVarExpoUniform


  #-------------------------------------------------------------------------
  # 3) arcsine transform & derivatives
  #-------------------------------------------------------------------------
  gArcsine      <- function(S)  2 * asin( sqrt(S) )
  gprimeArcsine <- function(S)  1 / sqrt( S * (1 - S) )

  #-------------------------------------------------------------------------
  # 4) effect size
  #-------------------------------------------------------------------------
  epsilon <- gArcsine(S1) - gArcsine(S0)

  #-------------------------------------------------------------------------
  # 5) Solve for n by uniroot:
  #    We want: z_{1-alpha}*tau0(n) + z_{1-beta}*tau1(n) = epsilon
  #-------------------------------------------------------------------------
  z_alpha <- qnorm(1 - alpha)  # one-sided
  z_beta  <- qnorm(power)      # for 1 - beta

  f_lhs <- function(nTest) {
    # Var under null
    var0 <- kmVarExpoUniform(lambda0, timePoint, nTest)
    # Var under alt
    var1 <- kmVarExpoUniform(lambda1, timePoint, nTest)

    tau0 <- sqrt( gprimeArcsine(S0)^2 * var0 )
    tau1 <- sqrt( gprimeArcsine(S1)^2 * var1 )
    lhs  <- z_alpha * tau0 + z_beta * tau1
    return(lhs)
  }

  solveFun <- function(nTest) f_lhs(nTest) - epsilon

  # uniroot on [2, 1e7]
  out  <- uniroot(solveFun, interval = c(2, 1e7), tol = 1e-9)
  nSol <- out$root
  nReq <- ceiling(nSol)

  return(nReq)
}
