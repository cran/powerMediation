########################################
# sample size per group for comparing mean change in longitudinal data with 2 time point
ssLong=function(es, rho=0.5, alpha=0.05, power=0.8)
{
  za=qnorm(1-alpha/2)
  zb=qnorm(power)
  n=4*(1-rho)*(za+zb)^2/es^2

  n.int=ceiling(n)
  return(n.int)
}

# delta = abs(mu1 - mu2)
#  m1 = underlying mean change over time t in group 1
#  m2 = underlying mean change over time t in group 2
# sigma1^2 = variance of baseline values within a treatment group 
# sigma2^2 = variance of follow-up values within a treatment group 
# rho = correlation coefficient between baseline and follow-up values
#   within a treatment group

# full parametrization
ssLongFull=function(delta, sigma1, sigma2, rho=0.5, alpha=0.05, power=0.8)
{
  sigmad2=sigma1^2+sigma2^2-2*rho*sigma1*sigma2
  za=qnorm(1-alpha/2)
  zb=qnorm(power)
  n=2*sigmad2*(za+zb)^2/delta^2

  n.int=ceiling(n)
  return(n.int)
}


# power for comparing mean change in longitudinal data with 2 time point
powerLong=function(es, n, rho=0.5, alpha=0.05)
{
  za=qnorm(1-alpha/2)
  power = pnorm(-za + abs(es)*sqrt(n/(1-rho))/2)

  return(power)
}

# full parametrization
powerLongFull=function(delta, sigma1, sigma2, n, rho=0.5, alpha=0.05)
{
  sigmad2=sigma1^2+sigma2^2-2*rho*sigma1*sigma2
  za=qnorm(1-alpha/2)

  power = pnorm(-za + delta*sqrt(n)/sqrt(2*sigmad2))

  return(power)
}

#ssLongFull(delta=5, sigma1=15, sigma2=15, rho=0.7, alpha=0.05, power=0.8)
#powerLongFull(delta=5, sigma1=15, sigma2=15, n=85, rho=0.7, alpha=0.05)

# Example 8.33 on page 336 of Rosner (2006)
# n=85
#ssLong(es=5/sqrt(225), rho=0.7, alpha=0.05, power=0.8)
# Example 8.33 on page 336 of Rosner (2006)
# power = 0.8
#powerLong(es=5/sqrt(225), n=85, rho=0.7, alpha=0.05)

# Example 8.34 on page 336 of Rosner (2006)
# power=0.75
#powerLong(es=5/sqrt(225), n=75, rho=0.7, alpha=0.05)


