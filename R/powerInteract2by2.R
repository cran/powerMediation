# two-way ANOVA model: 
#  y_{ijk}=mu + tau_i + beta_j + (tau*beta)_{ij} + epsilon_{ijk}
#  i=1, ..., a
#  j=1, ..., b
#  k=1, ..., n
#  epsilon_{ijk} i.i.d N(0, sigma^2)
#  sum_{i=1}^{a} tau_i=0
#  sum_{j=1}^{b} beta_j=0
#  sum_{i=1}^{a}sum_{j=1}^{b} (tau*beta)_{ij}=0
#
# a = 2: number of levels for factor 1
# b = 2: number of levels for factor 2
# n: number of subjects per group
# tauBetaSigma: (tau*beta)_{22}/sigma = theta/sigma
# sigma: standard deviation of random error
# alpha: family-wise type I error rate
# nTests is for high-throughput omics study, in which
#    we perform one-way ANOVA for each of 'nTests' probes.
#    We use Bonferroni correction to control for family-wise type I error rate
# verbose: logical flag indicating if intermediate results should be printed out

powerInteract2by2 = function(n, tauBetaSigma, alpha = 0.05, nTests = 1, verbose = FALSE)
{
  a = 2
  b = 2
  df1=(a-1)*(b-1)
  df2=a*b*(n-1)
  
  # control for family-wise type I error rate
  alpha2=alpha/nTests
  
  # obtain rejection region cutoff F0
  F0=qf(p=1-alpha2, df1=df1, df2=df2, ncp=0)
  
  # obtain delta
  
  delta = 4*tauBetaSigma^2
  
  # ncp under Ha
  ncp= n*delta
  
  power=1-pf(q=F0, df1=df1, df2=df2, ncp=ncp)
  
  if(verbose)
  {
    cat("df1=", df1, ", df2=", df2, ", F0=", F0, ", ncp=", ncp, 
        ", power=", power, "\n")
  }
  
  res=list(power=power, df1=df1, df2=df2, F0=F0, ncp=ncp)
  
  invisible(res)
}
