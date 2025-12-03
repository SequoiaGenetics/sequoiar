sg._cal_Fstats_binary = function(beta, eaf, n_case, n_control){
   var_X = n_case * n_control / (n_case + n_control) / (n_case + n_control)
   N = n_case + n_control
   R2 = 2*eaf*(1-eaf)*beta*beta/var_X
   K = length(eaf)
   Fstats = R2*(N-K-1)/K/(1-R2)
   return(Fstats)
}
