def inv_var_meta( studies : List[Tuple['Study','VariantData']], is_het_test = False):
  
  weights = []
  effs_size_org = []
  
  effs_inv_var = []
  sum_inv_var=0
  for s in studies:
    study = s[0]
  dat = s[1]
  if dat.se is None or dat.se==0:
    print("Standard error was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
  break
  var = (dat.se * dat.se)
  
  inv_var =  (1/var)
  sum_inv_var+=inv_var
  effs_inv_var.append( inv_var *  dat.beta )
  
  weights.append(inv_var)
  effs_size_org.append(dat.beta)
  
  beta_meta=sum(effs_inv_var)/ sum_inv_var
  if is_het_test:
    het_p=het_test(effs_size_org, weights, beta_meta)
  else:
    het_p=None
  return (beta_meta,
  math.sqrt(1/sum_inv_var),
  max(sys.float_info.min * sys.float_info.epsilon, 
  2 * scipy.stats.norm.sf(abs(sum(effs_inv_var) / math.sqrt(sum_inv_var) ))), het_p) if len(effs_inv_var)==len(studies) else None ')
  
  (beta_meta, standard_error, p_value, het_p)
  
  