# 
# Saturated zones (zero sd) with an empty root zone (?) allow evaporation at the full potential rate
#
SaturatedEvap <- function(groups, flows, stores, pe)
{  
  # evap from saturated zone
  ae <- pe*stores$srz/groups$srz_max
  esat <- pe-ae # might be zero if root zone full
  storediff <- groups$sd_max - stores$sd

  evaporating <- which(esat>0 & storediff > 0)
  if(length(evaporating)>0)
  {
    browser()
  # update moisture deficit
    stores[evaporating,]$sd <- stores[evaporating,]$sd + storediff[evaporating]/esat[evaporating]
  }
  return(stores)
}