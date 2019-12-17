# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 16:59:07 2019

@author: s4524462
"""
'''
liquid water sauturation as a function of matric potential based on Fayer(1955)
Fayer, M. J., & Simmons, C. S. (1995). Modified Soil Water 
Retention Functions for All Matric Suctions. Water Resources 
Research, 31(5), 1233-1238.  
      Input:
      psim (m) -- matric potential in meters. this values needs to be negative
      av   (m) -- capillary fringe in matric potential meters 
                  This value must be negative.
      nv   (m) -- pore size distribution coefficient
      psim0(m) -- the matric potential value where liquid water saturation becomes zero
                  this value needs to be negative as well

'''
import numpy as np
def SWCC_Fayer1995WRR(psim,av,nv,psi0,saturation_residual):
    beta_Fayer = (np.log(-psi0)-np.log(-psim))/np.log(-psi0)
    saturation_effective = ( 1+(-psim*av)**nv   )**(1/nv-1)
    saturation           = (1-beta_Fayer*saturation_residual)*saturation_effective +beta_Fayer*saturation_residual
    return saturation_effective, saturation