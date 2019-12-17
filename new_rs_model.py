# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:41:29 2019

From matlab to python
"""
# calculate surface resistance from soil surface

#variable
"""
evapo_rate_relative_capillary               % relative conductance of vapor through the soil air interface contributed to by the vaporization of the
                                              capillary water in the topmost soil layer
evapo_rate_relative_vapor                   % relative conductance of vapor through the soil air interface contributed to by vapor generated from the
                                              vaporization plane beneath the dry soil layer
R_0                   % radius of water saturated pores(bottom of funnel)
R_m                   % radius that corresponds to psi_m according to the Young_Laplace equation (m)
R_2                   % radius of water saturated pores(top of funnel)
R_3                   % radius of cylinder building block
R_c                   % characteristic radius  
Psi_m                 % matric potential
"""
#%% import parameters from json
import json
read = json.load(open('soil_parameter.json', 'r', encoding='utf-8'))

xi               = read['xi']                # Young_Laplace equation constant(assuming contact angle is zero)
psi_0_m          = read['psi_0_m']           # matric potential that corresponds to zero liquid water saturation
tortuosity_0     = read['tortuosity_0']      # tortuosity of liquid water in soils when the liquid water saturation is zero
diffusivity_m2Ps = read['diffusivity_m2Ps']  # diffusivity of vapor at 22 centigrade
free_path_gas_m  = read['free_path_gas_m']   # mean free path of gas molecules at 22 centigrade

#experimental conditions
thickness_NSL_m      = read['thickness_NSL_m']       # thickness of the near surface soil layer(NSL)
thickness_aero_edl_m = read['thickness_aero_edl_m']  # thickness of the external diffusive layer(EDL) by aerodynamics

#parameters about soil
psi_p_m             = read['psi_p_m']              # matric potential in the NSL corresponding to the initial liquid water saturation at early stage IV(m)
porosity            = read['porosity']            
saturation_residual = read['saturation_residual']  # residual liquid water saturation
n                   = read['n']                    # correction function between TSL and NSL
beta                = read['beta']                 # the characteristic angle of soil particle shape

#fitting parameter for the van Genuchten soil water rete3tion curve
#below are working parameters
av_Pm             = read['av_Pm'] 
nv                = read['nv']   
radius_particle_m = read['radius_particle_m']   # average particle size
#%% calculate
import numpy as np

psim_m_ay=np.hstack ((-np.arange(0.0001, 0.01, 0.0001), -np.arange(0.01, 0.1, 0.001),-np.arange(0.1, 1, 0.001),\
                      -np.arange(1, 10, 0.1), -np.arange(10, 100, 1), -np.arange(100, 1000, 10), -np.arange(1000, 50001, 500)))

r_m_ay   = xi/psim_m_ay
r_0_ay   = -xi*av_Pm*(-av_Pm*psim_m_ay)**(nv-1) * (  (1+ (-av_Pm*psim_m_ay)**-nv)  **(1-1/nv) -1)
r_c_m    = radius_particle_m/np.tan(beta);
r_2_ay   = r_0_ay+r_c_m;

saturation_effective_NSL_ay, saturation_NSL_ay=SWCC_Fayer1995WRR(psim_m_ay, av_Pm, nv, psi_0_m, saturation_residual)

thickness_funnel_ay            = radius_particle_m/np.exp(1/r_c_m*(r_m_ay-r_0_ay))
saturation_effective_TSL_ay    = saturation_effective_NSL_ay**(1+n)
relative_wetted_surface_ay     = saturation_effective_TSL_ay * (r_m_ay/(r_0_ay+r_c_m))**2
#relative_wetted_surface_ay     = saturation_effective_TSL_ay * (1-thickness_funnel_ay*(1-porosity)/radius_particle_m);

water_content_NSL_ay           = porosity * saturation_NSL_ay

r_3_ay                         = r_m_ay/relative_wetted_surface_ay**0.5;

evapo_rate_relative_capillary  = 1/(1  +  r_3_ay**2 * thickness_funnel_ay/r_2_ay**2/thickness_aero_edl_m + \
          r_m_ay/2/thickness_aero_edl_m/relative_wetted_surface_ay * \
         (  2*free_path_gas_m/r_m_ay  +  1/ (1+free_path_gas_m/r_m_ay ) - relative_wetted_surface_ay**0.5)  )

evapo_rate_relative_vapor      = (thickness_aero_edl_m+thickness_funnel_ay)*tortuosity_0*porosity*np.log(-psi_0_m) /\
                             thickness_NSL_m/np.log(-psim_m_ay-psi_p_m)
                
evapo_rate_relative            = evapo_rate_relative_capillary*(1-evapo_rate_relative_vapor)+evapo_rate_relative_vapor

#%% surface resistance
surface_resistance_new        = thickness_aero_edl_m/diffusivity_m2Ps*(1/evapo_rate_relative - 1)

surface_resistance_Griend_Owe    = 10*np.exp(35.63*(0.15-water_content_NSL_ay)) #surface resistance from model of van Genuchten

radius_pore_avg                  = xi*-av_Pm  #average pore size(from Lehmann-Soil Texture Effects on Surface Resistance to Bare-Soil Evaporation)
evapo_rate_relative_Sch          = 1/(1+2/np.pi*radius_pore_avg/(thickness_aero_edl_m/diffusivity_m2Ps*diffusivity_m2Ps)*(np.pi/4/water_content_NSL_ay)**0.5*((np.pi/4/water_content_NSL_ay)**0.5-1))
surface_resistance_Schlunder     = thickness_aero_edl_m/diffusivity_m2Ps*(1/evapo_rate_relative_Sch - 1) #surface resistance from model of van Genuchten

#%% plot
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5, 3))
#ax.plot(saturation_NSL_ay,surface_resistance_ay_sPm, color='b', marker='.', linestyle='solid')
ax.plot(saturation_NSL_ay,surface_resistance_new, 'r', label='new model' )
ax.plot(saturation_NSL_ay,surface_resistance_Griend_Owe, 'b', label='Griend & Owe'  )
ax.plot(saturation_NSL_ay,surface_resistance_Schlunder, 'g', label='Schlunder'  )

plt.yscale('symlog')
#ax.set_title('Combined debt growth over time')
ax.legend(loc='upper right')
ax.set_xlabel('Liquid water saturation (-)')
ax.set_ylabel('Surface resistance (s/m)')
ax.set_xlim(xmin=0, xmax=1)
ax.set_ylim(ymin=0.5, ymax=10000)
fig.tight_layout()
plt.show()
plt.savefig('Rs_v_Se.svg', format='svg', dpi=300)