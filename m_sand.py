# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:45:46 2019

@author: s4524462
"""
#parameters of medium sand

import json
import math

soil_parameter = {}

#constant
soil_parameter['xi']               =-1.469e-5  # Young_Laplace equation constant(assuming contact angle is zero)
soil_parameter['psi_0_m']          =-5e4       # matric potential that corresponds to zero liquid water saturation
soil_parameter['tortuosity_0']     = 0.66      # tortuosity of liquid water in soils when the liquid water saturation is zero
soil_parameter['diffusivity_m2Ps'] = 2.62e-5   # diffusivity of vapor at 22 centigrade
soil_parameter['free_path_gas_m']  = 0.6e-7    # mean free path of gas molecules at 22 centigrade

#experimental conditions
soil_parameter['thickness_NSL_m']      = 0.05  # thickness of the near surface soil layer(NSL)
soil_parameter['thickness_aero_edl_m'] = 5e-3  # thickness of the external diffusive layer(EDL) by aerodynamics

#parameters about soil
soil_parameter['psi_p_m']              = -10   # matric potential in the NSL corresponding to the initial liquid water saturation at early stage IV(m)
soil_parameter['porosity']             = 0.39
soil_parameter['saturation_residual']  = 0.09  # residual liquid water saturation
soil_parameter['n']                    = 0.0   # correction function between TSL and NSL
soil_parameter['beta']                 = math.pi/4 # the characteristic angle of soil particle shape

#fitting parameter for the van Genuchten soil water rete3tion curve
#below are working parameters
soil_parameter['av_Pm'] = 10.5
soil_parameter['nv']    = 5.5
soil_parameter['radius_particle_m']   = 3.5e-4 # average particle size

with open('soil_parameter.json', 'w') as json_file:
    json.dump(soil_parameter, json_file)

import new_rs_model


