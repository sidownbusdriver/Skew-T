import numpy as np
import matplotlib.pyplot as plt


# Constants 
C_to_K = 273.15
c_p_dry = 1005.7          # J/kg K
R_d = 287.04              # J/kg K
c_V_dry = c_p_dry - R_d   # J/kg K
eps = 0.622
k_dry = 0.2854

# Functions
def sat_vapor_pressure(T):
	""" Calculates saturation vapor pressure (mb) using temperature in degrees Celsius"""
	e_s = 6.112*np.exp((17.67*T)/(T+243.5))
	return e_s
	
def sat_vapor_temperature(e_s):
	""" Calculates temperature (C) using only sat. vapor pressure (mb) """
	T = ((243.5*np.log(e_s)) - 440.8)/(19.48 - np.log(e_s))
	return T

def sat_mixing_ratio(p, T):
	""" Calculate saturation mixing ratio (kg/kg) using temp and pressure 
	pressure is in mb and temperature is in degrees celsius """
	e_s =  sat_vapor_pressure(T)
	ws = (eps*e_s)/(p - e_s)
	return ws
	
def mixing_ratio_line(p, w_s):
	""" Calculates temperature in degrees Celsius using p (mb) and w_s (kg/kg) """
	e_s = (w_s*p)/(eps + w_s)
	T = sat_vapor_temperature(e_s)
	return T

def RH(T, p , w):
	""" Calculates relative humidity using temperature (C), pressure (mb), and mixing 
	ratio (kg/kg) """
	w_s = sat_mixing_ratio(p, T)
	rh = (w/w_s) * 100
	return rh
	
def T_LCL(T, rh):
	""" Calculates temperature (K) of the LCL using temperature (C) and RH (%) """
	T_K = T + C_to_K
	lcl_temp = (1/((1/(T_K-55.0)) - (np.log(rh/100.0)/2840.0))) + 55.0
	return lcl_temp

def theta_dry(theta, p, p_0=1000.0):
	""" Calculates theta (K) when air is dry using theta/temp (C) and pressure (mb) """
	dry = theta/((p_0/p)**k_dry)
	return dry
	
def pseudoeq_potential_T(T, p, w, p_0=1000.0):
	""" Calculate pseudoadiabatic equivalent potential temperature (K) using temperature
	(C), pressure (mb), and mixing ratio (kg/kg) """
	temp = T + C_to_K
	rh = RH(T, p, w)
	Tlcl = T_LCL(T, rh)
	w = w*1000.0
	part1 = temp*((p_0/p)**(.2854*(1-(.28*10**-3)*w)))
	exponent = np.exp(((3.376/Tlcl)-.00254)*w*(1+(.81*10**-3)*w))
	thetaEP = part1 * exponent
	return thetaEP

def theta_ep_field(T, p, p_0 = 1000.0):
	""" Calculates theta_ep (K) using temperature (C), and pressure (mb) """
	w_s = sat_mixing_ratio(p, T)
	theta_ep = pseudoeq_potential_T(T, p, w_s)
	return theta_ep








