import numpy as np
import matplotlib.pyplot as plt
import bolton as bt
import readsoundings as rd
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import FuncFormatter, Formatter
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

C_to_K = 273.15
skew_slope = 40.0
path = '/Users/ahardin/Documents/Graduate_School/Cloud_Physics/skew_T/FWD_20141014_00.txt'

# Functions
def x_from_Tp(T, p):
	""" Calculate x-coordinate using temperature (K) and pressure (mb). """
	x = T - skew_slope*np.log(p)
	return x

def y_from_p(p):
	""" Calculate y-coordinate only using pressure"""
	y = -np.log(p)
	return y

def T_from_xp(x,p):
	""" Calculate temperature using the x-coordinate and pressure. """
	temp = x + skew_slope*np.log(p)
	return temp

def p_from_y(y):
	""" Calculate pressure from y-coordinate. """
	p = np.exp(-y)
	return p

def to_thermo(x,y):
	""" Transform (x,y) coordinates to T in degrees Celsius and
	p in mb. """
	p = p_from_y(y)
	T_C = T_from_xp(x,p) - C_to_K
	return T_C, p

def from_thermo(T_C, p):
	""" Transform T_C (in degrees Celsius) and p (in mb) to (x,y). """
	y = y_from_p(p)
	x = x_from_Tp(T_C+C_to_K, p)
	return x, y
	
# Values along the bottom and left edges
p_bottom = 1050.0
p_top = 150.0
T_min = -40.0 + C_to_K
T_max = 50.0 + C_to_K

x_min = x_from_Tp(T_min, p_bottom)
x_max = x_from_Tp(T_max, p_bottom)
y_min = y_from_p(p_bottom)
y_max = y_from_p(p_top)

P_levels = np.arange(1000, 150-50, -50)       # mb
T_C_Levels = np.arange(-80, 40+10, 10)        # degrees C
T_Levels = T_C_Levels + C_to_K          # Convert temps to Kelvin
theta_levels = np.arange(-40,100+10,10) + C_to_K       # kelvin
theta_ep_levels = theta_levels.copy()
mixes = [.4,1,2,3,5,8,12,16,20]
mixing_ratios = np.asarray(mixes) /1000.0        # kg/kg

p_all = np.arange(p_bottom, p_top-1, -1)
y_p_levels = y_from_p(P_levels)
y_all_p = y_from_p(p_all)
x_T_levels = [x_from_Tp(Ti, p_all) for Ti in T_Levels]

x_thetas = [x_from_Tp(bt.theta_dry(theta_i, p_all), p_all) for theta_i in theta_levels]

x_mixing_ratios = [x_from_Tp((bt.mixing_ratio_line(p_all, w_i))+C_to_K, p_all) for w_i in mixing_ratios]

mesh_T, mesh_p = np.meshgrid(np.arange(-60.0, T_Levels.max()-C_to_K+.1, .1), p_all)
theta_ep_mesh = bt.theta_ep_field(mesh_T, mesh_p)

# Read in sounding data
sounding_data = rd.parse_SPC(path, skip_rows=7)    # changed skip_rows because Td and T had bad values in first row
snd_T = sounding_data['T']               # read in temp               
snd_Td = sounding_data['Td']             # read in dew point
snd_P = sounding_data['p']               # read in pressure

# all temperature values, deg. C, should be in this range
good_T = (snd_T > -100.0) & (snd_T < 60.0)
good_Td = (snd_Td > -100.0) & (snd_Td < 60.0)
good_P = (snd_P > 50.0) & (snd_P < 1005.0)

# Convert T, Td, and P to x-y coordinates
x_snd_T = x_from_Tp(snd_T+C_to_K, snd_P)
x_snd_Td = x_from_Tp(snd_Td+C_to_K, snd_P)
y_snd_p = y_from_p(snd_P)

# Plotting
skew_grid_helper = GridHelperCurveLinear((from_thermo, to_thermo))
fig = plt.figure()
ax = Subplot(fig, 1, 1, 1, grid_helper=skew_grid_helper)
fig.add_subplot(ax)

for yi in y_p_levels:
	ax.plot((x_min, x_max), (yi,yi), color=(1.0, .8, .8))
	
for x_T in x_T_levels:
	ax.plot(x_T, y_all_p, color=(1.0, .5, .6))

for x_theta in x_thetas:
	ax.plot(x_theta, y_all_p, color=(1.0, .7, .7))
	
for x_mixing_ratio in x_mixing_ratios:
	good = p_all >=600   # restrict mixing ratio lines to below 600mb
	ax.plot(x_mixing_ratio[good], y_all_p[good], color=(.8,.8,.6))

n_moist = len(theta_ep_levels)
moist_colors = ((.6,.9,.7),)*n_moist
ax.contour(x_from_Tp(mesh_T + C_to_K, mesh_p), y_from_p(mesh_p), theta_ep_mesh, theta_ep_levels, colors=moist_colors)

def format_coord(x, y):
	T, p = to_thermo(x, y)
	return "{0:5.1f C, {1:5.1f} mb".format(float(T), float(p))
ax.format_coord = format_coord
ax.plot(x_snd_Td, y_snd_p, linewidth=2, color='g')
ax.plot(x_snd_T, y_snd_p, linewidth=2, color='r')
ax.axis((x_min, x_max, y_min, y_max))
plt.title('Sounding FWD 10/14/2014 00Z')
plt.xlabel('Temperature (C)')
plt.ylabel('Pressure (mb)')
plt.show()
#plt.savefig('/Users/ahardin/Documents/Graduate_School/Cloud_Physics/skew_T/FWD_20141014_00.png')

##########################################################################################
# Number 8
'''One of the more confusing parts was number 6 part e. It took me a little while to 
figure out what to add into the code and where to put it. I would like to know how to 
calculate CAPE and CIN from these soundings and plot those two as well on the skew T. But
I did like that we were forced to use functions, which I'll be using more in the future. 
'''
	
	
	
	
	