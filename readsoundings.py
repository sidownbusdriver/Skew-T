import numpy as np

def parse_SPC(filename, skip_rows=6):
	""" Read in sounding data from the SPC """
	dtype = [ ('p', float),               # pressure, mb
	          ('z', float),               # altitude, m
	          ('T', float),               # temperature, C
	          ('Td', float),              # dewpoint, C
	          ('wind_dir', float),        # wind direction, degrees
	          ('wind_spd', float)         # wind speed, knots
	          ]
	data = np.genfromtxt(filename, dtype=dtype, skip_header=skip_rows, delimiter=',')
	return data