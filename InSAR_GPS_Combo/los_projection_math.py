# A few functions that help project into and out of LOS
# Mostly trig
# Example: 
# Dlos = [U_n sin(phi) - U_e cos(phi)]*sin(lamda) + U_u cos(lamda)
# [U_e, U_n, U_u] are the east, north, and up components of the deformation. 
# phi   = azimuth of satellite heading vector, positive clockwise from north.
# lamda = local incidence angle at the reflector.
# from Fialko et al., 2001. 

import numpy as np 


def project_ENU_to_LOS(U_e,U_n,U_u,flight_angle,incidence_angle):
	# Dlos = [U_n sin(phi) - U_e cos(phi)]*sin(lamda) + U_u cos(lamda)
	# [U_e, U_n, U_u] are the east, north, and up components of the deformation. 
	# phi   = azimuth of satellite heading vector, positive clockwise from north.
	# lamda = local incidence angle at the reflector (usually angle from the vertical).
	# Fialko 2001
	# Works for single values and for 1D arrays
	
	phi=np.deg2rad(flight_angle); 
	lamda=np.deg2rad(incidence_angle);

	if np.size(U_e)>1:  # processing a 1D array of values
		d_los = [0*i for i in U_e];
		for i in range(len(U_e)):
			d_los[i] = ( (U_n[i]*np.sin(phi) - U_e[i]*np.cos(phi) )*np.sin(lamda)  + U_u[i]*np.cos(lamda) );
	else:               # processing a single value
		d_los = (U_n*np.sin(phi) - U_e*np.cos(phi) )*np.sin(lamda)  + U_u*np.cos(lamda) ;

	# Flipping by negative one, because the convention is positive when moving away from the satellite. 
	d_los=np.array(d_los)*-1;

	return d_los
