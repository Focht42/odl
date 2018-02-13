import numpy as np # numpy is also used in the odl
# return the p and reg created in diri_IRGNM_CG.m

def diri_IRGNM_CG():
	p = dict()
	reg = dict()
	# Parameters of the forward problem
	p["op_name"] = "DirichletOP"
	p["syntheticdata_flag"] = True
	p["kappa"] = 3 # Wave Number
	# directions of incident waves
	N_inc = 1
	# original: t=2*pi*[0:N_inc-1]/N_inc
	# note: [a:b] is inclusive in matlab but exclusive in python
	t_range = np.array(list(range(N_inc)))
	t = 2*np.pi*t_range/N_inc
	p["inc_directions"] = np.vstack((np.cos(t),np.sin(t)))

	N_meas = 64
	# original: t=2*pi*[0:N_meas-1]/N_meas
	# note: [a:b] is inclusive in matlab but exclusive in python
	t_range = np.array(list(range(N_meas)))
	t = 2*np.pi*t_range/N_meas
	p["N_ieq"] = 64
	p["meas_directions"] = np.vstack((np.cos(t),np.sin(t)))
	p["plotWhat"] = dict()
	p["plotWhat"]["field"] = True
	p["plotWhat"]["ptsx"] = 60
	p["plotWhat"]["ptsy"] = 60
	p["plotWhat"]["plots"] = 'X'

	p["true_curve"] = 'nonsym_shape'
	p["noiselevel"] = 0.000
	p["N_FK"] = 64

	# parameters of the regularization method
	# regularization method

	reg["method"] = "IRGNM_CG"
	# parameter for abort of Newton iteration by discrepancy principle
	reg["tau"] = 2
	# maximum number of outer iterations
	reg["N_max_it"] = 20
	# starting reg parameter for IRGNM
	reg["alpha0"] = 0.001
	reg["plot_steps"] = np.array(list(range(0,reg["N_max_it"]+1)))
	reg["CG_TOL_rel_accY"] = 0.3
	reg["CG_TOL_rel_accX"] = 0.3
	reg["CG_TOL_resi_red"] = 0
	reg["verbose"] = 2

	return p, reg
