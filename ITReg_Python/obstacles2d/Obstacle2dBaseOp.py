import numpy as np

def Obstacle2dBaseOp_plot(F,x_k,x_start,y_k,y_obs,k):
    def fixme(F,x_k,x_start,y_k,y_obs,k):
        print("TODO: Implement this!")
    F = fixme
    #F = scat_plot(F,x_k,x_start,y_k,y_obs,k) # scat_plot is a matlab plotter
    return F

def Obstacle2dBaseOp_applyGramX(F,v):
    res = F["bd"]["applyGram"](F["bd"],v)
    return res

def Obstacle2dBaseOp_applyGramX_inv(F,v):
    res = F["bd"]["applyGram_inv"](F["bd"],v)
    return res

def Obstacle2dBaseOp_L2err(F,h):
    res = F["bd"]["L2err"](h,F["xdag"])
    return res

def Obstacle2dBaseOp_applyGramY(F,v):
    res = (2*np.pi/F["Ydim"])*v
    return res


def Obstacle2dBaseOp(p):

    pref["syntheticdata_flag"] = true
    pref["sobo_Index"] = 1.6 # Sobolev index of preimage space
    pref["init_guess"] = []
    pref["y_obs"] = []
    if "N_FK" in p.keys():
        pref["N_FK"] = p["N_FK"]
    else
        pref["N_FK"] = 64 #number of Fourier coefficients

    if 'bd_type' in p.keys():
        pref["bd_type"] = p["bd_type"]
    else
        pref["bd_type"] = 'StarTrig'

    if 'true_curve' in p.keys():
        pref["true_curve"] = p["true_curve"]
    else
        pref["true_curve"] = 'nonsym_shape'

    pref["noiselevel"] = 0.05

    ## compute precomputable quantities
    # initialize boundary curve structure

    if pref["syntheticdata_flag"]:
        F["bd_ex"] = GenCurve(pref["true_curve"])
        N_plot =128
        bd_ex = F["bd_ex"]["bd_eval"](F["bd_ex"],N_plot,0)
        F["xdag_pts_plot"] = bd_ex["z"]

    if pref["bd_type"] in ["StarTrig"]:
        F["bd"] = StarTrig(pref["N_FK"],pref["sobo_Index"])
        pref["init_guess"] = ones(pref["N_FK"],1)  #initial guess = unit circle
        F["Xdim"] = pref["N_FK"]
        if pref["syntheticdata_flag"]:
            F["xdag"] = F["bd_ex"]["radial"](F["bd_ex"],pref["N_FK"]).getH();

    elif pref["bd_type"] in ['GenTrig']
        F["bd"] = GenTrig(pref["N_FK"],pref["sobo_Index"])
        # remember python range differs from Matlab range!
        #TODO i suspect i will run into some matrix issues here!
        t_h = np.matrix(np.array(list(range(0,pref["N_FK"])))).getH() # i suspect the getH is wrong
        t = 2*pi*t_h/pref["N_FK"] #assuming here that N_FK is a single number
        ig_cos = np.cos(t)
        ig_sin = np.sin(t)
        pref["init_guess"] = np.vstack((ig_cos,ig_sin)) # initial guess = unit circle
        F["Xdim"] = 2*pref["N_FK"]
        F["xdag"] = np.matrix(np.zeros(F["Xdim"])).getH()
    else:
        print('unknown boundary curve type')

    F["applyGramX"] = Obstacle2dBaseOp_applyGramX
    F["applyGramX_inv"] = Obstacle2dBaseOp_applyGramX_inv
    F["applyGramY"] = Obstacle2dBaseOp_applyGramY
    F["other_X_err"] = [Obstacle2dBaseOp_L2err,'L2err']

return F, pref
