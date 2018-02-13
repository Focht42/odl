import Obstacle2dBaseOp
import ITReg.util.complete_parameter_structure

def ScatOp(p):
    F, pref = Obstacle2dBaseOp.Obstacle2dBaseOp(p)

    pref["kappa"] = 3 # wave number
    pref["N_ieq"] = 32 # 2*F_ieq is the number of discretization points
    pref["N_ieq_synth"] = 32 # 2*N_ieq is the number of discretization points for
    # the boundary integral equation when computing synthetic data (choose different
    # to N_ieq to avoid invere crime)
    pref["meas_directions"] = 64 # measurement directions
    pref["inc_directions"] = [1;0]
    # the characters of field plots encode the following:
    # X plot of reconstruction
    # x plot of initial guess/exact solution
    plotWhat["plots"] = 'XY'
    plotWhat["plot_start"] = true
    plotWhat["plot_arrows"] = false
    plotWhat["plot_only_real_farfield"] = true
    plotWhat["Nplot"] = 128
    plotWhat["field"] = false #plot total field for true obstacle?
    # If yes, the following 4 parameter are needed:
    plotWhat["ptsx"]=100 #number of points on x axis
    plotWhat["ptsy"]=100 #number of points on y axis
    plotWhat["xintv"] = [-2.5,2.5] #interval on x axis
    plotWhat["yintv"] = [-2.5,2.5] #interval on y axis

    if "plotWhat" in p.keys():
        pref["plotWhat"] = ITReg.util.complete_parameter_structure.complete_parameter_structure(p["plotWhat"],plotWhat)
    else:
        pref["plotWhat"] = plotWhat

    #TODO: adjust this function wrapper to the plotter used!
    def ScatOp_plot(F,x_k,x_start,y_k,y_obs,k):
        def NYI():
            print("not yet implemented!")

        F = NYI
        #F = scat_plot(F,x_k,x_start,y_k,y_obs,k)
        return F

    F["plot"] = ScatOp_plot

    return F, pref
