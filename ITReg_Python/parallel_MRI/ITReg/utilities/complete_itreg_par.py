import ITReg.utilities.complete_parameter_structure as cps

def complete_itreg_par(par):
    ref = dict()
    par_new = 0

    ref["method"] = par["method"]
    ref["stoprule"] = "discrepancy"
    # parameter of the discrepancy principle
    ref["tau"] = 2
    # parameter of the Lepskij principle
    ref["c1"] = 0.001;
    ref["c2"] = 50
    # maximum number of iterations
    ref["N_max_it"] = 50
    # numbers of the reconstructions to be plot_total_field
    ref["plot_steps"] = list(range(-1,ref["N_max_it"]+1,1))
    # numbers of reconstructions to be stored
    ref["store_rec"] = []
    # the higher the value, the more messages are given
    ref["verbose"] = 1

    inner_params = par.pop("inner_params",[]) # removes key. If not present returns []

    #TODO: just making sure
    if par == 0:
        print("Aborting: error in complete_itreg_par")
        exit(1)

    #this is the switch case statement of the matlab file
    if par["method"] in ["IRGNM_CG", "IRGNM", "constr_IRGNM_CG"]:
        #starting reg. parameter for IRGNM
        ref["alpha0"] = 1
        # decreasing step for the regularization parameter (IRGNM)
        ref["alpha_step"] = 2/3
        # maximum number of CG steps
        ref["max_CG_steps"] = 50
        # controls the relative accuracy of the Newton update in preimage space
        ref["CG_TOL_rel_accX"] = 0.3
        # controls the relative accuracy of the Newton update in data space
        ref["CG_TOL_rel_accY"] = 0.3;
        # controls the reduction of the residual
        ref["CG_TOL_resi_red"] = 1e-6;
        par_new = cps.complete_parameter_structure(par,ref)

    elif par["method"] in ['IRGNM_L1_fid']:
        print("not yet implemented")
        """
        % starting reg. parameter for IRGNM
        ref.alpha0 = 1;
        % decreasing step for the regulation parameter (IRGNM)
        ref.alpha_step = 2/3;
        reg.alphaL1 = 1e-4;
        par_new = complete_parameter_structure(par,ref);
        """

    elif par["method"] in ['Newton_CG']:
        print("not yet implemented")
        """
        % rho for abort criterion of CG iteration
        ref.rho = 0.8;
        % maximum number of CG steps
        ref.max_CG_steps = 50;
        par_new = complete_parameter_structure(par,ref);
        """
    elif par["method"] in ['IRNM_KL']:
        print("not yet implemented")
        """
            % starting reg. parameter for IRNM
            ref.alpha0 = 5e-6;
            % decreasing step for the regulation parameter (IRNM)
            ref.alpha_step = 2/3;
        	ref.inner_solver_name = 'SQP';
        	par_new = complete_parameter_structure(par,ref);
        	par_new.inner_solver = feval(par_new.inner_solver_name,inner_params);
            par_new.inner_solver.params.verbose = par_new.verbose;

        """
    elif par["method"] in ['IRNM_KL_Newton', 'IRNM_KL_Newton_point_eval']:
        print("not yet implemented")
        """
            % starting reg. parameter for IRNM
            ref.alpha0 = 2e-6;
            % decreasing step for the regulation parameter (IRNM)
            ref.alpha_step = 2/3;
            % maximum number of CG steps
            ref.max_CG_steps = 50;
            ref.offset0 =1e-4;
            % offset is reduced in each Newton step by a factor offset_step
            ref.offset_step = 0.8;
            % relative tolerance value for termination of inner Newton iteration
            ref.inner_residual = 1e-10;
            % max number of inner iterations to minimize the KL-functional
            ref.inner_it = 10;
            par_new = complete_parameter_structure(par,ref);

        """
    elif par["method"] in ['Landweber']:
        print("not yet implemented")
        """
        par_new = complete_parameter_structure(par,ref);
        """
    elif par["method"] in ['TikREG']:
        print("not yet implemented")
        """
        % reg. parameter
        ref.alpha0 = 5e-6;
    	ref.inner_solver_name = 'SQP';
    	par_new = complete_parameter_structure(par,ref);
    	par_new.inner_solver = feval(par_new.inner_solver_name,inner_params);
        par_new.inner_solver.params.verbose = par_new.verbose;
        """
    elif par["method"] in ['CONST_IRGNM']:
        print("not yet implemented")
        """
        ref.alpha0 = 1;
        ref.alpha_min = 10^-10;
        ref.alpha_step = 2/3;

        ref.it_const_solver = @ssNewton;
        % SSN parameters
        % maximum number of SSN-iterations
        ref.ssn.max_it = 100;
        % to compute the update of the new CG-tolerance
        ref.ssn.TolFac = 2;
        % inactive set initialization
        ref.ssn.inactive = [];

        % ssN stopping by duality gap estimation (1 on, 0 off)
        ref.ssn.duality_gap = 1;
        % the size of the duality gap when ssN will stop
        ref.ssn.eps = 10^-2;

        ref.it_inner_solver = @CG_selfadjoint;
        % CG-parameters
        % maximum number of CG steps
        ref.CG.max_CG_steps = 1000;
        % final CG stopping tolerance
        ref.CG.Tol_CG = 10^-14;
        % initial CG stopping tolerance
        ref.CG.Tol_CG_start = 10;
        par_new = complete_parameter_structure(par,ref);
        """
    else:
        #TODO: real error_stat_output
        print("Error: unknown regularization method")

    return par_new
