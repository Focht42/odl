
"""
function [x_recon, statistics] = main(parameter_filename)
parameter_dir = 'simulations';
if (nargin==0)
    parameter_filename = 'ksp3x2_irgnm';
end
addpath('../itreg/',parameter_dir);
addsubpaths('../itreg/');
close all;

%% get parameter structures
[problem_param,reg_par] = feval(parameter_filename);
reg_par = complete_itreg_par(reg_par);

%% create operator
[F,problem_param] = parallel_MRI(problem_param);

%% compute synthetic data
if F.syntheticdata_flag
    data = F.create_synthetic_data(F);
else
    % normalize data
    yscal = 100/norm(problem_param.data(:));
    data = yscal * problem_param.data(:);
end

%% perform inversion
tic
[x_recon, statistics] = feval(reg_par.method, problem_param.init_guess,...
    data,F,reg_par);
toc
end
"""

import simulations.ksp3x2_irgnm
import ITReg.utilities.complete_itreg_par
import ITReg.stopping_rules.maxit
import ITReg.IRGNM_CG
import parallel_MRI
import numpy as np

def main():
    #start parse command line parameter
    #TODO: Do this maybe at some point. For now we will default to ksp3x2_irgnm
    #end parse command line parameter

    #%% get parameter structures

    [problem_param,reg_par] = simulations.ksp3x2_irgnm.ksp3x2_irgnm()
    reg_par = ITReg.utilities.complete_itreg_par.complete_itreg_par(reg_par)
    #%% create operator
    [F,problem_param] = parallel_MRI.parallel_MRI(problem_param)
    #NOTE: replace strings with function handles!
    reg_par["method"] = ITReg.IRGNM_CG.IRGNM_CG
    reg_par["stop_rule"] = ITReg.stopping_rules.maxit.maxit
    #NOTE: Done replacing.



    #%% compute synthetic data
    if F["syntheticdata_flag"]:
        print("Main: Compute Synthetic Data")
        data = F["create_synthetic_data"](F)
    else:
        #% normalize data
        print("Main: normalize Data")
        yscal = 100/np.linalg.norm(problem_param["data"].reshape(problem_param["data"].size))
        data = yscal * problem_param["data"].reshape(problem_param["data"].size)

    #%% perform inversion
    #tic
    print("tic")
    [x_recon, statistics] = reg_par["method"](problem_param["init_guess"],data,F,reg_par)
    print("toc")
    #toc


if __name__ == "__main__":
    main()
