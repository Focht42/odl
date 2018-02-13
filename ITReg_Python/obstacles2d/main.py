# This is just a file to let experiments run.
# slightly altered to not need parameters.

# currently only works for diri_IRGNM_CG
import diri_IRGNM_CG
import DirichletOp
import ITReg.util.complete_itreg_par
import ITReg.IRGNM_CG

# get parameter structures
(p,reg) = diri_IRGNM_CG.diri_IRGNM_CG()

reg = ITReg.util.complete_itreg_par.complete_itreg_par(reg)

# create the operator
F = 0
if p["op_name"] == "DirichletOP":
    (F,p) = DirichletOp.DirichletOp(p)
else:
    print("Aborting: Unknown Operator")
    exit()
# show the new parameter structures, completed by default values
if reg["verbose"]>=1:
    parameter_filename = "diri_IRGNM_CG"
    print('The problem parameters specified in ',parameter_filename,' were completed as follows:')
    print(p)
    print('The regularization parameters specified in ',parameter_filename,' were completed as follows:')
    print(reg)

# save p and reg
# TODO: save([parameter_path, parameter_filename, '.par.mat'],'p','reg');

if F["syntheticdata_flag"]:
    # compute synthetic data
    #TODO: find this function
    [y_obs,F] = F.create_synthetic_data(F);
else:
     y_obs = F["y_obs"]


## perform inversion


if reg[method] == "IRGNM_CG":
    (x_recon, statistics) = ITReg.IRGNM_CG.IRGNM_CG(F["init_guess"],y_obs,F,reg)
else:
    print("Aborting: Unknown reg method")

# TODO: Save x_recon and statistics.
#save([parameter_path, parameter_filename, '.rec.mat'],'x_recon');
#save([parameter_path, parameter_filename, '.stat.mat'],'statistics');
print(x_recon)
print(statistics)

""" Original Matlab Code
function [x_recon, statistics,F,p] = main(parameter_filename)
if (nargin==0)
    parameter_filename = 'diri_IRGNM_CG';
end
parameter_path = 'experiments/';
addpath('../itreg/',parameter_path,'bd_curves/','int_op','plots');
addsubpaths('../itreg/');
close all;

%% get parameter structures
[p,reg] = feval(parameter_filename);

reg = complete_itreg_par(reg);

%% create the operator
[F,p] = feval(p.op_name,p);

% show the new parameter structures, completed by default values
if reg.verbose>=1
    disp(['The problem parameters specified in ' parameter_filename ' were completed as follows:']); p
    disp(['The regularization parameters specified in ' parameter_filename ' were completed as follows:']); reg
end
save([parameter_path, parameter_filename, '.par.mat'],'p','reg');

if F.syntheticdata_flag
    %% compute synthetic data
    [y_obs,F] = F.create_synthetic_data(F);
else
     y_obs = F.y_obs;
end

%% perform inversion
%profile on;
[x_recon, statistics] = feval(reg.method, F.init_guess,y_obs,F,reg);
%[x_recon, statistics] = aux_IRGNM_CG(p.init_guess,y_obs,F,reg);
%profile viewer;

save([parameter_path, parameter_filename, '.rec.mat'],'x_recon');
save([parameter_path, parameter_filename, '.stat.mat'],'statistics');
end

"""
