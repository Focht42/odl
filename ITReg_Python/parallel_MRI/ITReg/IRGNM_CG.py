
"""
function [x_k, stat,F] = IRGNM_CG(x_start,y_obs,F,par)
% iteratively regularized Gauss Newton method
% the regularized Newton equations are solved by the conjugate gradient
% method applied to the normal equation
%
% x_start: initial guess
% y_obs  : observed data
% F      : forward operator (see README.txt)
% par    : structure containing the following parameters:
%     N_max_it            maximum number of outer Newton iterations
%     max_CG_steps        maximum number of CG iteration
%     alpha0,alpha_step   regularization parameter in k-th Newton step
%                  is alpha0 * alpha_step^k
%     Update:

par.Update = @StandardUpdate;
verb =  par.verbose;
stoprule = feval(par.stop_rule,par);
par.it_step = 0;

x_k = x_start;
[y_k,F] = F.evaluate(F,x_start);

[stat,F] = error_stat_output(F,verb>=1,[],x_k,y_obs,y_k,x_start,par);
condPrint(verb>=1,', alpha=%1.3e\n',par.alpha0);
stat.CG_steps = [];

[mustExit, stoprule] = stoprule.stop(stoprule,F,x_k,y_obs,par);

while ~mustExit
    par.it_step = par.it_step+1;
    regpar = par.alpha0 * par.alpha_step^(par.it_step);

    [update, CG_steps] = CGNE_reg(F,y_obs-y_k,x_start-x_k,regpar,par);

    x_k = par.Update(x_k, update,F);

    [y_k,F] = F.evaluate(F,x_k);

    [stat,F] = error_stat_output(F,verb>=1,stat,x_k,y_obs,y_k,x_start,par);
    condPrint(verb>=1,', alpha=%1.1e, CGstp = %i\n',regpar,CG_steps);
    stat.CG_steps = [stat.CG_steps, CG_steps];
    condPrint(par.it_step==par.N_max_it & verb>=1,...
        'Maximum number of GN iterations reached.\n');
    [mustExit, stoprule] = stoprule.stop(stoprule,F,x_k,y_obs,par);
end

[stop_ind,x_k] = stoprule.select_index(stoprule,F);
condPrint(verb>=1,'The stopping rule has selected the index %u\n',stop_ind);
stat.chosen_rec = stop_ind;

%% plotting
if ~isempty(find(par.plot_steps==-1))
    F = F.plot(F,x_k,x_start,y_k,y_obs,-stop_ind);
end
end

%% ------------------------------------------------------------------------
function newx = StandardUpdate(oldx,h,F)
    newx = oldx+h;
end

%% -------------------------------------------------------------------------
function [h, CGstep] = CGNE_reg(F,y,xref,regpar,param)
% solves A h = b by CGNE with
%    A := G_X^{-1} F'* G_Y F' + regpar I
%    b := G_X^{-1}F'^* G_Y y + regpar xref
% A is self-adjoint with respect to the inner product <u,v> = u'G_X v
%
% G_X, G_Y -> F.applyGramX, F.applyGramY
% G_X^{-1} -> F.applyGramX_inv
% F'       -> F.derivative
% F'*      -> F.adjoint

verb=param.verbose;
max_CG_steps = param.max_CG_steps;
CG_TOL_rel_accX = param.CG_TOL_rel_accX;
CG_TOL_rel_accY = param.CG_TOL_rel_accY;
CG_TOL_resi_red = param.CG_TOL_resi_red;

CGstep =0;
kappa = 1;
% compute stilde = G_X b
ztilde = F.applyGramY(F,y);
stilde = F.adjoint(F,ztilde);
stilde = stilde + regpar * F.applyGramX(F,xref);

s = F.applyGramX_inv(F,stilde);
d = s;
dtilde = stilde;
norm_s = real(stilde'*s);
norm_s0 = norm_s;
norm_h = 0;

h = zeros(size(s));
Th = zeros(size(y));
Thtilde = zeros(size(y));

while (sqrt(norm_s/norm_h/kappa)/regpar > CG_TOL_rel_accX/(1+CG_TOL_rel_accX) && ...
       sqrt(norm_s/real(Thtilde'*Th)/kappa/regpar) > CG_TOL_rel_accY/(1+CG_TOL_rel_accY) && ...
       sqrt(norm_s/norm_s0/kappa)> CG_TOL_resi_red && ...
       CGstep <= max_CG_steps)
    z = F.derivative(F,d);
    ztilde = F.applyGramY(F,z);
    gamma = norm_s / real(regpar*dtilde'*d+ztilde'*z);
    h      = h + gamma * d;
    Th     = Th + gamma * z;
    Thtilde = Thtilde + gamma* ztilde;
    stilde = stilde - gamma * (F.adjoint(F,ztilde) + regpar*dtilde);
    s      = F.applyGramX_inv(F,stilde);
    norm_s_old = norm_s;
    norm_s = real(stilde'*s);
    beta   = norm_s/norm_s_old;
    d      = s + beta * d;
    dtilde = stilde + beta * dtilde;
    norm_h = h' * F.applyGramX(F,h);
    kappa = 1+beta * kappa;
    CGstep= CGstep + 1;
    condPrint(verb>=2,'  %i: resi-red. %3.1e (%3.1e),  rel.accY %3.1e (%3.1e),  rel.accX %3.1e (%3.1e), sqrt kappa = %3.1f \n',...
        CGstep, sqrt(norm_s/norm_s0), sqrt(norm_s/norm_s0/kappa), ...
        sqrt(norm_s/real(Thtilde'*Th)/regpar), sqrt(norm_s/real(Thtilde'*Th)/kappa/regpar), ...
        sqrt(norm_s/norm_h)/regpar, sqrt(norm_s/norm_h/kappa)/regpar, sqrt(kappa));
end
end
"""
import ITReg.utilities.error_stat_output as eso
import ITReg.utilities.condPrint as cndP
import numpy as np

def IRGNM_CG(x_start,y_obs,F,par):
    """
    % iteratively regularized Gauss Newton method
    % the regularized Newton equations are solved by the conjugate gradient
    % method applied to the normal equation
    %
    % x_start: initial guess
    % y_obs  : observed data
    % F      : forward operator (see README.txt)
    % par    : structure containing the following parameters:
    %     N_max_it            maximum number of outer Newton iterations
    %     max_CG_steps        maximum number of CG iteration
    %     alpha0,alpha_step   regularization parameter in k-th Newton step
    %                  is alpha0 * alpha_step^k
    %     Update:
    """
    par["Update"] = StandardUpdate
    verb =  par["verbose"]
    stoprule = par["stop_rule"](par)
    par["it_step"] = 0

    x_k = x_start
    [y_k,F] = F["evaluate"](F,x_start)

    [stat,F] = eso.error_stat_output(F,verb>=1,dict(),x_k,y_obs,y_k,x_start,par)
    cndP.condPrint(verb>=1,', alpha=%1.3e\n',par["alpha0"])
    stat["CG_steps"] = []

    [mustExit, stoprule] = stoprule["stop"](stoprule,F,x_k,y_obs,par)

    while not mustExit:
        par["it_step"] = par["it_step"]+1
        regpar = par["alpha0"] * par["alpha_step"]**(par["it_step"])

        [update, CG_steps] = CGNE_reg(F,y_obs-y_k,x_start-x_k,regpar,par)

        x_k = par["Update"](x_k, update,F)

        [y_k,F] = F["evaluate"](F,x_k)

        [stat,F] = eso.error_stat_output(F,verb>=1,stat,x_k,y_obs,y_k,x_start,par)
        cndP.condPrint(verb>=1,', alpha=%1.1e, CGstp = %i\n',regpar,CG_steps)
        stat["CG_steps"] = [stat["CG_steps"], CG_steps]
        cndP.condPrint(par["it_step"]==par["N_max_it"] & verb>=1,'Maximum number of GN iterations reached.\n')
        [mustExit, stoprule] = stoprule["stop"](stoprule,F,x_k,y_obs,par)


    [stop_ind,x_k] = stoprule["select_index"](stoprule,F)
    cndP.condPrint(verb>=1,'The stopping rule has selected the index %u\n',stop_ind)
    stat["chosen_rec"] = stop_ind

    #%% plotting
    #TODO
    """
    if ~isempty(find(par["plot_steps"]==-1))
        F = F["plot"](F,x_k,x_start,y_k,y_obs,-stop_ind);
    end
    """

    return x_k, stat,F

#%% ------------------------------------------------------------------------
def StandardUpdate(oldx,h,F):
    newx = oldx+h
    return newx

#%% -------------------------------------------------------------------------
def CGNE_reg(F,y,xref,regpar,param):
    """
    % solves A h = b by CGNE with
    %    A := G_X^{-1} F'* G_Y F' + regpar I
    %    b := G_X^{-1}F'^* G_Y y + regpar xref
    % A is self-adjoint with respect to the inner product <u,v> = u'G_X v
    %
    % G_X, G_Y -> F.applyGramX, F.applyGramY
    % G_X^{-1} -> F.applyGramX_inv
    % F'       -> F.derivative
    % F'*      -> F.adjoint
    """

    verb=param["verbose"]
    max_CG_steps = param["max_CG_steps"]
    CG_TOL_rel_accX = param["CG_TOL_rel_accX"]
    CG_TOL_rel_accY = param["CG_TOL_rel_accY"]
    CG_TOL_resi_red = param["CG_TOL_resi_red"]

    CGstep =0
    kappa = 1
    #% compute stilde = G_X b
    ztilde = F["applyGramY"](F,y)
    stilde = F["adjoint"](F,ztilde)
    stilde = stilde + regpar * F["applyGramX"](F,xref)

    s = F["applyGramX_inv"](F,stilde)
    d = s
    dtilde = stilde
    print(s.shape)
    print(stilde.shape)
    norm_s = np.real(stilde.conj().transpose().dot(s))
    norm_s0 = norm_s
    norm_h = 0

    h = np.zeros(s.shape) #TODO: verify that we need a shaped matrix and not just an array
    Th = np.zeros(y.shape) #TODO: verify that we need a shaped matrix and not just an array
    Thtilde = np.zeros(y.shape) #TODO: verify that we need a shaped matrix and not just an array


    while (np.sqrt(norm_s/norm_h/kappa)/regpar > CG_TOL_rel_accX/(1+CG_TOL_rel_accX) and
            np.sqrt(norm_s/np.real(Thtilde.conj().transpose().dot(Th))/kappa/regpar) > CG_TOL_rel_accY/(1+CG_TOL_rel_accY) and
            np.sqrt(norm_s/norm_s0/kappa)> CG_TOL_resi_red and
            CGstep <= max_CG_steps):

        #TODO: Check if those * are .dot() or element wise
        z = F["derivative"](F,d)
        ztilde = F["applyGramY"](F,z)
        gamma = norm_s / np.real(regpar * dtilde.conj().transpose().dot(d) + ztilde.conj().transpose().dot(z))
        h      = h + gamma * d
        Th     = Th + gamma * z
        Thtilde = Thtilde + gamma* ztilde
        stilde = stilde - gamma * (F["adjoint"](F,ztilde) + regpar*dtilde)
        s      = F["applyGramX_inv"](F,stilde)
        norm_s_old = norm_s
        norm_s = np.real(stilde.conj().transpose().dot(s))
        beta   = norm_s/norm_s_old
        d      = s + beta * d
        dtilde = stilde + beta * dtilde
        norm_h = h.conj().transpose().dot(F["applyGramX"](F,h))
        kappa = 1+beta * kappa
        CGstep= CGstep + 1
        cndP.condPrint(verb>=2,'  %i: resi-red. %3.1e (%3.1e),  rel.accY %3.1e (%3.1e),  rel.accX %3.1e (%3.1e), sqrt kappa = %3.1f \n', CGstep, np.sqrt(norm_s/norm_s0), np.sqrt(norm_s/norm_s0/kappa), np.sqrt(norm_s/np.real(Thtilde.conj().transpose().dot(Th))/regpar), np.sqrt(norm_s/np.real(Thtilde.conj().transpose().dot(Th))/kappa/regpar), np.sqrt(norm_s/norm_h)/regpar, np.sqrt(norm_s/norm_h/kappa)/regpar, np.sqrt(kappa))

    return h, CGstep
