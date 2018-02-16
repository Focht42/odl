"""
function [problem_param,reg] = ksp3x2_irgnm()
%% ========================================================================
%           parameter defining the operator

load('ksp3x2.mat');

[Nx,Ny,nrcoils] = size(Y);

samplingIndx = find(sum(abs(Y),3)>0>0);

problem_param = struct('Nx',Nx,'Ny',Ny,'nrcoils',nrcoils,'samplingIndx',samplingIndx);
problem_param.syntheticdata_flag = false;

N = length(samplingIndx);
problem_param.data = zeros(N*nrcoils,1);
for j=1:nrcoils
    aux = Y(:,:,j);
    problem_param.data((j-1)*N+1:j*N) = aux(samplingIndx);
end

%% =========================================================================
%               Parameters for plotting
% in what the letters encode the following
%    'r' real part
%    'i' imaginary part
%    'a' absolute value
%    'p' phase
%    'F' Fourier coefficients
% where encodes the corresponding figure numbers
plotWhat.rho_what='pa';
plotWhat.rho_where = [11 12];

%plotWhat.coils_where = [1 2];
% the entries of coils_which are coil numbers. 0 stands for the sum
%plotWhat.coils_which{1} = [1:12];   %coil profiles plotted in figure coils_where(1)
%plotWhat.coils_which{2} = [1:12];   %coil profiles plotted in figure coils_where(1)
%plotWhat.coils_what{1}='aaaaaaaaaaaa'; %plot what in figure coils_where(1)
%plotWhat.coils_what{2}='pppppppppppp'; %plot what in figure coils_where(2)
% figure coils_where(l) has n1(l) x n2(l) subplots
%plotWhat.coils_n1 = [3 3];
%plotWhat.coils_n2 = [4 4];
plotWhat.coils_where = [];
plotWhat.coils_which = [];
plotWhat.coils_what = [];
plotWhat.coils_n1 = [1];
plotWhat.coils_n2 = [1];
problem_param.plotWhat = plotWhat;


%% =========================================================================
%               Parameters for the regularization method

% regularization method
reg.method = 'IRGNM_CG';
%'IRGNM_CG' % 'Landweber'; %'Newton_CG';

% maximum number of iterations
reg.N_max_it =8;
% maximum number of CG iterations
reg.max_CG_steps = 400;
% starting reg. parameter for IRGNM
reg.alpha0 = 1;
% decreasing step for the regulation parameter (IRGNM)
reg.alpha_step = 1/3;
reg.stop_rule = 'maxit';
reg.CG_TOL_rel_accY = 0.3;
reg.CG_TOL_rel_accX = 0.3;
reg.verbose = 2;
end
"""

from scipy.io import loadmat
import numpy as np

def ksp3x2_irgnm():
    #%% ========================================================================
    #%           parameter defining the operator

    matdata = loadmat('simulations/ksp3x2.mat') #NOTE: hope this is the right file path!
    Y = matdata['Y']
    [Nx,Ny,nrcoils] = Y.shape

    # samplingIndx = find(sum(abs(Y),3)>0>0); NOTE: I don't get the >0>0 here :(
    samplingIndx = np.nonzero(sum(abs(Y))>0)

    problem_param = {"Nx":Nx, "Ny":Ny, "nrcoils": nrcoils, "samplingIndx": samplingIndx}
    problem_param["syntheticdata_flag"] = False

    N = np.size(samplingIndx,1) # numpy stores x and y as 2 arrays!
    problem_param["data"] = np.zeros(N*nrcoils, dtype=np.complex_)
    for j in range(nrcoils):
        aux = Y[:,:,j]
        problem_param["data"][j*N:(j+1)*N] = aux[samplingIndx]

    '''
    %% =========================================================================
    %               Parameters for plotting
    % in what the letters encode the following
    %    'r' real part
    %    'i' imaginary part
    %    'a' absolute value
    %    'p' phase
    %    'F' Fourier coefficients
    % where encodes the corresponding figure numbers
    '''
    plotWhat = dict()
    plotWhat["rho_what"]='pa'
    plotWhat["rho_where"] = np.array([11,12])

    '''
    %plotWhat.coils_where = [1 2];
    % the entries of coils_which are coil numbers. 0 stands for the sum
    %plotWhat.coils_which{1} = [1:12];   %coil profiles plotted in figure coils_where(1)
    %plotWhat.coils_which{2} = [1:12];   %coil profiles plotted in figure coils_where(1)
    %plotWhat.coils_what{1}='aaaaaaaaaaaa'; %plot what in figure coils_where(1)
    %plotWhat.coils_what{2}='pppppppppppp'; %plot what in figure coils_where(2)
    % figure coils_where(l) has n1(l) x n2(l) subplots
    %plotWhat.coils_n1 = [3 3];
    %plotWhat.coils_n2 = [4 4];
    '''
    plotWhat["coils_where"] = []
    plotWhat["coils_which"] = []
    plotWhat["coils_what"] = []
    plotWhat["coils_n1"] = [1]
    plotWhat["coils_n2"] = [1]
    problem_param["plotWhat"] = plotWhat

    '''
    %% =========================================================================
    %               Parameters for the regularization method
    '''
    reg = dict()
    #% regularization method
    reg["method"] = 'IRGNM_CG'
    #%'IRGNM_CG' % 'Landweber'; %'Newton_CG';

    #% maximum number of iterations
    reg["N_max_it"] =8
    #% maximum number of CG iterations
    reg["max_CG_steps"] = 400
    #% starting reg. parameter for IRGNM
    reg["alpha0"] = 1
    #% decreasing step for the regulation parameter (IRGNM)
    reg["alpha_step"] = 1/3
    reg["stop_rule"] = 'maxit'
    reg["CG_TOL_rel_accY"] = 0.3
    reg["CG_TOL_rel_accX"] = 0.3
    reg["verbose"] = 2

    return problem_param,reg
