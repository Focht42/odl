function [problem_param,reg_par] = test_param()
%% =========================================================================
%               Parameters for operator
% =========================================================================

% param contains the parameters of the operator parallel_MRI;
% the meaning of these parameters is explained there

Nx = 128;
Ny = 96;
nrcoils = 4;

undersampling = 2;
nrcenterlines =  16;

% create sampling pattern P
P = zeros(Nx,Ny);
P(:,1:undersampling:end) = 1;
P(:,(Ny-nrcenterlines)/2:(Ny+nrcenterlines)/2) = 1;
samplingIndx = find(P>0);
clear P;
problem_param = struct('Nx',Nx,'Ny',Ny,'nrcoils',nrcoils,'samplingIndx',samplingIndx);
problem_param.syntheticdata_flag = true;

% simulate spin density rho
x = linspace(-1,1,Nx);
y = linspace(-1,1,Ny);
[Y,X] = meshgrid(y,x);
rho = (X.^2+Y.^2<0.8^2);

cj = simulate_coil_profiles(Nx,Ny,nrcoils);
cjFK = zeros(size(cj));
for j=1:nrcoils
	 cjFK(:,:,j) = fftshift(fft2(fftshift(cj(:,:,j)))) / sqrt(Nx*Ny);
end
problem_param.xdag = [rho(:);cjFK(:)];

%% =========================================================================
%               Parameters for plotting
% in what the letters encode the following
%    'r' real part
%    'i' imaginary part
%    'a' absolute value
%    'F' Fourier coefficients
% where encodes the corresponding figure numbers
plotWhat.rho_what='ap';
plotWhat.rho_where = [11 12];

plotWhat.coils_where = [1 2 3 4];
% the entries of coils_which are coil numbers. 0 stands for the sum
plotWhat.coils_which{1} = [1:4];   %coil profiles plotted in figure coils_where(1)
plotWhat.coils_which{2} = [1:4];   %coil profiles plotted in figure coils_where(2)
plotWhat.coils_which{3} = [0];      %coil profiles plotted in figure coils_where(3)
plotWhat.coils_which{4} = [0];      %coil profiles plotted in figure coils_where(4)
plotWhat.coils_what{1}='aaaa'; %plot what in figure coils_where(1)
plotWhat.coils_what{2}='pppp'; %plot what in figure coils_where(2)
plotWhat.coils_what{3}='a';           %plot what in figure coils_where(3)
plotWhat.coils_what{4}='p';           %plot what in figure coils_where(4)
% figure coils_where(l) has n1(l) x n2(l) subplots
plotWhat.coils_n1 = [2 2 1 1];
plotWhat.coils_n2 = [2 2 1 1];

problem_param.plotWhat = plotWhat;

%% =========================================================================
%               Parameter defining the regularization method

% regularization method
reg_par.method = 'IRGNM_CG';
%'IRGNM_CG' % 'Landweber'; %'Newton_CG'; 

% maximum number of iterations
reg_par.N_max_it =10;
% maximum number of CG iterations
reg_par.N_CG = 500;
% starting reg. parameter for IRGNM
reg_par.alpha0 = 1;
% decreasing step for the regulation parameter (IRGNM)
reg_par.alpha_step = 1/3;
end
