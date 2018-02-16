import numpy as np
import ITReg.utilities.complete_parameter_structure as cps

'''
Python translation of parallel_MRI.m with Matlab code in comments
Also includes the plotter functions which aren't yet implemented
'''
#function [F,p]=parallel_MRI(p)
def parallel_MRI(p):

    '''
    % Reference: M. Uecker, T. Hohage, K.T. Block and J. Frahm.
    % Image Reconstruction by Regularized Nonlinear Inversion -
    % Joint Estimation of Coil Sensitivities and Image Content.
    % Magnetic Resonance in Medicine, 60:674-682, 2008.
    '''
    F = dict() # Define F as a dict
    #F.op_name = 'parallel_MRI';
    F["op_name"] = 'parallel_MRI'
    #F.syntheticdata_flag = false;
    F["syntheticdata_flag"] = False

    #% number of receiver coils
    #F.nrcoils = 12;
    F["nrcoils"] = 12
    #% the unknown spin density is an image of size Nx*Ny
    #F.Nx = 128;
    F["Nx"] = 128
    #F.Ny = 96;
    F["Ny"] = 96
    #% vector of indices of Fourier coefficients, which are sampled
    #F.samplingIndx = [];
    F["samplingIndx"] = [] #TODO: Determine if this is the correct data type

    #if p.syntheticdata_flag:
        #F.xdag = []
    if p["syntheticdata_flag"]:
        F["xdag"] = []
    #F.data = []

    F["data"] = [] #TODO: Determine if this is the correct data type
    '''
    % in what the letters encode the following
    %    'r' real part
    %    'i' imaginary part
    %    'a' absolute value
    %    'F' Fourier coefficients
    % where encodes the corresponding figure numbers
    '''
    plotWhat = dict()
    #plotWhat.rho_what='a';
    #plotWhat.rho_where = [11];

    plotWhat["rho_what"]='a'
    plotWhat["rho_where"] = [11]

    #plotWhat.coils_where = [1];
    plotWhat["coils_where"] = [1]

    #% the entries of coils_which are coil numbers. 0 stands for the sum
    #plotWhat.coils_which{1} = [1];   %coil profiles plotted in figure coils_where(1)
    #plotWhat.coils_what{1}='a'; %plot what in figure coils_where(1)
    #% figure coils_where(l) has n1(l) x n2(l) subplots
    #plotWhat.coils_n1 = [1];
    #plotWhat.coils_n2 = [1];
    #F.plotWhat = plotWhat;
    #if isfield(p,'plotWhat')
    #    p.plotWhat = complete_parameter_structure(p.plotWhat,F.plotWhat);

    #NOTE: I believe the {} scalar field above is just used as a list/dict
    #plotWhat["coils_which"] = {1:[1]}  #%coil profiles plotted in figure coils_where(1)#TODO: this causes errors
    #plotWhat["coils_what"] = {1:'a'}  #%plot what in figure coils_where(1) #TODO: this causes errors
    #% figure coils_where(l) has n1(l) x n2(l) subplots
    plotWhat["coils_n1"] = [1]
    plotWhat["coils_n2"] = [1]
    F["plotWhat"] = plotWhat
    if 'plotWhat' in p.keys():
        p["plotWhat"] = cps.complete_parameter_structure(p["plotWhat"],F["plotWhat"]) #TODO: Locate Function



    #F = complete_parameter_structure(p,F); #TODO: Locate Function
    #p=F;
    #rmfield(F,'data');

    F = cps.complete_parameter_structure(p,F) #TODO: Locate Function
    p=F.copy()
    F.pop('data')


    '''
    if isempty(F.samplingIndx)
        undersampling = 2;
        nrcenterlines =  16;
        % create sampling pattern P
        P = zeros(Nx,Ny);
        P(:,1:undersampling:end) = 1;
        P(:,(Ny-nrcenterlines)/2:(Ny+nrcenterlines)/2) = 1;
        samplingIndx = find(P>0);
        clear P;
    '''
    #TODO: Try/Catch stuff here
    if np.size(F["samplingIndx"]) <= 0:
        undersampling = 2
        nrcenterlines =  16
        #% create sampling pattern P
        P = np.zeros((F["Nx"],F["Ny"])) #Nx and Ny are not defined but I assume they are F["Nx"] and F["Ny"]
        P[:,0::undersampling] = 1
        P[:,int((Ny-nrcenterlines)/2-1):int((Ny+nrcenterlines)/2-1)] = 1 #-1 because of Matlab-1-indexing, also cast to int. TODO:Check if this actually selects the right elements
        samplingIndx = np.nonzero(P>0) # this now returns a tuple of two arrays that contain the x and y indices TODO: Check if this causes problems
        del P #P not needed anymore

    '''
    if ~isfield(p,'init_guess')
        N = p.Nx*p.Ny;
        p.init_guess = zeros(N*(1+p.nrcoils),1);
        p.init_guess(1:N) = 1;
    '''

    if not 'init_guess' in p.keys():
        N = p["Nx"] * p["Ny"]
        p["init_guess"] = np.zeros((N*(1+p["nrcoils"]),1))
        p["init_guess"][0:N,0] = 1 #NOTE: I feel like i need to reshape this at some point.

    '''
    F.Xdim = F.Nx * F.Ny * (F.nrcoils+1);
    F.Ydim = length(F.samplingIndx) * F.nrcoils;
    x = [-F.Nx/2:1:F.Nx/2-0.01]/F.Nx;
    y = [-F.Ny/2:1:F.Ny/2-0.01]/F.Ny;
    [Y,X] = meshgrid(y,x);
    '''
    F["Xdim"] = F["Nx"] * F["Ny"] * (F["nrcoils"]+1)
    F["Ydim"] = np.size(F["samplingIndx"],1) * F["nrcoils"] # need to get the correct size for samplingIndx
    x = np.arange(-F["Nx"]/2,F["Nx"]/2-0.01,1)/F["Nx"]
    y = np.arange(-F["Ny"]/2,F["Ny"]/2-0.01,1)/F["Ny"]
    [Y,X] = np.meshgrid(y,x)

    '''
    %F.Fourier_weights = (1. + 20*(X.^2+Y.^2)).^10;
    F.Fourier_weights = zeros(F.Nx,F.Ny);
    for l = [1:F.Nx],
        for j = [1:F.Ny],
            d = ((l - 1) / F.Nx - 0.5)^2 + ((j - 1) / F.Ny - 0.5)^2;
            F.Fourier_weights(l, j) = (1. + 220. * d)^32;
    '''


    # %F.Fourier_weights = (1. + 20*(X.^2+Y.^2)).^10;
    F["Fourier_weights"] = np.zeros((F["Nx"],F["Ny"]))
    for l in range(F["Nx"]):
        for j in range(F["Ny"]):
            d = (l / F["Nx"] - 0.5)**2 + (j  / F["Ny"] - 0.5)**2 # modified due to index shift between Matlab and Python
            F["Fourier_weights"][l, j] = (1 + 220 * d)**32



    '''
    F.evaluate = @parallel_MRI_evaluate;
    F.create_synthetic_data = @parallel_MRI_create_synthetic_data;
    F.derivative = @parallel_MRI_derivative;
    F.adjoint = @parallel_MRI_adjoint;
    F.plot = @parallel_MRI_plot;
    F.applyGramX = @parallel_MRI_applyGramX;
    F.applyGramX_inv = @parallel_MRI_applyGramX_inv;
    F.applyGramY = @parallel_MRI_applyGramY;
    '''
    # Function Handles
    F["evaluate"] = parallel_MRI_evaluate #TODO: Find Function
    F["create_synthetic_data"] = parallel_MRI_create_synthetic_data
    F["derivative"] = parallel_MRI_derivative
    F["adjoint"] = parallel_MRI_adjoint #TODO: Find Function
    F["plot"] = parallel_MRI_plot #TODO: Find Function
    F["applyGramX"] = parallel_MRI_applyGramX #TODO: Find Function
    F["applyGramX_inv"] = parallel_MRI_applyGramX_inv #TODO: Find Function
    F["applyGramY"] = parallel_MRI_applyGramY #TODO: Find Function

    return F, p

'''
function data  = parallel_MRI_create_synthetic_data(F)
    data = F.evaluate(F,F.xdag);
    % normalize data
    yscal = 100/norm(data(:));
    data = yscal*data;
end
'''
def parallel_MRI_create_synthetic_data(F):
    data = F["evaluate"](F,F["xdag"])
    #% normalize data
    yscal = 100/np.linalg.norm(data.flatten()) #NOTE: flatten forces a copy. Better way?
    data = yscal*data
    return data

'''
function [data,F] = parallel_MRI_evaluate(F,rho_and_coils)
    nrcoils = F.nrcoils;
    samplingIndx = F.samplingIndx;
    N = F.Nx*F.Ny;

    % current spin density and coil sensitivities are stored for
    % later evaluations of F.derivative and F.adjoint:
    % The derivative will be taken at the point (rho,coils)
    F.rho = reshape(rho_and_coils(1:N),F.Nx,F.Ny,1);
    F.coils = reshape(rho_and_coils(N+1:end),F.Nx,F.Ny,nrcoils);

    data = zeros(length(samplingIndx),nrcoils);

    for j = 1:nrcoils
        F.coils(:,:,j) = myifft(F.coils(:,:,j));
        aux = myfft(F.rho .* F.coils(:,:,j));
        data(:,j) = aux(samplingIndx);
    end
    data = data(:);
end
'''
def parallel_MRI_evaluate(F,rho_and_coils):
    nrcoils = F["nrcoils"]
    samplingIndx = F["samplingIndx"]
    N = F["Nx"]*F["Ny"]

    #% current spin density and coil sensitivities are stored for
    #% later evaluations of F.derivative and F.adjoint:
    #% The derivative will be taken at the point (rho,coils)
    #TODO: This creates copies. It sure can be improved!
    F["rho"] = rho_and_coils.flatten()[:N].reshape(F["Nx"],F["Ny"]) # not ,1
    F["coils"] = rho_and_coils.flatten()[N:].reshape(F["Nx"],F["Ny"],nrcoils)
    data = np.zeros((np.size(samplingIndx,1),nrcoils), dtype = np.complex_)

    for j in range(nrcoils):
        F["coils"][:,:,j] = myifft(F["coils"][:,:,j])
        aux = myfft(F["rho"] * F["coils"][:,:,j])
        data[:,j] = aux[samplingIndx] #NOTE: selects elements from samplingIndx

    data = data.flatten()
    return data,F

'''
function d_K = parallel_MRI_derivative(F,h)
    nrcoils = F.nrcoils;
    samplingIndx = F.samplingIndx;
    N = F.Nx*F.Ny;

    d_K = zeros(length(samplingIndx),nrcoils);
    d_rho = reshape(h(1:N),F.Nx,F.Ny,1);
    d_coils = reshape(h(N+1:end),F.Nx,F.Ny,nrcoils);

    for j = 1:nrcoils
        aux = myfft(F.rho .* myifft(d_coils(:,:,j)) + d_rho .* F.coils(:,:,j));
        d_K(:,j) = aux(samplingIndx);
    end
    d_K = d_K(:);
end
'''
def parallel_MRI_derivative(F,h):
    nrcoils = F["nrcoils"]
    samplingIndx = F["samplingIndx"]
    N = F["Nx"]*F["Ny"]

    d_K = np.zeros((np.size(samplingIndx,1),nrcoils))
    d_rho = h.flatten()[:N].reshape(F["Nx"],F["Ny"],1)
    d_coils = h.flatten()[N:].reshape(F["Nx"],F["Ny"],nrcoils)

    for j in range(nrcoils):
        aux = myfft(F["rho"] * myifft(d_coils[:,:,j]) + d_rho * F["coils"][:,:,j])
        d_K[:,j] = aux[samplingIndx]

    d_K = d_K.flatten()
    return d_K

'''
function d_rho_and_coils = parallel_MRI_adjoint(F, d_K)
    nrcoils = F.nrcoils;
    samplingIndx = F.samplingIndx;
    Nx = F.Nx;
    Ny = F.Ny;
    M = length(samplingIndx);

    aux = zeros(Nx,Ny);
    d_rho = zeros(Nx,Ny);
    d_coils = zeros(Nx,Ny,nrcoils);
    for j = 1:nrcoils
        aux(samplingIndx) =  d_K((j-1)*M+1:j*M);
        aux2 = myifft(aux);
        d_rho = d_rho + aux2 .* conj(F.coils(:,:,j));
        d_coils(:,:,j) = myfft(aux2 .* conj(F.rho));
    end
    d_rho_and_coils = [d_rho(:);d_coils(:)];
end
'''
def parallel_MRI_adjoint(F, d_K):
    nrcoils = F["nrcoils"]
    samplingIndx = F["samplingIndx"]
    Nx = F["Nx"]
    Ny = F["Ny"]
    M = np.size(samplingIndx,1)

    aux = np.zeros((Nx,Ny))
    d_rho = np.zeros((Nx,Ny))
    d_coils = np.zeros((Nx,Ny,nrcoils))
    for j in range(nrcoils):

        aux[samplingIndx] =  d_K[j*M:(j+1)*M] #NOTE: This should work if d_k is flat
        aux2 = myifft(aux)
        d_rho = d_rho + aux2 * np.conj(F["coils"][:,:,j])
        d_coils[:,:,j] = myfft(aux2 * np.conj(F["rho"]))

    d_rho_and_coils = np.concatenate((d_rho.reshape((d_rho.shape[0],d_rho.shape[1],1)),d_coils),2) #NOTE: really flaten??
    return d_rho_and_coils

'''
function res = parallel_MRI_applyGramX(F,rho_and_coils)
   Nx = F.Nx;
   Ny = F.Ny;
   N = Nx*Ny;

   res = rho_and_coils;
   for j=1:F.nrcoils
       idx_range = j*N+1:(j+1)*N;
       res(idx_range) = rho_and_coils(idx_range).*F.Fourier_weights(:);
   end
end
'''
def parallel_MRI_applyGramX(F,rho_and_coils):
    Nx = F["Nx"]
    Ny = F["Ny"]
    N = Nx*Ny
    sp = F["Fourier_weights"].shape
    res = rho_and_coils
    #NOTE: The following should work as a work-around for the Memory issue:
    res.resize((sp[0],sp[1],rho_and_coils.size//(sp[0]*sp[1])))
    rho_and_coils.resize(res.shape)
    for j in range(F["nrcoils"]):
        #idx_range = np.arange((j+1)*N,(j+2)*N)
        #res[idx_range] = rho_and_coils[idx_range].reshape(sp) *F["Fourier_weights"]
        res[:,:,j+1] = rho_and_coils[:,:,j+1] * F["Fourier_weights"]

    return res

'''
function res = parallel_MRI_applyGramX_inv(F,rho_and_coils)
   Nx = F.Nx;
   Ny = F.Ny;
   N = Nx*Ny;

   res = rho_and_coils;
   for j=1:F.nrcoils
       idx_range = j*N+1:(j+1)*N;
       res(idx_range) = rho_and_coils(idx_range)./F.Fourier_weights(:);
   end
end
'''
def parallel_MRI_applyGramX_inv(F,rho_and_coils):
    Nx = F["Nx"]
    Ny = F["Ny"]
    N = Nx*Ny
    sp = F["Fourier_weights"].shape
    rho_and_coils.resize((sp[0],sp[1],rho_and_coils.size//(sp[0]*sp[1])))

    res = rho_and_coils
    for j in range(F["nrcoils"]):
        #idx_range = np.arange((j+1)*N,(j+2)*N)
        res[:,:,j+1] = rho_and_coils[:,:,j+1]/F["Fourier_weights"]

    return res

'''
function res = parallel_MRI_applyGramY(F,im)
   res = im;
end
'''
def parallel_MRI_applyGramY(F,im):
   res = im
   return res

"""
function F=parallel_MRI_plot(F,x_k,x_start,y_k,y_obs,k)
What = F.plotWhat;
if k>=0
    titmsg = sprintf(', it. %i',k);
else
    titmsg = sprintf(', final it. %i',abs(k));
end
if k==inf
    titmsg = ', true';
end

Nx=F.Nx; Ny=F.Ny; N=Nx*Ny;
rho_where = What.rho_where;
coils_where = What.coils_where;
if k>=-1
    if k==-1
        rho_where = rho_where+100;
        coils_where = coils_where+100;
    end
    sumcoils = myifft(sum(reshape(x_k(N+1:end),Nx,Ny,F.nrcoils),3));
    plot_rho(F,x_k,What.rho_what,rho_where,titmsg,sumcoils);
    plot_coils(F,x_k,What.coils_n1,What.coils_n2, ...
                What.coils_which, What.coils_what, coils_where,titmsg,sumcoils);
else
    sumcoils = myifft(sum(reshape(x_start(N+1:end),Nx,Ny,F.nrcoils),3));
    plot_rho(F,x_start,What.rho_what,rho_where+200,titmsg,sumcoils);
    plot_coils(F,x_start,What.coils_n1,What.coils_n2, ...
                What.coils_which, What.coils_what, coils_where+200,titmsg,sumcoils);
end
end
"""
def parallel_MRI_plot(F,x_k,x_start,y_k,y_obs,k):
    print("I should plot something!")
    return F

#%%%%%%% local functions %%%%%%%%%%%%%%%%

"""
function plot_coils(F,x,n1,n2,which,what,where,titmsg,sumcoils)
Nx = F.Nx; Ny = F.Ny;
N=Nx*Ny;
maxsum = max(max(abs(sumcoils)));
gray_aug = gray;% [0 0 1;gray; 1 0 0];

for fig = where
    figure(where(fig));
    counter=0;
    for j1=1:n1(fig)
        for j2=1:n2(fig)
            counter=counter+1;
            nr = which{fig}(counter);
            if nr>0
                c = myifft(reshape(x(nr*N+1:(nr+1)*N),Nx,Ny));%./F.Fourier_weights;
                c = c./sumcoils;
                txt1 = sprintf('coil %i',nr);
            else if (nr==0)
                    c = sumcoils;
                    txt1 = sprintf('sum coils');
                else
                    c= zeros(Nx,Ny);
                end
            end
            txt1 = [txt1,titmsg];
            thisWhat =what{fig}(counter);
            subplot(n1(fig),n2(fig),counter);
            thisWhat = what{fig}(counter);
            if (thisWhat == 'r')
                imagesc(real(c)); colorbar;
                colormap(gray_aug);
                caxis(maxsum*[-2 2]);
                title(['Re ',txt1]);
            end
            if (thisWhat == 'i')
                imagesc(imag(c)); colorbar;
                colormap(gray_aug);
                caxis(maxsum*[-2 2]);
                title(['Im ',txt1]);
            end
            if (thisWhat == 'a')
                colormap(gray_aug);
                imagesc(abs(c)); colorbar;
                caxis(maxsum*[-0.1 2]);
                title(['|.| ',txt1]);
            end
            if (thisWhat == 'p')
                colormap(hsv);
                imagesc(angle(c)); colorbar;
                title(['phase ',txt1]);
            end
            if (thisWhat == 'F')
                imagesc(log(abs(cFK))/log(10)); colorbar;
                colormap(jet);
                %caxis(cscale);
                title(['FK ',txt1]);
            end
            axis off;
        end
    end
    drawnow;
end
end
"""
#TODO: Write me!
def plot_coils(F,x,n1,n2,which,what,where,titmsg,sumcoils):
    """
    just a dummy function for now!
    """
    print("I am assumed to plot coils")

'''
function plot_rho(F,x,what,where,titmsg,sumcoils)
if (length(what)>0)
    Nx = F.Nx; Ny = F.Ny;
    N=Nx*Ny;
    rho = reshape(x(1:N),Nx,Ny) .*sumcoils;
    nsp = length(what);

    for counter = 1:nsp
        figure(where(counter)); hold off;
%        subplot(nsp,1,counter);
        switch what(counter)
            case 'r'
                imagesc(real(rho)); colormap(jet);colorbar;
                title(['real part of rho' titmsg]);
            case 'i'
                imagesc(imag(rho)); colormap(jet);colorbar;
                title(['imaginary part of rho' titmsg]);
            case 'a'
                imagesc(abs(rho));  colormap(gray); colorbar;
                title(['absolute value of rho' titmsg]);
            case 'p'
                imagesc(angle(rho)); colormap(hsv); colorbar;
                title(['phase of rho' titmsg]);
            otherwise
                error('unknown character in plot_rho');
        end
        axis off;
    end
end
drawnow;
end
'''
#TODO: Write me!
def plot_rho(F,x,what,where,titmsg,sumcoils):
    """
    just a dummy function for now
    """
    print("I am assumed to plot rho")

'''
function data = myfft(rho)
    [Nx,Ny] = size(rho);
    data = fftshift(fft2(fftshift(rho))) / sqrt(Nx*Ny);
end
'''
def myfft(rho):
    [Nx,Ny] = rho.shape
    data = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(rho))) / np.sqrt(Nx*Ny)
    return data

'''
function rho = myifft(data)
    [Nx,Ny] = size(data);
    rho = fftshift(ifft2(fftshift(data))) * sqrt(Nx*Ny);
end
'''

def myifft(data):
    [Nx,Ny] = data.shape
    rho = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(data))) * np.sqrt(Nx*Ny)
    return rho
