%% The class GenTrig describes boundaries of domains in R^2 which are 
%% parameterized by 
%%      z(t) = [z_1(t),z_2(t)]      0<=t<=2pi
%% where z_1 and z_2 are trigonometric polynomials with N coefficient. 
%% Here N must be even, so the highest order monomial is cos(t*N/2), 
%% but sin(t*N/2) does not occur. 
%% z and its derivatives are sampled at n equidistant points. 
%% Application of the Gramian matrix and its inverse w.r.t. the 
%% Sobolev norm ||z||_{H^s} are implemented. 

function curve = GenTrig(N,s)
if mod(N,2)==1
    error('N should be even');
end
curve.type = 'GenTrig';
curve.coeff = zeros(2*N,1);
curve.sobo_index = s;

curve.compute_FK = @GenTrig_compute_FK;
curve.bd_eval = @GenTrig_bd_eval;
curve.der_normal = @GenTrig_der_normal;
curve.adjoint_der_normal    = @GenTrig_adjoint_der_normal;
curve.arc_length_der = @GenTrig_arc_length_der;
curve.applyGram = @GenTrig_applyGram;
curve.applyGram_inv = @GenTrig_applyGram_inv;
curve.L2err = @GenTrig_L2err;
curve.coeff2Curve = @GenTrig_coeff2Curve;
end

function coeffhat2= GenTrig_compute_FK(val,n)
% computes n Fourier coeffients to the point values given by by val.
% such that ifft(fftshift(coeffhat)) is an interpolation of val

if mod(n,2)==1
    error('length of t should be even')
end
N = length(val);
coeffhat = fft(val);
coeffhat2 = zeros(n,1);
if (n>=N)
    coeffhat2(1:N/2)= coeffhat(1:N/2);
    coeffhat2(n-N/2+2:n) = coeffhat(N/2+2:N);
    if (n>N)
        coeffhat2(N/2+1) = 0.5*coeffhat(N/2+1);
        coeffhat2(n-N/2+1) = 0.5*coeffhat(N/2+1);
    else %n==N
        coeffhat2(N/2+1) = coeffhat(N/2+1);
    end
else
    coeffhat2(1:n/2) = coeffhat(1:n/2);
    coeffhat2(n/2+2:n) = coeffhat(N-n/2+2:N);
    coeffhat2(n/2+1) = 0.5*(coeffhat(n/2+1)+coeffhat(N-n/2+1));
end;
coeffhat2 = n/N*fftshift(coeffhat2);
end

function curve = GenTrig_bd_eval(curve,n,der)
% evaluates the first der derivatives of the parametrization of
% the curve on n equidistant time points
N = length(curve.coeff)/2;

coeffhat = [curve.compute_FK(curve.coeff(1:N),n),...
    curve.compute_FK(curve.coeff(N+1:2*N),n)].';
curve.z = [real(ifft(fftshift( coeffhat(1,:))));...
    real(ifft(fftshift( coeffhat(2,:))))];
if der>=1
    curve.zp = [real(ifft(fftshift( (i*[-n/2:n/2-1]) .* coeffhat(1,:))));...
        real(ifft(fftshift( (i*[-n/2:n/2-1]) .* coeffhat(2,:))))];
    curve.zpabs = sqrt(curve.zp(1,:).^2 + curve.zp(2,:).^2);
    %outer normal vector
    curve.normal = [curve.zp(2,:);...
        -curve.zp(1,:)];
end
if der>=2
    curve.zpp = [real(ifft(fftshift( (i*[-n/2:n/2-1]).^2 .* coeffhat(1,:))));...
        real(ifft(fftshift( (i*[-n/2:n/2-1]).^2 .* coeffhat(2,:))))];
end
if der>=3
    curve.zppp = [real(ifft(fftshift( (i*[-n/2:n/2-1]).^3 .* coeffhat(1,:))));...
        real(ifft(fftshift( (i*[-n/2:n/2-1]).^3 .* coeffhat(2,:))))];
end
if der>3
    error('only derivatives up to order 3 implemented');
end
end

function der = GenTrig_der_normal(curve,h)
%computes the normal part of the perturbation of the curve caused by
%perturbing the coefficient vector curve.coeff in direction h
N= length(h)/2;
n = size(curve.z,2);
if N==n
    hn= [h(1:n)';h(n+1:2*n)'];
else
    h_hat = [curve.compute_FK(h(1:N),n),...
        curve.compute_FK(h(N+1:2*N),n)].';
    hn = [real(ifft(fftshift(h_hat(1,:))));...
        real(ifft(fftshift(h_hat(2,:))))];
end
der = (sum(hn.*curve.normal,1)./curve.zpabs)';
end

function adj = GenTrig_adjoint_der_normal(curve,g)
%applies the adjoint of the linear mapping h->der_normal(curve,h) to g
N = length(curve.coeff)/2;
n= length(g);
adj_n = repmat(g'./curve.zpabs,2,1) .* curve.normal;
if N==n
    adj = [adj_n(1,:)';adj_n(2,:)'];
else
    adj_hat = [curve.compute_FK(adj_n(1,:),N),...
        curve.compute_FK(adj_n(2,:),N)]*n/N;
    adj = [ifft(fftshift(adj_hat(:,1)));...
        ifft(fftshift(adj_hat(:,2)))];
end
end

function dhds = GenTrig_arc_length_der(curve,h)
% computes the derivative of h with respect to arclength
n=size(curve.q,2);
dhds = ifft(fftshift((i*[-n/2:n/2-1]') .*curve.compute_FK(h,n)))./curve.zpabs';
end

function res = GenTrig_applyGram(curve,v)
n=length(v)/2;
vhatx = fftshift(fft(v(1:n)));
vhaty = fftshift(fft(v(n+1:2*n)));
res = [real(ifft(fftshift((1+[-n/2:n/2-1]'.^2).^curve.sobo_index.*vhatx)));...
    real(ifft(fftshift((1+[-n/2:n/2-1]'.^2).^curve.sobo_index.*vhaty)))];
end

function res = GenTrig_applyGram_inv(curve,v)
n=length(v)/2;
vhatx = fftshift(fft(v(1:n)));
vhaty = fftshift(fft(v(n+1:2*n)));
res = [real(ifft(fftshift((1+[-n/2:n/2-1]'.^2).^(-curve.sobo_index).*vhatx)));...
    real(ifft(fftshift((1+[-n/2:n/2-1]'.^2).^(-curve.sobo_index).*vhaty)))];
end

function res = GenTrig_L2err(q1,q2)
res = norm(q1-q2)/sqrt(length(q1));
end

function pts = GenTrig_coeff2Curve(coeff,n)
N = length(coeff)/2;
coeffhat = [GenTrig_compute_FK(coeff(1:N),N),...
    GenTrig_compute_FK(coeff(N+1:2*N),N)].';
pts = [real(ifft(fftshift( coeffhat(1,:))));...
    real(ifft(fftshift( coeffhat(2,:))))];
end
