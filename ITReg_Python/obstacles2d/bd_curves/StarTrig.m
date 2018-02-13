%% The class StarTrig describes boundaries of domains in R^2 which are 
%% star-shaped w.r.t. the origin and which are parameterized by 
%%      z(t) = q(t)[cos(t);sin(t)]      0<=t<=2pi
%% where q is a trigonometric polynomial with N coefficient. 
%% Here N must be even, so the highest order monomial is cos(t*N/2), 
%% but sin(t*N/2) does not occur. 
%% z and its derivatives are sampled at n equidistant points. 
%% Application of the Gramian matrix and its inverse w.r.t. the 
%% Sobolev norm ||q||_{H^s} are implemented. 

function curve = StarTrig(N,s)
curve.type = 'StarTrig';
if mod(N,2)==1
    error('N should be even');
end
curve.coeff = zeros(N,1);
curve.sobo_index = s;

curve.compute_FK = @StarTrig_compute_FK;
curve.bd_eval = @StarTrig_bd_eval;
curve.der_normal = @StarTrig_der_normal;
curve.adjoint_der_normal    = @StarTrig_adjoint_der_normal;
curve.arc_length_der = @StarTrig_arc_length_der;
curve.applyGram = @StarTrig_applyGram;
curve.applyGram_inv = @StarTrig_applyGram_inv;
curve.L2err = @StarTrig_L2err;
curve.coeff2Curve = @StarTrig_coeff2Curve;
%%%% The following two methods are not needed for operators depending
%%%% only on the curve, but not on its parametrization.
%%%% They are included for test purposes.
curve.derivative = @StarTrig_derivative;
curve.adjoint_derivative    = @StarTrig_adjoint_derivative;
end

function coeffhat2= StarTrig_compute_FK(val,n)
% computes n Fourier coeffients to the point values given by by val
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

function curve = StarTrig_radial(curve,n,der)
% evaluates all derivatives of the radial function up to order der
% at n equidistant points

coeffhat = curve.compute_FK(curve.coeff,n);
for d=0:der
    curve.q(d+1,:) = real(ifft(fftshift( (i*[-n/2:n/2-1]').^d .* coeffhat)))';
end
end

function curve = StarTrig_bd_eval(curve,n,der)
curve = StarTrig_radial(curve,n,der);
q=curve.q;
t = 2*pi*[0:n-1]/n;
cost = cos(t);
sint = sin(t);

curve.z = [q(1,:).*cost;
    q(1,:).*sint];
if der>=1
    curve.zp = [q(2,:).*cost - q(1,:).*sint;
        q(2,:).*sint + q(1,:).*cost];
    curve.zpabs = sqrt(curve.zp(1,:).^2 + curve.zp(2,:).^2);
    %outer normal vector
    curve.normal = [curve.zp(2,:);
        -curve.zp(1,:)];
end
if der>=2
    curve.zpp = [q(3,:).*cost - 2*q(2,:).*sint - q(1,:).*cost;
        q(3,:).*sint + 2*q(2,:).*cost - q(1,:).*sint];
end
if der>=3
    curve.zppp = [q(4,:).*cost - 3*q(3,:).*sint - 3*q(2,:).*cost + q(1,:).*sint;
        q(4,:).*sint + 3*q(3,:).*cost - 3*q(2,:).*sint - q(1,:).*cost];
end
if der>3
    error('only derivatives up to order 3 implemented');
end
end

function der = StarTrig_der_normal(curve,h)
%computes the normal part of the perturbation of the curve caused by
%perturbing the coefficient vector curve.coeff in direction h

n=size(curve.q,2);
h_long = ifft(fftshift(curve.compute_FK(h,n)));
der =  curve.q(1,:)' .* h_long ./ curve.zpabs';
end

function adj = StarTrig_adjoint_der_normal(curve,g)
%applies the adjoint of the linear mapping h->der_normal(curve,h) to g

N = length(curve.coeff);
n = length(g);
adj_long = g.*curve.q(1,:)' ./ curve.zpabs';
adj = ifft(fftshift(curve.compute_FK(adj_long,N))) * n/N;
end

function der = StarTrig_derivative(curve,h)
%computes the perturbation of the curve caused by perturbing the coefficient
% vector curve.coeff in direction h
n=size(curve.q,2);
h_long = ifft(fftshift(curve.compute_FK(h,n)))';
t = 2*pi*[0:n-1]/n;
der =  [h_long.*cos(t);h_long.*sin(t)];
end

function adj = StarTrig_adjoint_derivative(curve,g)
%applies the adjoint of the linear mapping h->derivative(curve,h) to g
N = length(curve.coeff);
n = length(g);
t = 2*pi*[0:n-1]/n;
adj_long = g(1,:).*cos(t)+g(2,:).*sin(t);
adj = ifft(fftshift(compute_FK(adj_long,N))) * n/N;
end

function dhds = StarTrig_arc_length_der(curve,h)
% computes the derivative of h with respect to arclength
n=size(curve.q,2);
dhds = ifft(fftshift((i*[-n/2:n/2-1]') .*compute_FK(h,n)))./curve.zpabs';
end

function res = StarTrig_applyGram(curve,v)
n=length(v);
vhat = fftshift(fft(v));
res = real(ifft(fftshift((1+[-n/2:n/2-1]'.^2).^curve.sobo_index.*vhat)));
end

function res = StarTrig_applyGram_inv(curve,v)
n=length(v);
vhat = fftshift(fft(v));
res = real(ifft(fftshift((1+[-n/2:n/2-1]'.^2).^(-curve.sobo_index).*vhat)));
end

function res = StarTrig_L2err(q1,q2)
res = norm(q1-q2)/sqrt(length(q1));
end

function pts = StarTrig_coeff2Curve(coeff,n)
radial = ifft(fftshift(StarTrig_compute_FK(coeff,n)));
t = 2*pi/n * [0:n-1];
pts = [radial'.*cos(t);
    radial'.*sin(t)];
end