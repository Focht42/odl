function T = op_T(bd,dat)
%set up matrix representing the normal derivative of the double-layer potential
%(T\phi)(z(t)):=2|z'(t)|\frac{\partial}{\partial\nu(z(t))}
%\int_0^{2\pi}{\frac{\partial\Phi(z(t),z(s))}{\partial\nu(z(s))}|z'(s)|%\phi(z(s))ds
%Our implementation is based on the paper "On the numerical solution of a
%hypersingular integral equation in scattering theory" by Rainer Kress,
%J. Comp. Appl. Math. 61:345-360, 1995, and uses the formula
%T = \frac{d}{ds}S\frac{d\phi}{ds} + \kappa^2\nu\cdot S(\nu\phi)

dim=size(bd.z,2);
z = bd.z;
zp = bd.zp;
zpp = bd.zpp;
zppp = bd.zppp;
zpabs = bd.zpabs;
kappa = dat.kappa;

N_tilde = kappa*(z'*zp -  ones(dim,1)*sum(z.*zp)) ./ dat.kdist;
N_tilde = -N_tilde'.*N_tilde;
Nker = 1i/2*N_tilde.*( kappa^2*dat.bessH0 - 2*kappa^2*dat.bessH1quot) ...
    +1i*kappa^2/2*(zp'*zp) .* dat.bessH1quot ...
    + toeplitz([pi/2 1/(4*pi)*sin(pi*[1:dim-1]/dim).^(-2)]);
N1  = -1/(2*pi)*N_tilde .* (kappa^2*real(dat.bessH0)-2*kappa^2*real(dat.bessH1quot))...
    - kappa^2/(2*pi)* (zp'*zp) .* real(dat.bessH1quot);
N2 = Nker - N1.*dat.logsin;
for j=1:dim
    N1(j,j) = -kappa^2*zpabs(j)^2/(4*pi);
    N2(j,j) = kappa^2*zpabs(j)^2/(4*pi) ...
        * ( pi*1i -1 -2*dat.Euler_gamma - 2*log(kappa*zpabs(j)/2) ) ...
        + 1/12/pi +  1/(2*pi) * sum(zp(:,j).*zpp(:,j)).^2 ./ zpabs(j).^4 ...
         - 1/(4*pi) * sum(zpp(:,j).^2)                    ./ zpabs(j).^2 ...
         - 1/(6*pi) * sum(zp(:,j).*zppp(:,j))             ./ zpabs(j).^2;
end;

T_weights = zeros(1,dim);
T_weights(2:2:dim) = (1/dim) * sin(pi*[1:2:dim-1]/dim).^(-2);
T_weights(1) = -dim/4;

T =  toeplitz(T_weights) ... 
    - 2*pi*( N1.*dat.logsin_weights + N2/dim )...
    + kappa^2.*op_S(bd,dat).* (zp'*zp) ./ (zpabs'*zpabs); 
end
