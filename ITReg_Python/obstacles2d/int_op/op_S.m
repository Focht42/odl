function S = op_S(bd,dat)
% set up matrix representing the single layer potential operator
% (S\phi)(z(t)):=2|z'(t)|\int_0^{2\pi} \Phi(z(t),z(s)) |z'(s)| \phi(s) ds
% Our implementation is based section 3.5 of the monograph "Inverse Acoustic
% and Electomagnetic Scattering Theory" by R. Kress and D. Colton.
% As opposed to this reference we multiplied by |z'(t)| to obtain a 
% complex symmetric matrix 

dim = size(bd.z,2);
M1 = -1/(2*pi)*real(dat.bessH0);
M = i/2*dat.bessH0;
M1_logsin = M1.* dat.logsin;
M2 = M - M1_logsin;
for j = 1:dim
    M2(j,j) = (1i/2 - dat.Euler_gamma/pi - 1/pi*log(dat.kappa/2*bd.zpabs(j)));
end;
S =  2*pi* (M1 .* dat.logsin_weights + M2/dim) .* (bd.zpabs'*bd.zpabs);
end
