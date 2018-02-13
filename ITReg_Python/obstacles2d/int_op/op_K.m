function K = op_K(bd,dat)
% set up matrix representing the double layer potential operator
% (K\phi)(z(t)):=2|z'(t)|\int_0^{2\pi}{\frac{\partial\Phi(z(t),z(s))}{\partial\nu(z(s))}
% |z'(s)| \phi(z(s)) ds
% Our implementation is based section 3.5 of the monograph "Inverse Acoustic
% and Electomagnetic Scattering Theory" by R. Kress and D. Colton.
% As opposed to this reference we multiplied by |z'(t)|.

dim = size(bd.z,2);
kappa = dat.kappa;

aux = bd.z' * bd.normal - ones(dim,2) * (bd.normal.*bd.z);
H = 0.5i*kappa^2* aux .* dat.bessH1quot;
H1 = -kappa^2/(2*pi)*aux .* real(dat.bessH1quot);
H2 = H - H1.*dat.logsin;
for j = 1:dim
    H1(j,j) = 0;
    H2(j,j) = 1/(2*pi) * (bd.normal(:,j)'*bd.zpp(:,j)) / bd.zpabs(j)^2;
end;
K =  (2*pi) * spdiags(bd.zpabs',0,dim,dim) * ( H1 .* dat.logsin_weights + H2/dim );
end
