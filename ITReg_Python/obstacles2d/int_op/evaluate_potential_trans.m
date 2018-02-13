function field = evaluate_potential_trans(bd,phi,zz,kappa,weightSL, weightDL)
% Evaluates weightDL * double layer potential of the first half of bd.phi 
% + weightSL * single layer potential of the second half of bd.phi

Nbd =size(bd.z,2);
zz_len = size(zz,2);
field = zeros(zz_len,1);
nphi = size(phi,1)/2;
for j=1:zz_len
    kdist = kappa * sqrt((sum( (repmat(zz(:,j),1,Nbd)-bd.z).^2,1)));
    fieldSL = (0.25*1i*2*pi/Nbd)*(besselh(0,1,kdist)*(bd.zpabs'.*phi((nphi+1):2*nphi,:)));
    fieldDL = (0.25*1i*kappa^2*2*pi/Nbd) ...
         * ((zz(:,j)' * bd.normal - sum(bd.normal.*bd.z,1)) ...
        .*  (besselh(1,1,kdist)./kdist)) * phi(1:nphi,:);
    field(j)  = weightSL*fieldSL + weightDL*fieldDL;
end
end
