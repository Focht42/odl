function field = evaluate_potential(bd,phi,zz,kappa,weightSL, weightDL)
% Evaluates the combined layer potential with density phi at points zz and 
% wave length kappa.
% The weights of the single and double layer potential are given by
% the input parameters weightSL and weightDL, respectively.
Nbd =size(bd.z,2);
zz_len = size(zz,2);
field = zeros(zz_len,1);
for j=1:zz_len
    kdist = kappa * sqrt((sum( (repmat(zz(:,j),1,Nbd)-bd.z).^2,1)));
    fieldSL = (0.25*i*2*pi/Nbd)*(besselh(0,1,kdist)*(bd.zpabs'.*phi));
    fieldDL = (0.25*i*kappa^2*2*pi/Nbd) ...
         * ((zz(:,j)' * bd.normal - sum(bd.normal.*bd.z,1)) ...
        .*  (besselh(1,1,kdist)./kdist)) * phi;
    field(j)  = weightSL*fieldSL + weightDL*fieldDL;
end
end
