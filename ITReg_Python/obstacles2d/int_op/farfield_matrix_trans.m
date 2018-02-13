function FFmat = farfield_matrix_trans(bd,dir,kappa,weightSL,weightDL)
%set up matrix corresponding to the far field evaluation of a combined
%single and double layer potential with weigths weightSL and weightDL, rsp.
%according to eq. 4.12 by exploiting the linearity of the integral.
FFmat_a = 2*pi / (size(bd.z,2)*sqrt(8*pi*kappa)) * exp(1i*pi/4) ...
        * (-1i*weightDL*kappa*dir'*bd.normal) ...
        .* exp(-1i*kappa* (dir' * bd.z));
 
FFmat_b = 2*pi / (size(bd.z,2)*sqrt(8*pi*kappa)) * exp(1i*pi/4) ...
        * (weightSL*repmat(bd.zpabs,size(dir,2),1)) ...
        .* exp(-1i*kappa* (dir' * bd.z));
 
FFmat = [FFmat_a FFmat_b];
end
