function FFmat = farfield_matrix(bd,dir,kappa,weightSL,weightDL)
%set up matrix corresponding to the far field evaluation of a combined
%single and double layer potential with weigths weightSL and weightDL, rsp.
FFmat =  pi / (size(bd.z,2)*sqrt(8*pi*kappa)) * exp(-1i*pi/4) ...
     * (weightDL*kappa*dir'*bd.normal +i*weightSL*repmat(bd.zpabs,size(dir,2),1)) ...
     .* exp(-i*kappa* (dir' * bd.z));
end
