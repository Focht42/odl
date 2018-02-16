function [h, evals] = CGNE_reg(F,y,xref,regpar,param)
% solves A h = b by CGNE with
%    A := G_X^{-1} F'* G_Y F' + regpar I
%    b := G_X^{-1}F'^* G_Y y + regpar xref
% A is self-adjoint with respect to the inner product <u,v> = u'G_X v
%
% G_X, G_Y -> F.applyGramX, F.applyGramY
% G_X^{-1} -> F.applyGramX_inv
% F'       -> F.derivative
% F'*      -> F.adjoint

verb = param.verbose;
N_CG = param.N_CG;
epsCG = 1e-2; %controls relative error of approximate solution to linear eq.

% compute rtilde = G_X b
auxY = F.applyGramY(F,y);
rtilde = F.adjoint(F,auxY);
rtilde = rtilde + regpar * F.applyGramX(F,xref);

r = F.applyGramX_inv(F,rtilde); 
d = r;
norm_r = real(rtilde'*r);
norm_r0 = norm_r;
norm_h = 0;

h = zeros(size(r));
CGstep =1;
while (sqrt(norm_r/norm_r0) > epsCG && CGstep <= N_CG)
    %(sqrt(norm_r/norm_h) > regpar*epsCG  & CGstep <= N_CG)
    % apply (G_X A) to d
    auxY = F.derivative(F,d);
    auxY = F.applyGramY(F,auxY);
    Adtilde   = F.adjoint(F,auxY);
    Adtilde   = Adtilde + regpar * F.applyGramX(F,d);

    Ada = real(Adtilde' * d);
    alpha = norm_r / real(Adtilde' * d);
    h      = h + alpha * d;
    rtilde = rtilde - alpha * Adtilde;
    r      = F.applyGramX_inv(F,rtilde);
    norm_r_old = norm_r;
    norm_r = real(rtilde'*r);
    beta   = norm_r/norm_r_old;
    d      = r + beta * d;
    norm_h = h' * F.applyGramX(F,h);
    condPrint(verb>=2,'  %i: %5.4f %5.4f;  \n',...
            CGstep,sqrt(norm_r/norm_r0),sqrt(norm_r/norm_h)/regpar);
    CGstep= CGstep + 1; 
end;
evals = 2*(CGstep-1)+1;
condPrint(isequal(CGstep,N_CG) && verb>=2, 'Maximum number of CG iterations reached.\n')
end
