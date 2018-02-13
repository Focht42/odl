import ScatOp
import ITReg.util.complete_parameter_structure
import numpy as np
import int_op.setup_iop_data
import int_op.op_S
import int_op.op_K
import int_op.farfield_matrix

""" 2 dimensional obstacle scattering problem with Dirichlet boundary condition
% see sec. 4 in T. Hohage "Logarithmic convergence rates of the iteratively
% regularized Gauss-Newton method for an inverse potential
% and an inverse scattering problem" Inverse Problems 13 (1997) 1279�1299
"""

def DirichletOp(p):
    #merge parameters with reference parameters
    F, pref = ScatOp.ScatOp(p)
    pref["op_name"] = "DirichletOp"
    # use a mixed single and double layer potential ansatz with
    # weights wSL and wDL
    pref["wSL"] = -pref["kappa"]j
    pref["wDL"] = 1
    pnew = ITReg.util.complete_parameter_structure.complete_parameter_structure(p,pref)

    for k in pnew.keys():
        F[k] = pnew[k]

    # create parameters and function handels required by regularization methods
    F["Ydim"] = 2* F["meas_directions"].shape[1] * F["inc_directions"].shape[1]
    F["evaluate"] = DirichletOp_evaluate
    F["create_synthetic_data"] = DirichletOp_create_synthetic_data
    F["derivative"] = DirichletOp_derivative
    F["adjoint"] = DirichletOp_adjoint

    return F, pnew


def DirichletOp_create_synthetic_data(F):
    bd = F["bd_ex"]["bd_eval"](F["bd_ex"],2*F["N_ieq_synth"],2)

    # compute the grid points of the exact boundary and derivatives of the
    # parametrization and save these quantities as members of bd_ex

    #set up the boudary integral operator
    Iop_data = int_op.setup_iop_data.setup_iop_data(bd,F["kappa"])
    if not F.wSL == 0:
        Iop = F["wSL"]*int_op.op_S.op_S(bd,Iop_data)
    else:
        Iop = np.zeros((bd["z"].shape[1],bd["z"].shape[1]))

    if not F.wDL == 0:
        Iop = Iop + F["wDL"]*(np.diag(bd["zpabs"])+ int_op.op_K.op_K(bd,Iop_data))

    # set up the matrix mapping the density to the far field pattern

    FF_combined = int_op.farfield_matrix.farfield_matrix(bd,F["meas_directions"],F["kappa"],F["wSL"],F["wDL"])
    farfield = []

    for l in range(1,F["inc_directions"].shape[1]+1):
        rhs = - 2*np.exp(1j*F["kappa"]*F["inc_directions"][,0].getH() *bd["z"])*bd["zpabs"]

        phi = np.linalg.solve(Iop,rhs)
        phi = phi.getH()

        complex_farfield = np.matmul(FF_combined , phi) # this is a matrix multiplication!
        farfield = np.vstack((farfield,np.real(complex_farfield),np.imag(complex_farfield)))

        noise = np.random.standard_normal(farfield.shape)

        #TODO check if i implemented this line correctly. is / a right matrix division?
        #MATLAB: data = farfield + F["noiselevel"] * noise/sqrt(noise'*F.applyGramY(F,noise));
        data = farfield + F["noiselevel"] * noise/sqrt(noise.getH()*F["applyGramY"](F,noise))

        # TODO: try catch this
        if F["plotWhat"]["field"]:
            #TODO: Translate this
            # plot_total_field(F,bd,phi,F.inc_directions(:,end),F.wSL,F.wDL);
            pass
    return data,F

#TODO: Continue
function [farfield,F]=DirichletOp_evaluate(F,coeff)
% solve the forward Dirichlet problem for the obstacle parameterized by
% coeff. Quantities needed again for the computation of derivatives and
% adjoints are stored as members of F.

F.bd.coeff = coeff;
%compute the grid points of the boundary parameterized by coeff and derivatives
%of the parametrization and save these quantities as members of F.bd
F.bd =  F.bd.bd_eval(F.bd,2*F.N_ieq,2);
Iop_data = setup_iop_data(F.bd,F.kappa);
%Iop = diag(F.bd.zpabs)+ op_K(F.bd,Iop_data) - i*F.eta*op_S(F.bd,Iop_data);
if F.wSL~=0
    Iop = F.wSL*op_S(F.bd,Iop_data);
else
    Iop = zeros(size(F.bd.z,2),size(F.bd.z,2));
end
if F.wDL~=0
	Iop = Iop + F.wDL*(diag(F.bd.zpabs)+ op_K(F.bd,Iop_data));
end
F.dudn = zeros(2*F.N_ieq,size(F.inc_directions,2));
FF_SL = farfield_matrix(F.bd,F.meas_directions,F.kappa,-1.,0.);
[F.L, F.U,F.perm] = lu(Iop,'vector');
%F.Iop=Iop;
F.FF_combined = farfield_matrix(F.bd,F.meas_directions,F.kappa,F.wSL,F.wDL);
farfield = [];

for l=1:size(F.inc_directions,2)
    rhs = 2*exp(i*F.kappa*F.inc_directions(:,l)'*F.bd.z).* ...
       (F.wDL*i*F.kappa*F.inc_directions(:,l)'*F.bd.normal +F.wSL.*F.bd.zpabs);
    %same as
    %F.dudn(:,l) = F.Iop.' \ rhs.';
    F.dudn(:,l) = (F.L.') \ ((F.U.') \ rhs(F.perm).');
    complex_farfield = FF_SL * F.dudn(:,l);
    farfield = [farfield;real(complex_farfield);imag(complex_farfield)];
end
end

function der = DirichletOp_derivative(F,h)
der = [];
for l=1:size(F.inc_directions,2)
    rhs = - 2*F.dudn(:,l) .* F.bd.der_normal(F.bd,h) .* F.bd.zpabs';
    %same as
    %phi = F.Iop \ rhs;
    phi = F.U \ (F.L \ rhs(F.perm));
    complex_farfield = F.FF_combined * phi;
    der = [der;real(complex_farfield);imag(complex_farfield)];
end
end

function adj = DirichletOp_adjoint(F,g)
res = zeros(2*F.N_ieq,1);
rhs = zeros(2*F.N_ieq,1);
N_FF = size(F.meas_directions,2);
for l=1:size(F.inc_directions,2);
    g_complex = g(2*(l-1)*N_FF+[1:N_FF]) + i*g(2*(l-1)*N_FF+[N_FF+1:2*N_FF]);
    phi = F.FF_combined'*g_complex;
    %rhs = F.Iop' \ phi;
    %rhs(F.perm) = F.L' \ (F.U' \ phi);
    rhs(F.perm) = ((phi'/F.U)/F.L)';
    res = res -2*real(rhs.*conj(F.dudn(:,l)));
end
adj = F.bd.adjoint_der_normal(F.bd,res .* F.bd.zpabs');
end




    return F,pnew
