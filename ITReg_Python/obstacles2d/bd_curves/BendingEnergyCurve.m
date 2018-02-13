function curve = BendingEnergyCurve(N)
curve.bd_eval = @BendingEnergyCurve_bd_eval;
curve.arc_length_der = @BendingEnergyCurve_arc_length_der;
curve.der_normal = @BendingEnergyCurve_der_normal;
curve.adjoint_der_normal = @BendingEnergyCurve_adjoint_der_normal;
curve.applyHess = @BendingEnergyCurve_applyHess;
curve.applyHess_inv = @BendingEnergyCurve_applyHess_inv;
curve.update =  @BendingEnergyCurve_update;
end

function curve = BendingEnergyCurve_bd_eval(curve,coeff,der)
  curve.z = 0;  %%%% TODO
  curve.zp = 0;
  curve.zpabs = sqrt(curve.zp(1,:).^2 + curve.zp(2,:).^2);
  curve.normal = [curve.zp(2,:);
                   -curve.zp(1,:)];
  if der>=2
      curve.zpp = 0;
  end
  if der>=3
      curve.zppp = 0;
  end
  curve.DerMatrix = 0;  %%%% TODO
  curve.Hess = 0;  %%%% TODO
end

function res = BendingEnergyCurve_arc_length_der(curve,val)
%%% NOT TESTED!
N = length(val);
if mod(N,2)==1
    error('length of t should be even')
end
coeffhat = fftshift(fft(val))
dhds = ifft(fftshift((i*[-n/2:n/2-1]') .*coeffhat./curve.zpabs';
end

function res = BendingEnergyCurve_der_normal(curve,h)
  res = curve.DerMatrix * h;
end

function res = BendingEnergyCurve_der_normal_adjoint(curve,h)
  res = curve.DerMatrix' * h;
end

function res = BendingEnergyCurve_applyHess(curve,h)
  res = curve.Hess * h;
end

function res = BendingEnergyCurve_applyHess_inv(curve,h)
  res = curve.Hess\ h;
end

function coeff = BendingEnergyCurve_update(curve,direcion)
%%%% TODO
end
