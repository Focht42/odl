% test operator for boundary representations 
% F(q) := P\circ z_q where the zero level set of P:R^2\to R 
% is peanut shaped

function [F,p]=PeanutFit(p)
%% merge parameters with reference parameters
F.op_name = 'PeanutFit';
F.syntheticdata_flag = true;
F.alpha = 0.25;
F.scaling = 2;
F.sobo_index =  0.6;  % Sobolev index of preimage space
F.N_FK = 32;     % dimension of the trigonometric subspace describing curves
F.init_guess = 0.7*ones(F.N_FK,1);   %initial guess is a circle
F.true_curve = 'peanut';
F.plotWhat.plots = 'X';
if isfield(p,'plotWhat')
    F.plotWhat = complete_parameter_structure(p.plotWhat,F.plotWhat);
end
F.noiselevel = 0.0;
F = complete_parameter_structure(p,F); 
p=F;

%% compute precomputable quantities
% initialize boundary curve structure
F.bd = star_trig(p.N_FK,p.sobo_index);
if isfield(p,'true_curve')
F.bd_exact = star_curve(p.true_curve);
  F.xdag = F.bd_exact.radial(F.bd_exact,p.N_FK)';
    N_plot =128;
    t_plot = 2*pi/N_plot*[0:N_plot];
    rad_plot =  F.bd_exact.radial(F.bd_exact,N_plot);
    rad_plot = [rad_plot rad_plot(1)];
    F.xdag_pts_plot = [rad_plot.*cos(t_plot);
        rad_plot.*sin(t_plot)];
end
               
F.Xdim = F.N_FK;
F.Ydim = F.N_FK;
F.evaluate = @PeanutFit_evaluate;
F.create_synthetic_data = @PeanutFit_create_synthetic_data;
F.derivative = @PeanutFit_derivative;
F.adjoint = @PeanutFit_adjoint;
F.plot = @PeanutFit_plot;
F.applyGramX = @PeanutFit_applyGramX;
F.applyGramX_inv = @PeanutFit_applyGramX_inv;
F.applyGramY = @PeanutFit_applyGramY;
F.other_X_err = {@PeanutFit_L2err;'L2err'};
end

function [data,F]=PeanutFit_create_synthetic_data(F)
noise = randn(F.N_FK,1);
data = F.noiselevel * noise/sqrt(noise'*F.applyGramY(F,noise));
end

function [VAL,F]=PeanutFit_evaluate(F,coeff)
F.bd = F.bd;
F.bd.coeff = coeff;
F.bd = F.bd.bd_eval(F.bd,F.Ydim,1);

% X = F.bd.z(1,:);
% Y = F.bd.z(2,:);
% rad = sqrt(X.^2+0.25*Y.^2);
% ind = find(rad==0); rad(ind)= 1e-10;
% P = 2*(X./rad).^2-1;
% c = 0.25/(1+sqrt(F.alpha+1));
% R = c*(P+sqrt(F.alpha+P.^2));
% VAL = (rad.^2/(F.scaling^2)-R)';
% Py = -1*X.^2.*Y./(rad.^4);
% Px = 4*X./(rad.^2).*(1-X.^2./rad.^2);
% F.grad_VAL = [(2./F.scaling^2)*X- c*(1+P/sqrt(F.alpha+P.^2)).*Px;...
%     (2/F.scaling^2)*Y-c*(1+P/sqrt(F.alpha+P.^2)).*Py];

X = F.bd.z(1,:);
Z = F.bd.z(2,:);
rad = sqrt(X.^2+Z.^2);
ind = find(rad==0); rad(ind)= 1e-10;
P = 2*(Z./rad).^2-1;
c = 0.25/(1+sqrt(F.alpha+1));
R = c*(P+sqrt(F.alpha+P.^2));
VAL = (rad.^2/(F.scaling^2)-R)';
Px = -4*Z.^2.*X./(rad.^4);
Pz = 4*Z./(rad.^2).*(1-Z.^2./rad.^2);
F.grad_VAL = [(2/F.scaling^2)*X-c*(1+P./sqrt(F.alpha+P.^2)).*Px;...
    (2./F.scaling^2)*Z- c*(1+P./sqrt(F.alpha+P.^2)).*Pz];
end

function der=PeanutFit_derivative(F,h_coeff)
h = F.bd.derivative(F.bd,h_coeff);
der = sum(h.*F.grad_VAL,1)';
end

function adj=PeanutFit_adjoint(F,g)
adj = (g*[1 1])'.*F.grad_VAL;
adj = F.bd.adjoint_derivative(F.bd,adj);
end

function F = PeanutFit_plot(F,x_k,x_start,y_k,y_obs,k)
nr_plots = length(F.plotWhat.plots);
if isfield(F,'fig_rec')
    figure(F.fig_rec);
else
    scrsz = get(0,'ScreenSize');
    fac=0.5;
    F.fig_rec = figure('Name',['peanut fit'],'NumberTitle','off',...
        'Position',[1 (0.9-fac)*scrsz(4)  fac*length(F.plotWhat.plots)*scrsz(4) fac*scrsz(4)])
end
plot_count = 0;
for this_plot = F.plotWhat.plots
    plot_count = plot_count+1;
    switch this_plot
        case 'X'
            subplot(1,nr_plots,plot_count);
            bd=F.bd;
            pts_xk = bd.coeff2Curve(x_k,128);
            pts_xk = [pts_xk, pts_xk(:,1)];
            if isfield(F,'rec_plot_handle')
                delete(F.rec_plot_handle);
                F.rec_plot_handle = plot(pts_xk(1,:),pts_xk(2,:),'b-');
            else
                pts_start = bd.coeff2Curve(x_start,128);
                pts_start = [pts_start, pts_start(:,1)];
                %subplot(1,2,1)
                hold off;
                F.rec_plot_handle = plot(pts_xk(1,:),pts_xk(2,:),'b-');
                hold on
                plot(F.xdag_pts_plot(1,:), F.xdag_pts_plot(2,:),'r-');
                plot(pts_start(1,:),pts_start(2,:),'k-');
                legend('reconstrucion','exact solution','initial guess','measurement points');
            end
            if k>=0
                tit = sprintf('source reconstruction at step %i',k);
            else
                tit = sprintf('final source reconstruction at step %i',abs(k));
            end
            title(tit);
            drawnow;
            pause(1)
         otherwise
            error('Unknown plot command');
    end
end
end

function res = PeanutFit_applyGramX(F,v)
  res = F.bd.applyGram(F.bd,v);
end

function res = PeanutFit_applyGramX_inv(F,v)
  res = F.bd.applyGram_inv(F.bd,v);
end

function res = PeanutFit_L2err(F,h)
res = F.bd.L2err(h,F.xdag);
end

function res = PeanutFit_applyGramY(F,g)
  res = (2*pi/F.Ydim) * g;
end
