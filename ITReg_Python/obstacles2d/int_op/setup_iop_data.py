function dat = setup_iop_data(bd,kappa)
% computes data needed to set up the boundary integral matrices
% to avoid repeated computations
dim = length(bd.z);
dat.kappa = kappa;

%compute matrix of distances of grid points
t1=repmat(bd.z(1,:)',1,dim)-repmat(bd.z(1,:),dim,1);
t2=repmat(bd.z(2,:)',1,dim)-repmat(bd.z(2,:),dim,1);
dat.kdist = kappa*sqrt(t1.^2 + t2.^2);
bessj0_kdist = bessj0(dat.kdist);
dat.bessH0 = bessj0_kdist+i*bessy0(dat.kdist,bessj0_kdist);
%bessH0 = besselh(0,1,dat.kdist);

bessj1_kdist= bessj1(dat.kdist);
dat.bessH1quot = (bessj1_kdist+ i*bessy1(dat.kdist,bessj1_kdist))./dat.kdist;
%bessH1quot = besselh(1,1,dat.kdist) ./ dat.kdist;
for j=1:dim
    dat.bessH0(j,j)=1.;
end

%set up prototyp of the singularity of boundary integral operators
t=2*pi*[1:dim-1]/dim;
dat.logsin = toeplitz([1, log(4*sin(t/2).^2)]);

%quadrature weight for weight function log(4*(sin(t-tau)).^2)
sign=ones(1,dim); sign(2:2:dim)=-1;
t = 2*pi*[0:dim-1]/dim;
s=0;
for m=1:dim/2-1
    s=s+cos(m*t)/m;
end;
dat.logsin_weights = toeplitz(-2*(s + sign/dim)/dim);

%euler constant 'eulergamma'
dat.Euler_gamma =  0.577215664901532860606512;
end
