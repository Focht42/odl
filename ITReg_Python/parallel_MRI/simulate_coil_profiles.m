function cj = simulate_coil_profiles(Nx,Ny,nrcoils)
% creates simulated coil profiles, which are Gaussians centered
% at equidistant points on a circle 
    cj = zeros(Nx,Ny,nrcoils);
    x = linspace(-1,1,Nx);
    y = linspace(-1,1,Ny);
    [Y,X] = meshgrid(y,x);
    for j=1:nrcoils
        angle = pi*(2*j-1)/nrcoils;
        x0 = sqrt(2)*cos(angle);
        y0 = sqrt(2)*sin(angle);
        cj(:,:,j) = exp(-0.25*(X-x0).^2-0.25*(Y-y0).^2);
    end
end