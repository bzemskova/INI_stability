%% Test and plot for a particular set of values

% Specify values: Ro, Ek, Nx (grid points in x), H
%       also m (vertical wavenumber)
clear
Ro=-1.2;
Ek=1e-3;
Nx=512; 
H=1;
m = 2*pi; 

% Run the eigenvalue calculation
[mval,eigVals, ind, val, v,x]=INI_1d(m,Ro,Ek,Nx,H);

% Plot vertical velocity
z = linspace(-1,0);
psi1 = v(1*length(x)+1:2*length(x),ind(1)); % this is streamfunction \psi of the fastest growing mode
psix = gradient(psi1,x);  % d\psi/dx

w1 = exp(1i*m*z); % calculate vertical modes
w = 1i*psix/m; % calculate w = i\psi/m (from continuity equation)
% Now compute 2D vertical velocity field w(x,z)
for j = 1:length(w)
    W(j,:) = w(j)*w1;
end

figure; pcolor(x,z,imag(W)'); shading interp
colorbar
clim([-max(imag(W(:))) max(imag(W(:)))])
xlabel('x','Interpreter','latex')
ylabel('z','Interpreter','latex')
title('$w: Ek = 1\times 10^{-3},\, m=2\pi$','Interpreter','latex')

