%% Run a sweep over values of Ro, calculate fastest growing mode

% Specify values: Ro, Ek, Nx (grid points in x), H
%       also m (vertical wavenumber)
clear
Ek=1e-3;
Nx=512; 
H=1;
m = 2*pi;  

Ro_vec = linspace(-1,-5);
mval_sweep = 0*Ro_vec;

% Loop over values of Ro
for jj = 1:length(Ro_vec)
    [mval,eigVals, ind, val, v]=INI_1d(m,Ro_vec(jj),Ek,Nx,H);
    mval_sweep(jj) = mval;
end
mval_sweep(mval_sweep<0) = 0; %mask out values with negative growth rate

figure; 
plot(Ro_vec,mval_sweep)
xlabel('Ro','Interpreter','latex')
ylabel('Fastest growth rate','Interpreter','latex')
title('$Ek = 1\times 10^{-3},\, m=2\pi$','Interpreter','latex')
set(gca, 'XDir','reverse')