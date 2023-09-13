function [mval,eigVals, ind, val, v,x]=INI_1d(m,Ro,Ek,Nx,H)

%% Compute the largest growth rate and structure for the fastest growing
%% mode for the given paramaters

%% Define the domain paramaters
Ldom = 8*H;
% define x, Dxx and I
x = linspace(-Ldom,Ldom,Nx);
I = eye(Nx);
% Rigid lid matrices used for w
D2_D=ddz2(x,0);                     %2nd derivative matrix with Dirichlet BCs                          
Lmat_D = D2_D - I*m^2;

% Free slip matrices used for v
D2_N=ddz2(x,1);                     %2nd derivative matrix with Neumann BCs                                
Lmat_N = D2_N - I*m^2;

% Define Basic State
v0 = -(Ro)*exp(0.5);
V = -v0.*exp(-(x.^2)/2);
Vp = V.*(-x);
i = sqrt(-1);
Qmat = diag(Vp+1);
%figure(1); plot(x,diag(Qmat))
% figure(2); plot(x,Vp)
% figure(3); plot(x,V)
B = diag((i*m*I))./diag(Lmat_D);
 %% Build matrix
A = [[Ek*Lmat_D, i*m*Qmat]; [-diag(B), Ek*Lmat_N]];

% Solve Eigenvalue Problem Directly
[v,e]=eig(A);
eigVals=diag(e);
[val,ind]=sort(real(eigVals),'descend');

% Here we will filter out the boundary modes.
% These are modes that are affected by our choice of wall BCs in x,
%     so they have large values near the ends.
% We will look at the gravest mode (largest eigenvalue) and see if 
%     the eigenvectors are largest near the walls or not.
% If largest near the wall, we will toss the first two modes (they are
% paired).
v1 = real(v(0*length(x)+1:1*length(x),ind(1)));
v1_end1 = max(abs(v1(1:round(0.1*Nx))));
v1_end2 = max(abs(v1(end-round(0.1*Nx):end)));
v1_middle = max(abs(v1(round(0.25*Nx):round(0.75*Nx))));

if v1_middle < v1_end1 || v1_middle < v1_end2
    val = val(3:end);
    ind = ind(3:end);
    eigVals = eigVals(3:end);
end
mval = val(1);
end