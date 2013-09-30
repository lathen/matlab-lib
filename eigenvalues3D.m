function [eig1 eig2 eig3] = eigenvalues3D(t1, t2, t3, t4, t5, t6)
% Parallel computation of eigenvalues in 3D. This follows appendix G in
% Gunnar Farnebäck's PhD thesis "Polynomial Expansion for Orientation and
% Motion Estimation"
%
% Elements are passed as arguments t1-t6 according to the layout:
%     | t1 t2 t3 |
% T = | t2 t4 t5 | 
%     | t3 t5 t6 |
%
% Author: Gunnar Farnebäck
%	  Computer Vision Laboratory
%	  Linköping University, Sweden
%	  gf@isy.liu.se

tr = (t1 + t4 + t6) / 3;

a = t1 - tr;
b = t4 - tr;
c = t6 - tr;

d = t2;
e = t3;
f = t5;

p = a.*b + a.*c + b.*c - d.^2 - e.^2 - f.^2;
q = a.*f.^2 + b.*e.^2 + c.*d.^2 - 2*d.*e.*f - a.*b.*c;
p(p>0) = 0;
beta = sqrt(-4/3*p);
phi = p.*beta;
phi(phi==0) = 1;                       % To avoid numerical problems.
gamma = 3*q./phi;
gamma(gamma>1) = 1;                    % To avoid numerical problems.
gamma(gamma<-1) = -1;                  % To avoid numerical problems.
alpha = acos(gamma)/3;

eig1 = tr + beta.*cos(alpha);
eig2 = tr + beta.*cos(alpha-2*pi/3);
eig3 = tr + beta.*cos(alpha+2*pi/3);
