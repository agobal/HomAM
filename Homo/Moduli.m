% E = 150000000000; % Young's modulus
% nu = 0.3;         % Poisson's ratio

% G = E/(2*(1 + nu)); % Shear modulus
Ex = 1.461;
nunu = 0.13;
matmat = [1/Ex -nunu/Ex -nunu/Ex 0 0 0;
    -nunu/Ex 1/Ex -nunu/Ex 0 0 0;
    -nunu/Ex -nunu/Ex 1/Ex 0 0 0;
    0 0 0 (1+nunu)/Ex 0 0;
    0 0 0 0 (1+nunu)/Ex 0;
    0 0 0 0 0 (1+nunu)/Ex];
matmat = inv(matmat);
C11 = matmat(1, 1);
C12 = matmat(1, 2);
C44 = matmat(6, 6);

c = zeros(3, 3, 3, 3);
c(1, 1, 1, 1) = C11; c(2, 2, 2, 2) = C11; c(3, 3, 3, 3) = C11;
c(1, 1, 2, 2) = C12; c(1, 1, 3, 3) = C12; c(2, 2, 1, 1) = C12;
c(2, 2, 3, 3) = C12; c(3, 3, 1, 1) = C12; c(3, 3, 2, 2) = C12;
c(1, 2, 1, 2) = C44; c(2, 1, 2, 1) = C44; c(1, 3, 1, 3) = C44; 
c(3, 1, 3, 1) = C44; c(2, 3, 2, 3) = C44; c(3, 2, 3, 2) = C44;
c(1, 2, 2, 1) = C44; c(2, 1, 1, 2) = C44; c(1, 3, 3, 1) = C44;
c(3, 1, 1, 3) = C44; c(2, 3, 3, 2) = C44; c(3, 2, 2, 3) = C44;

c2 = zeros(3, 3, 3, 3);
Ex = 0.7406;
nunu = 0.1;
matmat = [1/Ex -nunu/Ex -nunu/Ex 0 0 0;
    -nunu/Ex 1/Ex -nunu/Ex 0 0 0;
    -nunu/Ex -nunu/Ex 1/Ex 0 0 0;
    0 0 0 (1+nunu)/Ex 0 0;
    0 0 0 0 (1+nunu)/Ex 0;
    0 0 0 0 0 (1+nunu)/Ex];
matmat = inv(matmat);
c2_11 = matmat(1, 1);
c2_12 = matmat(1, 2);
c2_44 = matmat(6, 6);
c2(1, 1, 1, 1) = c2_11; c2(2, 2, 2, 2) = c2_11; c2(3, 3, 3, 3) = c2_11;
c2(1, 1, 2, 2) = c2_12; c2(1, 1, 3, 3) = c2_12; c2(2, 2, 1, 1) = c2_12;
c2(2, 2, 3, 3) = c2_12; c2(3, 3, 1, 1) = c2_12; c2(3, 3, 2, 2) = c2_12;
c2(1, 2, 1, 2) = c2_44; c2(2, 1, 2, 1) = c2_44; c2(1, 3, 1, 3) = c2_44; 
c2(3, 1, 3, 1) = c2_44; c2(2, 3, 2, 3) = c2_44; c2(3, 2, 3, 2) = c2_44;
c2(1, 2, 2, 1) = c2_44; c2(2, 1, 1, 2) = c2_44; c2(1, 3, 3, 1) = c2_44;
c2(3, 1, 1, 3) = c2_44; c2(2, 3, 3, 2) = c2_44; c2(3, 2, 2, 3) = c2_44;

C = c;
c = c2;
c0 = C;

%% Finding G and its Fourier transform

delta_c = c - c0;
DELTAC = [delta_c(1, 1, 1, 1) delta_c(1, 1, 2, 2) delta_c(1, 1, 3, 3) delta_c(1, 1, 1, 2) delta_c(1, 1, 1, 3) delta_c(1, 1, 2, 3);
    delta_c(2, 2, 1, 1) delta_c(2, 2, 2, 2) delta_c(2, 2, 3, 3) delta_c(2, 2, 1, 2) delta_c(2, 2, 1, 3) delta_c(2, 2, 2, 3);
    delta_c(3, 3, 1, 1) delta_c(3, 3, 2, 2) delta_c(3, 3, 3, 3) delta_c(3, 3, 1, 2) delta_c(3, 3, 1, 3) delta_c(3, 3, 2, 3);
    delta_c(1, 2, 1, 1) delta_c(1, 2, 2, 2) delta_c(1, 2, 3, 3) delta_c(1, 2, 1, 2) delta_c(1, 2, 1, 3) delta_c(1, 2, 2, 3);
    delta_c(1, 3, 1, 1) delta_c(1, 3, 2, 2) delta_c(1, 3, 3, 3) delta_c(1, 3, 1, 2) delta_c(1, 3, 1, 3) delta_c(1, 3, 2, 3);
    delta_c(2, 3, 1, 1) delta_c(2, 3, 2, 2) delta_c(2, 3, 3, 3) delta_c(2, 3, 1, 2) delta_c(2, 3, 1, 3) delta_c(2, 3, 2, 3)];

CBIG = [c(1, 1, 1, 1) c(1, 1, 2, 2) c(1, 1, 3, 3) c(1, 1, 1, 2) c(1, 1, 1, 3) c(1, 1, 2, 3);
    c(2, 2, 1, 1) c(2, 2, 2, 2) c(2, 2, 3, 3) c(2, 2, 1, 2) c(2, 2, 1, 3) c(2, 2, 2, 3);
    c(3, 3, 1, 1) c(3, 3, 2, 2) c(3, 3, 3, 3) c(3, 3, 1, 2) c(3, 3, 1, 3) c(3, 3, 2, 3);
    c(1, 2, 1, 1) c(1, 2, 2, 2) c(1, 2, 3, 3) c(1, 2, 1, 2) c(1, 2, 1, 3) c(1, 2, 2, 3);
    c(1, 3, 1, 1) c(1, 3, 2, 2) c(1, 3, 3, 3) c(1, 3, 1, 2) c(1, 3, 1, 3) c(1, 3, 2, 3);
    c(2, 3, 1, 1) c(2, 3, 2, 2) c(2, 3, 3, 3) c(2, 3, 1, 2) c(2, 3, 1, 3) c(2, 3, 2, 3)];

CZERO = [c0(1, 1, 1, 1) c0(1, 1, 2, 2) c0(1, 1, 3, 3) c0(1, 1, 1, 2) c0(1, 1, 1, 3) c0(1, 1, 2, 3);
    c0(2, 2, 1, 1) c0(2, 2, 2, 2) c0(2, 2, 3, 3) c0(2, 2, 1, 2) c0(2, 2, 1, 3) c0(2, 2, 2, 3);
    c0(3, 3, 1, 1) c0(3, 3, 2, 2) c0(3, 3, 3, 3) c0(3, 3, 1, 2) c0(3, 3, 1, 3) c0(3, 3, 2, 3);
    c0(1, 2, 1, 1) c0(1, 2, 2, 2) c0(1, 2, 3, 3) c0(1, 2, 1, 2) c0(1, 2, 1, 3) c0(1, 2, 2, 3);
    c0(1, 3, 1, 1) c0(1, 3, 2, 2) c0(1, 3, 3, 3) c0(1, 3, 1, 2) c0(1, 3, 1, 3) c0(1, 3, 2, 3);
    c0(2, 3, 1, 1) c0(2, 3, 2, 2) c0(2, 3, 3, 3) c0(2, 3, 1, 2) c0(2, 3, 1, 3) c0(2, 3, 2, 3)];



% For all the query points inside the structure
N = 27;
X = [1 1 1; 2 1 1; 3 1 1; 1 2 1; 2 2 1; 3 2 1; 1 3 1; 2 3 1; 3 3 1;
     1 1 2; 2 1 2; 3 1 2; 1 2 2; 2 2 2; 3 2 2; 1 3 2; 2 3 2; 3 3 2;
     1 1 3; 2 1 3; 3 1 3; 1 2 3; 2 2 3; 3 2 3; 1 3 3; 2 3 3; 3 3 3];
 for i = 1:3
     for j = 1:27
         X(j, i) = (X(j, i) - 1)*10 + 5;
     end
 end

 Comp = zeros(6*N, 6*N);
 Comp0 = zeros(6*N, 6*N);
 
 B = zeros(6*N, 6*N);
 DC = zeros(6*N, 6*N);
 for I = 1:N
     for J = 1:N
         if (I ~= J)
             x1 = X(I, :); % Location at which we want G
             x2 = X(J, :);
             x = x1 - x2;
             
             v1 = (cross(x, [1 1 1]));
             if v1 == 0
                 v1 = (cross(x,[0 1 0]));
             end
             v1 = v1/norm(v1);
             v2 = (cross(x, v1));
             v2 = v2/norm(v2);
             
             G_theta = zeros(3, 3);
             nint = 40;
             theta = [0:nint - 1]/nint*(2*pi);
             dtheta = (2*pi)/nint;
             
             for m = 1:nint
                 z = v1*cos(theta(m)) + v2*sin(theta(m));
                 
                 zz = zeros(3, 3);
                 for i = 1:3
                     for j = 1:3
                         for k = 1:3
                             for l = 1:3
                                 zz(j, k) = c(i, j, k, l)*z(i)*z(l) + zz(j, k);
                             end
                         end
                     end
                 end
                 zzinv = inv(zz);
                 
                 G_theta = zzinv + G_theta;
             end
             
             G = G_theta/8/pi^2/norm(x)*dtheta;
             
             Gamma_theta = zeros(3, 3, 3, 3);
             for mm = 1:nint
                 z = v1*cos(theta(mm)) + v2*sin(theta(mm));
                 
                 zz1 = zeros(3, 3);
                 for p = 1:3
                     for j = 1:3
                         for m = 1:3
                             for l = 1:3
                                 zz1(m, p) = c(p, j, m, l)*z(j)*z(l) + zz1(m, p);
                             end
                         end
                     end
                 end
                 zz1inv = inv(zz1);
                 
                 zz2 = zeros(3, 3);
                 for p = 1:3
                     for j = 1:3
                         for n = 1:3
                             for l = 1:3
                                 zz2(n, p) = c(p, j, n, l)*z(j)*z(l) + zz2(n, p);
                             end
                         end
                     end
                 end
                 zz2inv = inv(zz2);
                 
                 zz3 = zeros(3, 3);
                 for q = 1:3
                     for j = 1:3
                         for m = 1:3
                             for l = 1:3
                                 zz3(m, q) = c(q, j, m, l)*z(j)*z(l) + zz3(m, q);
                             end
                         end
                     end
                 end
                 zz3inv = inv(zz3);
                 
                 zz4 = zeros(3, 3);
                 for q = 1:3
                     for j = 1:3
                         for n = 1:3
                             for l = 1:3
                                 zz4(n, q) = c(q, j, n, l)*z(j)*z(l) + zz4(n, q);
                             end
                         end
                     end
                 end
                 zz4inv = inv(zz4);
                 
                 gamma = zeros(3, 3, 3, 3);
                 for m = 1:3
                     for n = 1:3
                         for p = 1:3
                             for q = 1:3
                                 gamma(m, n, p, q) = zz1inv(m, p)*z(n)*z(q) + zz2inv(n, p)*z(m)*z(q) + zz3inv(m, q)*z(n)*z(p) + zz4inv(n, q)*z(m)*z(p) + gamma(m, n, p, q);
                             end
                         end
                     end
                 end
                 gamma = gamma*(1/4);
                 Gamma_theta = gamma + Gamma_theta;
             end
             
             
             Gamma = Gamma_theta/8/pi^2/norm(x)*dtheta;
             
             GAM = [Gamma(1, 1, 1, 1) Gamma(1, 1, 2, 2) Gamma(1, 1, 3, 3) Gamma(1, 1, 1, 2) Gamma(1, 1, 1, 3) Gamma(1, 1, 2, 3);
                 Gamma(2, 2, 1, 1) Gamma(2, 2, 2, 2) Gamma(2, 2, 3, 3) Gamma(2, 2, 1, 2) Gamma(2, 2, 1, 3) Gamma(2, 2, 2, 3);
                 Gamma(3, 3, 1, 1) Gamma(3, 3, 2, 2) Gamma(3, 3, 3, 3) Gamma(3, 3, 1, 2) Gamma(3, 3, 1, 3) Gamma(3, 3, 2, 3);
                 Gamma(1, 2, 1, 1) Gamma(1, 2, 2, 2) Gamma(1, 2, 3, 3) Gamma(1, 2, 1, 2) Gamma(1, 2, 1, 3) Gamma(1, 2, 2, 3);
                 Gamma(1, 3, 1, 1) Gamma(1, 3, 2, 2) Gamma(1, 3, 3, 3) Gamma(1, 3, 1, 2) Gamma(1, 1, 1, 3) Gamma(1, 3, 2, 3);
                 Gamma(2, 3, 1, 1) Gamma(2, 3, 2, 2) Gamma(2, 3, 3, 3) Gamma(2, 3, 1, 2) Gamma(2, 3, 1, 3) Gamma(2, 3, 2, 3)];
             
             B(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) = GAM;
         elseif (I == J)
             B(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) = inv(DELTAC);
         end
         Comp0(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) = Comp0(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) + CZERO;
         Comp(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) = Comp(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) + CBIG;
         if (I == J)
             DC(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) = DELTAC;
         end
     end
 end
 
 ep1 = zeros(6*N, 1);
 ep2 = zeros(6*N, 1);
 ep3 = zeros(6*N, 1);
 ep4 = zeros(6*N, 1);
 ep5 = zeros(6*N, 1);
 ep6 = zeros(6*N, 1);
 f1 = zeros(6*N, 1);
 f2 = zeros(6*N, 1);
 f3 = zeros(6*N, 1);
 f4 = zeros(6*N, 1);
 f5 = zeros(6*N, 1);
 f6 = zeros(6*N, 1);
 f1(1, 1) = 1; f1(7, 1) = 1; f1(19, 1) = 1; f1(25, 1) = 1; f1(31, 1) = 1; f1(37, 1) = 1; f1(43, 1) = 1; f1(49, 1) = 1;
 for i = 1:6
     f1(6*(i-1)+1, 1) = 1;
     f2(6*(i-1)+2, 1) = 1;
     f3(6*(i-1)+3, 1) = 1;
     f4(6*(i-1)+4, 1) = 1;
     f5(6*(i-1)+5, 1) = 1;
     f6(6*(i-1)+6, 1) = 1;
 end
 for i = 1:N
     ep1(6*(i-1)+1:6*i, 1) = CBIG\f1(6*(i-1)+1:6*i, 1);
     ep2(6*(i-1)+1:6*i, 1) = CBIG\f2(6*(i-1)+1:6*i, 1);
     ep3(6*(i-1)+1:6*i, 1) = CBIG\f3(6*(i-1)+1:6*i, 1);
     ep4(6*(i-1)+1:6*i, 1) = CBIG\f4(6*(i-1)+1:6*i, 1);
     ep5(6*(i-1)+1:6*i, 1) = CBIG\f5(6*(i-1)+1:6*i, 1);
     ep6(6*(i-1)+1:6*i, 1) = CBIG\f6(6*(i-1)+1:6*i, 1);
 end

 temp = B*DC;
 eps1 = inv(temp)*ep1;
 eps2 = inv(temp)*ep2;
 eps3 = inv(temp)*ep3;
 eps4 = inv(temp)*ep4;
 eps5 = inv(temp)*ep5;
 eps6 = inv(temp)*ep6;
 sigma1 = Comp*eps1;
 sigma2 = Comp*eps2;
 sigma3 = Comp*eps3;
 sigma4 = Comp*eps4;
 sigma5 = Comp*eps5;
 sigma6 = Comp*eps6;
 
epep = zeros(6, 6);
sigsig = zeros(6, 6);
for i = 1:6
    for j = 1:N
        epep(i, 1) = epep(i, 1) + eps1(j*(i - 1) + 1)/N;
        sigsig(i, 1) = sigsig(i, 1) + sigma1(j*(i - 1) + 1)/N;
        epep(i, 2) = epep(i, 2) + eps2(j*(i - 1) + 1)/N;
        sigsig(i, 2) = sigsig(i, 2) + sigma2(j*(i - 1) + 1)/N;
        epep(i, 3) = epep(i, 3) + eps3(j*(i - 1) + 1)/N;
        sigsig(i, 3) = sigsig(i, 3) + sigma3(j*(i - 1) + 1)/N;
        epep(i, 4) = epep(i, 4) + eps4(j*(i - 1) + 1)/N;
        sigsig(i, 4) = sigsig(i, 4) + sigma4(j*(i - 1) + 1)/N;
        epep(i, 5) = epep(i, 5) + eps5(j*(i - 1) + 1)/N;
        sigsig(i, 5) = sigsig(i, 5) + sigma5(j*(i - 1) + 1)/N;
        epep(i, 6) = epep(i, 6) + eps6(j*(i - 1) + 1)/N;
        sigsig(i, 6) = sigsig(i, 6) + sigma6(j*(i - 1) + 1)/N;
    end
end

ee = [eps1 eps2 eps3 eps4 eps5 eps6];
ss = [sigma1 sigma2 sigma3 sigma4 sigma5 sigma6];

cfinal = ss*pinv(ee);
%Done