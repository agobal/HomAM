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
Ex = 0.7406
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
             delta_c = c - c0;
             DELTAC = [delta_c(1, 1, 1, 1) delta_c(1, 1, 2, 2) delta_c(1, 1, 3, 3) delta_c(1, 1, 1, 2) delta_c(1, 1, 1, 3) delta_c(1, 1, 2, 3);
                 delta_c(2, 2, 1, 1) delta_c(2, 2, 2, 2) delta_c(2, 2, 3, 3) delta_c(2, 2, 1, 2) delta_c(2, 2, 1, 3) delta_c(2, 2, 2, 3);
                 delta_c(3, 3, 1, 1) delta_c(3, 3, 2, 2) delta_c(3, 3, 3, 3) delta_c(3, 3, 1, 2) delta_c(3, 3, 1, 3) delta_c(3, 3, 2, 3);
                 delta_c(1, 2, 1, 1) delta_c(1, 2, 2, 2) delta_c(1, 2, 3, 3) delta_c(1, 2, 1, 2) delta_c(1, 2, 1, 3) delta_c(1, 2, 2, 3);
                 delta_c(1, 3, 1, 1) delta_c(1, 3, 2, 2) delta_c(1, 3, 3, 3) delta_c(1, 3, 1, 2) delta_c(1, 1, 1, 3) delta_c(1, 3, 2, 3);
                 delta_c(2, 3, 1, 1) delta_c(2, 3, 2, 2) delta_c(2, 3, 3, 3) delta_c(2, 3, 1, 2) delta_c(2, 3, 1, 3) delta_c(2, 3, 2, 3)];
             B(6*(I - 1) + 1:6*I, 6*(J - 1) + 1:6*J) = DELTAC;
         end
     end
 end