clear all
clc

% Defining the geometry of the container

x_container = 5;
y_container = 5;
z_container = 5;

% Defining radii of the discrete elements and an initial guess of their
% count

r_avg = 0.5;
r_stddev = 0.1;
den = 0.65;

% Packing the circles inside the domain

[x, y, z, r, n, neighbors, Neigh] = PackGen(x_container, y_container, z_container, r_avg , r_stddev, den)

% Finding the number of elements and their nodes
q = 1;
for i = 1:n
    for j = i:n
        if abs(sqrt((x(j, 1) - x(i, 1))^2 + (y(j, 1) - y(i, 1))^2 + (z(j, 1) - z(i, 1))^2) - r(j) - r(i)) < 0.1
            element(q, 1) = i;
            element(q, 2) = j;
            l(q, 1) = r(i) + r(j);
            q = q + 1;
        end
    end
end
element_count = q - 1;

% Relationships for plane stress problem

E = 1.461;  %modulus of elasticity
nu = 0.13;  %poisson's ratio

% Relationship between stiffness values and elasticity properties

% For plane stress problem
% k_n = E/(sqrt(3)*(1 - nu));
% k_s = E*(1 - 3*nu)/(sqrt(3)*(1 - nu^2));

% For plane strain problem
k_n = E/(sqrt(3)*(1-2*nu)*(1+nu));
k_s = E*(1-4*nu)/(sqrt(3)*(1-2*nu)*(1+nu));
k_t = k_s;

%% Assembly of the stiffness matrix
K_local = zeros(12, 12);
K_global = zeros(6*n, 6*n);


for i = 1:element_count
    kp11 = [k_n 0 0 0 0 0; 0 k_s 0 0 0 k_s*l(i); 0 0 k_t 0 -k_t*l(i) 0; 0 0 0 0 0 0; 0 0 -k_t*l(i) 0 k_t*l(i)^2 0; 0 k_s*l(i) 0 0 0 k_s*l(i)^2];
    kp12 = [-k_n 0 0 0 0 0; 0 -k_s 0 0 0 k_s*l(i); 0 0 -k_t 0 -k_t*l(i) 0; 0 0 0 0 0 0; 0 0 k_t*l(i) 0 k_t*l(i)^2 0; 0 -k_s*l(i) 0 0 0 k_s*l(i)^2];
    kp21 = [-k_n 0 0 0 0 0; 0 -k_s 0 0 0 -k_s*l(i); 0 0 -k_t 0 k_t*l(i) 0; 0 0 0 0 0 0; 0 0 -k_t*l(i) 0 k_t*l(i)^2 0; 0 k_s*l(i) 0 0 0 k_s*l(i)^2];
    kp22 = [k_n 0 0 0 0 0; 0 k_s 0 0 0 -k_s*l(i); 0 0 k_t 0 k_t*l(i) 0; 0 0 0 0 0 0; 0 0 k_t*l(i) 0 k_t*l(i)^2 0; 0 -k_s*l(i) 0 0 0 k_s*l(i)^2];

    kp = [kp11 kp12; kp21 kp22];
    
    [kp2] = rot(kp, i, element, l, x, y, z, r, kp11, kp12, kp21, kp22);
    
    K_local = kp2;
    
    a = element(i, 1);
    b = element(i, 2);
    K_global((6*(a - 1) + 1):6*a, (6*(b - 1) + 1):6*b) = K_global((6*(a - 1) + 1):6*a, (6*(b - 1) + 1):6*b) + K_local(1:6, 7:12);
    K_global((6*(a - 1) + 1):6*a, (6*(a - 1) + 1):6*a) = K_global((6*(a - 1) + 1):6*a, (6*(a - 1) + 1):6*a) + K_local(1:6, 1:6);
    K_global((6*(b - 1) + 1):6*b, (6*(a - 1) + 1):6*a) = K_global((6*(b - 1) + 1):6*b, (6*(a - 1) + 1):6*a) + K_local(7:12, 1:6);
    K_global((6*(b - 1) + 1):6*b, (6*(b - 1) + 1):6*b) = K_global((6*(b - 1) + 1):6*b, (6*(b - 1) + 1):6*b) + K_local(7:12, 7:12);
end

%% Solving a sample problem (elastic compression)

% F = zeros(3*circle_count, 1);
% u = zeros(3*circle_count, 1);
% 
% % Choosing the top layer nodes for applying force
% force_nodes = [];
% for i = 1:circle_count
%     if y(i) > (y_container - r)
%         force_nodes = [force_nodes i];
%     end
% end
% 
% for i = 1:length(force_nodes)
%     j = 3*(force_nodes(i) - 1) + 2;
%     F(j, 1) = -10/length(force_nodes);
% end
% 
% % Choosing the bottom layer as the fixed nodes
% fixed_nodes = [];
% fixed_dofs = [];
% for i = 1:circle_count
%     if y(i) == 0
%         fixed_nodes = [fixed_nodes i];
%         fixed_dofs = [fixed_dofs (3*i-2) (3*i-1) (3*i)];
%     end
% end
% free_dofs = [1:(3*circle_count)];
% free_dofs(fixed_dofs) = [];
% 
% u(free_dofs) = K_global(free_dofs, free_dofs)\F(free_dofs);
% 
%  for i = 1:3*(circle_count - 26)
% l1 = linspace(x(element(i, 1)), x(element(i, 2)), 100);
% l2 = linspace(y(element(i, 1)), y(element(i, 2)), 100);
% plot(l1, l2);
% hold on
% end
% hold on
% for i = 1:3*(circle_count-26)
% l1 = linspace(x(element(i, 1)) + u(3*element(i, 1) - 2), x(element(i, 2)) + u(3*element(i, 2) - 2), 100);
% l2 = linspace(y(element(i, 1)) + u(3*element(i, 1) - 1), y(element(i, 2)) + u(3*element(i, 2) - 1), 100);
% plot(l1, l2, 'r-');
% hold on
% end