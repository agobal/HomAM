function [xb, yb, zb, rb, n, neighbors, Neigh] = PackGen(xb_container, yb_container, zb_container, aveg,stddev,rden)
tic
%n = 150;%44643
phi = rden;
V_particles = 0;
for i = 1:20000
%     rb(i, 1) = (0.000005*rand + 0.0000075);
    rb(i, 1) = (stddev*rand + aveg);
    V_particles = V_particles + (4/3)*pi*(rb(i,1)^3);
    if V_particles >= (xb_container*yb_container*zb_container*phi)
        n = i - 1;
        break;
    end
end

xb = zeros(n, 1);
yb = zeros(n, 1);
zb = zeros(n, 1);
overlap = zeros(n, 1);
% phi = 0.64;
% V_container = V_particles/phi;
% zb_container = (V_container)^(1/3);
% xb_container = zb_container;
% yb_container = zb_container;
neighbors = zeros(n,10);
%Random positioning of the particles in cylindrical domain
for i = 1:n
    xb(i, 1) = rand*(xb_container - 2*rb(i, 1)) + rb(i, 1);
    yb(i, 1) = rand*(yb_container - 2*rb(i, 1)) + rb(i, 1);
    zb(i, 1) = rand*(zb_container - 2*rb(i, 1)) + rb(i, 1);
end
%Relocation
for k = 1:600
for i = 1:n
    q = 0;
    xn = 0;
    yn = 0;
    zn = 0;
    for j = 1:n
        if i ~= j
            overlap(j, 1) = (rb(i, 1) + rb(j, 1) - sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2 + (zb(i, 1) - zb(j, 1))^2))/(rb(i, 1) + rb(j, 1));
            if overlap(j, 1) > 0
                q = q + 1;
                xn = xn + xb(j, 1) + (xb(i, 1) - xb(j, 1))*(rb(i, 1) + rb(j, 1))/sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2 + (zb(i, 1) - zb(j, 1))^2);
                yn = yn + yb(j, 1) + (yb(i, 1) - yb(j, 1))*(rb(i, 1) + rb(j, 1))/sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2 + (zb(i, 1) - zb(j, 1))^2);
                zn = zn + zb(j, 1) + (zb(i, 1) - zb(j, 1))*(rb(i, 1) + rb(j, 1))/sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2 + (zb(i, 1) - zb(j, 1))^2);
            end
        end
    end
    if q >= 1
        xb(i, 1) = xn/q;
        yb(i, 1) = yn/q;
        zb(i, 1) = zn/q;
        %Dont accept if lie outside the boundary
        if xb(i, 1) >= (xb_container - rb(i,1))
            xb(i, 1) = (xb_container - rb(i,1));
        end
        if yb(i, 1) >= (yb_container - rb(i,1))
            yb(i, 1) = (yb_container - rb(i,1));
        end
        if zb(i, 1) >= (zb_container - rb(i,1))
            zb(i, 1) = (zb_container - rb(i,1));
        end
        if xb(i, 1) <= rb(i,1)
            xb(i, 1) = rb(i,1);
        end
        if yb(i, 1) <= rb(i,1)
            yb(i, 1) = rb(i,1);
        end
        if zb(i, 1) <= rb(i,1)
            zb(i, 1) = rb(i,1);
        end
    end
end
end

for i = 1:n
    q = 1;
    for j = 1:n
        if (i ~= j) && ((sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2 + (zb(i, 1) - zb(j, 1))^2) - (rb(i, 1) + rb(j, 1))) <= 0.000001)
            neighbors(i, q) = j;
            q = q + 1;
        end
    end
end

% q = 0;
% for i = 1:n
%     if ((xb(i, 1) >= rb(i,1)) && (xb(i, 1) <= (xb_container - rb(i,1))) && (zb(i, 1) >= rb(i,1)) && (zb(i, 1) <= (zb_container - rb(i,1))) && (yb(i, 1) >= rb(i,1)) && (yb(i, 1) <= (yb_container - rb(i,1))))
%         q = q + 1;
%         xf(q, 1) = xb(i, 1);
%         yf(q, 1) = yb(i, 1);
%         zf(q, 1) = zb(i, 1);
%         rf(q, 1) = rb(i, 1);
%     end
% end

max_overlap = 0;
for i = 2:n
    for j = 1:(i-1)
        if (sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2 + (zb(i, 1) - zb(j, 1))^2) < (rb(i, 1) + rb(j, 1)))
            max_overlap = ((rb(i, 1) + rb(j, 1)) - sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2 + (zb(i, 1) - zb(j, 1))^2))/(rb(i, 1) + rb(j, 1));
        end
    end
end

[X,Y,Z] = sphere(25);
for i=1:n
    surf(X*rb(i)+xb(i),Y*rb(i)+yb(i),Z*rb(i)+zb(i));
    hold on
end

sizen = size(neighbors);
Neigh = zeros(sizen(2), 1);
for i = 1:n
for j = 1:sizen(2)
if neighbors(i, j) == 0
Neigh(j, 1) = Neigh(j, 1) + 1;
break;
end
end
end

dlmwrite('particle_x.txt', xb);
dlmwrite('particle_y.txt', yb);
dlmwrite('particle_z.txt', zb);
dlmwrite('particle_r.txt', rb);
dlmwrite('Neighbors.txt', neighbors);

toc