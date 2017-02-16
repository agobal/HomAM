function [xb, yb, rb, n, neighbors] = TwoDPackGen(xb_container, yb_container, aveg, stddev, rden)

phi = rden;         %The desired packing density
V_particles = 0;    %Total volume of the particles inside the box

%Creating the desired number of particles with the desired sizes
for i = 1:2000
    rb(i, 1) = (stddev*rand + aveg);
    V_particles = V_particles + pi*(rb(i,1)^2);
    if V_particles >= (xb_container*yb_container*phi)
        n = i - 1;
        break;
    end
end

%initialize the location of particles
xb = zeros(n, 1);
yb = zeros(n, 1);
overlap = zeros(n, 1);
neighbors = zeros(n,10);

%Random positioning of the particles in a box
for i = 1:n
    xb(i, 1) = rand*(xb_container - 2*rb(i, 1)) + rb(i, 1);
    yb(i, 1) = rand*(yb_container - 2*rb(i, 1)) + rb(i, 1);
end

%Relocation of the particles with respect to the desired maximum overlap
for k = 1:600
for i = 1:n
    q = 0;
    xn = 0;
    yn = 0;
    for j = 1:n
        if i ~= j
            overlap(j, 1) = (rb(i, 1) + rb(j, 1) - sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2))/(rb(i, 1) + rb(j, 1));
            if overlap(j, 1) > 0
                q = q + 1;
                xn = xn + xb(j, 1) + (xb(i, 1) - xb(j, 1))*(rb(i, 1) + rb(j, 1))/sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2);
                yn = yn + yb(j, 1) + (yb(i, 1) - yb(j, 1))*(rb(i, 1) + rb(j, 1))/sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2);
            end
        end
    end
    if q >= 1
        xb(i, 1) = xn/q;
        yb(i, 1) = yn/q;
        
        %Dont accept if lie outside the boundary
        if xb(i, 1) >= (xb_container - rb(i,1))
            xb(i, 1) = (xb_container - rb(i,1));
        end
        if yb(i, 1) >= (yb_container - rb(i,1))
            yb(i, 1) = (yb_container - rb(i,1));
        end
        if xb(i, 1) <= rb(i,1)
            xb(i, 1) = rb(i,1);
        end
        if yb(i, 1) <= rb(i,1)
            yb(i, 1) = rb(i,1);
        end
    end
end
end

%Finding the neighbors of each particle
for i = 1:n
    q = 1;
    for j = 1:n
        if (i ~= j) && ((sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2) - (rb(i, 1) + rb(j, 1))) <= 0.000001)
            neighbors(i, q) = j;
            q = q + 1;
        end
    end
end

%Finding the maximum overlap existing in the packing
max_overlap = 0;
for i = 2:n
    for j = 1:(i-1)
        if (sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2) < (rb(i, 1) + rb(j, 1)))
            max_overlap = ((rb(i, 1) + rb(j, 1)) - sqrt((xb(i, 1) - xb(j, 1))^2 + (yb(i, 1) - yb(j, 1))^2))/(rb(i, 1) + rb(j, 1));
        end
    end
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

% dlmwrite('particle_x.txt', xb);
% dlmwrite('particle_y.txt', yb);
% dlmwrite('particle_r.txt', rb);
% dlmwrite('Neighbors.txt', neighbors);