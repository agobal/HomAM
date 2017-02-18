xb_container = 200;
yb_container = 200;

rden = 0.83;

aveg = 15;
stddev = 0;

[xbb, ybb, rbb, nb, neighborsb] = TwoDPackGen(xb_container, yb_container, aveg, stddev, rden);

% for i = 1:nb
%     x = xbb(i);
%     y = ybb(i);
%     r = 25;
%     ang = 0:0.01:2*pi;
%     xp = r*cos(ang);
%     yp = r*sin(ang);
%     plot(x+xp,y+yp);
%     hold on
% end

time_points = 10000;
dt = 0.01;

Tb = zeros(nb, time_points);
k = 1;

%Boundary conditions
for i = 1:nb
    if xbb(i) < 20
        Tb(i, :) = 1;
    end
    if xbb(i) > 180
        Tb(i, :) = 0;
    end
end

%Solving the heat transfer problem
for t = 1:(time_points - 1)
    for i = 1:nb
        Q = 0;
        for j = 1:10
            if neighborsb(i, j) ~= 0
                Q = Q + (Tb(neighborsb(i, j), t) - Tb(i, t))*k;
            end
        end
        Tb(i, t + 1) = Tb(i, t) + Q*dt;
        if xbb(i) < 20
            Tb(i, :) = 1;
        end
        if xbb(i) > 180
            Tb(i, :) = 0;
        end
    end
end
 plot(xbb(:, 1), Tb(:, 9999), 'ro');
%fit data to temperatures with respect to distance from the point of
%temperature and find out a formula then apply it to find the temperature
%of points from the neighboring points. the multipliers here must follow
%the partition of unity quality so we pressure them to do that.
               
% 
% sigma = 8;
% mu = 0;
% f = 20/(sigma*sqrt(2*pi))*exp(-0.5*((x-mu)/sigma).^2);

multiplier = zeros(n, 1);
sfn = zeros(n, 1);
shapef = zeros(n, nb);
for i = 1:n
    for j = 1:nb
        d(i, j) = sqrt((xb(i, 1) - xbb(j, 1))^2 + (yb(i, 1) - ybb(j, 1))^2);
        if d(i, j) < 21
            shapef(i, j) = f(1)*d(i, j)^4 + f(2)*d(i, j)^3 + f(3)*d(i, j)^2 + f(4)*d(i, j) + f(5);
            sfn(i, 1) = sfn(i, 1) + shapef(i, j);
        else
            shapef(i, j) = 0;
            sfn(i, 1) = sfn(i, 1) + shapef(i, j);
        end
    end
    multiplier(i, 1) = 1/sfn(i, 1);
end
hold on

Ttat = zeros(n, 1);
for i = 1:n
    for j = 1:nb
        if d(i, j) < 25
            Ttat(i, 1) = Ttat(i, 1) + multiplier(i, 1)*shapef(i, j)*Tb(j, 9999);
        end
    end
end

%     for i = 1:num
%     x = xb(i);
%     y = yb(i);
%     r = 1.5;
%     ang = 0:0.01:2*pi;
%     xp = r*cos(ang);
%     yp = r*sin(ang);
%     plot(x+xp,y+yp);
% hold on;
% text(xb(i), yb(i), num2str(i));
% hold on
% end
% list = [7, 56, 38, 49, 24, 101, 51, 27, 85];
% for i = 1:9
% x = xb(list(i));
% y = yb(list(i));
% r = 25;
% ang = 0:0.01:2*pi;
% xp = r*cos(ang);
% yp = r*sin(ang);
% plot(x+xp,y+yp,'r');
% hold on
% end