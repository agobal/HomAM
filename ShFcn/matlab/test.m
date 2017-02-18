xb_container = 200;
yb_container = 200;

rden = 0.83;

aveg = 5;
stddev = 0;

[xb, yb, rb, n, neighbors] = TwoDPackGen(xb_container, yb_container, aveg, stddev, rden);

% for i = 1:n
%     x = xb(i);
%     y = yb(i);
%     r = rb(i);
%     ang = 0:0.01:2*pi;
%     xp = r*cos(ang);
%     yp = r*sin(ang);
%     plot(x+xp,y+yp,'r');
%     hold on
% end

time_points = 10000;
dt = 0.01;

T = zeros(n, time_points);
k = 1;

%Boundary conditions
for i = 1:n
    if xb(i) < 20
        T(i, :) = 1;
    end
    if xb(i) > 180
        T(i, :) = 0;
    end
end

%Solving the heat transfer problem
for t = 1:(time_points - 1)
    for i = 1:n
        Q = 0;
        for j = 1:10
            if neighbors(i, j) ~= 0
                Q = Q + (T(neighbors(i, j), t) - T(i, t))*k;
            end
        end
        T(i, t + 1) = T(i, t) + Q*dt;
        if xb(i) < 20
            T(i, :) = 1;
        end
        if xb(i) > 180
            T(i, :) = 0;
        end
    end
end
plot(xb(:, 1), T(:, 9999), 'o');
%fit data to temperatures with respect to distance from the point of
%temperature and find out a formula then apply it to find the temperature
%of points from the neighboring points. the multipliers here must follow
%the partition of unity quality so we pressure them to do that.
                
        
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