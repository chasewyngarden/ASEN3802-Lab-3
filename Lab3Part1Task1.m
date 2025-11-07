% ASEN 3802 Lab 3 Aerodynamics Part 1 Task 1 NACA 4-Digit Generator

clear;
clc;
close all;

% c = chord length
% x = position along the chord from 0 to c
% y_t = the half thickness at a given x (mean camber line to surface)
% t = max thickness in percent chord
% m = maximum camber in percent chord
% p = location of max camber


function [x_b, y_b] = NACA_Airfoil_gen (c, t, m, p, x)

    y_t = ((t*c)/0.2) * (0.2969*sqrt(x/c) - 0.126*(x/c) - 0.3516*(x/c).^2 + 0.2843*(x/c).^3 - 0.1036*(x/c).^4); % thickness distribution of airfoil normal to the mean camber line
    
    y_c = 0; % formula for the mean camber line
    
    for i = 1:length(x)   % looping through x vector
        if i >= 0 && i < p*c
        
            y_c = ((m*x)/(p^2)) .* (2*p - (x/c));
        
        end
        
        if i >= p*c && i <= c
        
            y_c = ((m*(c-x)) / ((1-p)^2)) .* (1 + (x/c) - 2*p);
        
        end
    end
    
    dy_c_dx = diff(y_c) / diff(x);
    
    xi = atan(dy_c_dx);  % local angle
    
    % coordinates of upper(U) and lower(L) airfoil surface
    x_U = x - y_t.*sin(xi);  % i changes sign directions to go clockwise
    x_L = x + y_t.*sin(xi);

    x_b = [ x_U(end:-1:1), x_L(2:end) ]; % boundary matrix which reverses direction so starts at TE and goes cw
    
    y_U = y_c - y_t.*cos(xi);
    y_L = y_c + y_t.*cos(xi);

    y_b = [ y_U(end:-1:1), y_L(2:end) ];

end


c = 1;
t = 0.15;
m = 0.04;
p = 0.4;

x = linspace(0,c, 1000); % x is between 0 and c


[x_b, y_b] = NACA_Airfoil_gen (c, t, m, p, x);


figure;
hold on
plot (x_b, y_b, 'k');
% plot (x_L, y_L, 'k');
plot (x_b(1), y_b(1), 'ro');
plot (x_b(100), y_b(100), 'ro');
axis equal


