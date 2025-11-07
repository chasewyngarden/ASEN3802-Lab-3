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


c = 1;
x = linspace(0,c, 1000); % x is between 0 and c

% Givens for NACA 4415 Airfoil
t = 0.15;
m = 0.04;
p = 0.4;


[x_b_0006, y_b_0006] = NACA_Airfoil_gen (c, 0.06, 0, 0, x); % NACA 0006 Airfoil

[x_b_0012, y_b_0012] = NACA_Airfoil_gen (c, 0.12, 0, 0, x); % NACA 0006 Airfoil

[x_b_0018, y_b_0018] = NACA_Airfoil_gen (c, 0.18, 0, 0, x); % NACA 0018 Airfoil

[x_b_2418, y_b_2418] = NACA_Airfoil_gen (c, 0.18, 0.02, 0.4, x); % NACA 2418 Airfoil

CL = Vortex_Panel(x_b_0012, y_b_0012, 5);



% FUNCTIONS

function [x_b, y_b] = NACA_Airfoil_gen (c, t, m, p, x) % function to plot NACA 4 digit airfoils

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

    x_b = [ x_U(end:-1:1), x_L(2:end) ]; % boundary matrix which reverses direction so starts at TE and goes cw (both Upper and lower)
    
    y_U = y_c - y_t.*cos(xi);
    y_L = y_c + y_t.*cos(xi);

    y_b = [ y_U(end:-1:1), y_L(2:end) ];

    
    figure;
    hold on
    plot (x_b, y_b, 'k', 'LineWidth', 1.5);
    if m ~= 0 && p ~= 0
        plot (x, y_c, 'r', 'LineWidth', 1.5) % plots camber line
    end
    % plot (x_b(1), y_b(1), 'ro'); % Test plots to ensure proper direction, 
    % plot (x_b(100), y_b(100), 'ro'); % Clockwise starting from TE
    axis equal
    xlabel('x / c');
    ylabel('y / c');
    title(['Plot of NACA ', num2str(m*100), num2str(p*10), num2str(t*100), ' Airfoil']);
    hold off;


end




function [CL] = Vortex_Panel(XB,YB,ALPHA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                           %
%                                  %
% XB  = Boundary Points x-location %
% YB  = Boundary Points y-location %
% ALPHA = AOA in degrees           %
%                                  %
% Output:                          %
%                                  %
% CL = Sectional Lift Coefficient  %
% improves efficiency by preallocating matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Convert to Radians %
%%%%%%%%%%%%%%%%%%%%%%

ALPHA = ALPHA*pi/180;

%%%%%%%%%%%%%%%%%%%%%
% Compute the Chord %
%%%%%%%%%%%%%%%%%%%%%

CHORD = max(XB)-min(XB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Number of Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = max(size(XB,1),size(XB,2))-1;
MP1 = M+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate Matrices for Efficiency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(1,M);
Y = zeros(1,M);
S = zeros(1,M);
THETA = zeros(1,M);
SINE = zeros(1,M);
COSINE = zeros(1,M);
RHS = zeros(1,M);
CN1 = zeros(M);
CN2 = zeros(M);
CT1 = zeros(M);
CT2 = zeros(M);
AN = zeros(M);
AT = zeros(M);
V = zeros(1,M);
CP = zeros(1,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-Panel Relationships:                                  %
%                                                             %
% Determine the Control Points, Panel Sizes, and Panel Angles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    IP1 = I+1;
    X(I) = 0.5*(XB(I)+XB(IP1));
    Y(I) = 0.5*(YB(I)+YB(IP1));
    S(I) = sqrt( (XB(IP1)-XB(I))^2 +( YB(IP1)-YB(I))^2 );
    THETA(I) = atan2( YB(IP1)-YB(I), XB(IP1)-XB(I) );
    SINE(I) = sin( THETA(I) );
    COSINE(I) = cos( THETA(I) );
    RHS(I) = sin( THETA(I)-ALPHA );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:             %
%                                        %
% Determine the Integrals between Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    for J = 1:M
        if I == J
            CN1(I,J) = -1.0;
            CN2(I,J) = 1.0;
            CT1(I,J) = 0.5*pi;
            CT2(I,J) = 0.5*pi;
        else
            A = -(X(I)-XB(J))*COSINE(J) - (Y(I)-YB(J))*SINE(J);
            B = (X(I)-XB(J))^2 + (Y(I)-YB(J))^2;
            C = sin( THETA(I)-THETA(J) );
            D = cos( THETA(I)-THETA(J) );
            E = (X(I)-XB(J))*SINE(J) - (Y(I)-YB(J))*COSINE(J);
            F = log( 1.0 + S(J)*(S(J)+2*A)/B );
            G = atan2( E*S(J), B+A*S(J) );
            P = (X(I)-XB(J)) * sin( THETA(I) - 2*THETA(J) ) ...
              + (Y(I)-YB(J)) * cos( THETA(I) - 2*THETA(J) );
            Q = (X(I)-XB(J)) * cos( THETA(I) - 2*THETA(J) ) ...
              - (Y(I)-YB(J)) * sin( THETA(I) - 2*THETA(J) );
            CN2(I,J) = D + 0.5*Q*F/S(J) - (A*C+D*E)*G/S(J);
            CN1(I,J) = 0.5*D*F + C*G - CN2(I,J);
            CT2(I,J) = C + 0.5*P*F/S(J) + (A*D-C*E)*G/S(J);
            CT1(I,J) = 0.5*C*F - D*G - CT2(I,J);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:           %
%                                      %
% Determine the Influence Coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    AN(I,1) = CN1(I,1);
    AN(I,MP1) = CN2(I,M);
    AT(I,1) = CT1(I,1);
    AT(I,MP1) = CT2(I,M);
    for J = 2:M
        AN(I,J) = CN1(I,J) + CN2(I,J-1);
        AT(I,J) = CT1(I,J) + CT2(I,J-1);
    end
end
AN(MP1,1) = 1.0;
AN(MP1,MP1) = 1.0;
for J = 2:M
    AN(MP1,J) = 0.0;
end
RHS(MP1) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the gammas %
%%%%%%%%%%%%%%%%%%%%%%%%

GAMA = AN\RHS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Tangential Veloity and Coefficient of Pressure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    V(I) = cos( THETA(I)-ALPHA );
    for J = 1:MP1
        V(I) = V(I) + AT(I,J)*GAMA(J);
    end
    CP(I) = 1.0 - V(I)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Sectional Coefficient of Lift %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIRCULATION = sum(S.*V);
CL = 2*CIRCULATION/CHORD;

end









