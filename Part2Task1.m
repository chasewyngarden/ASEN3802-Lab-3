
clear;
clc;
close all;


% Define constants
N = 5;
b = 100;
a0_t = 2*pi;
a0_r = 2*pi;
c_t = 8;
c_r = 10;
aero_t = 0;
aero_r = 0;
geo_t = 5*pi/180;
geo_r = 5*pi/180;


%Call funtion
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);


function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

% solving for A1 and An 

    m = 1:N;

    theta_i = (m*pi)/(2*N); %define theta 


    % Variables depending on span
    ctheta_i = c_r -(cos(theta_i)*(c_r - c_t));

    a0_theta_i  = a0_r - (cos(theta_i)*(a0_r - a0_t));

    alphaL0_theta_i = aero_r - (cos(theta_i)*(aero_r - aero_t));

    AOA_theta_i = geo_r - (cos(theta_i)*(geo_r - geo_t));

    % Effective alpha

    di = (AOA_theta_i - alphaL0_theta_i)';

    %
    
    for i = 1:N

      for j = 1:N

        Mij(i,j) = ((4*b)/(a0_theta_i(i)*ctheta_i(i)))*sin((2*j - 1)*theta_i(i)) + (2*j - 1)*((sin((2*j - 1)*theta_i(i)))/sin(theta_i(i)));

      end  
       
    end 

% solve matrix using \
A = Mij\di;


% find AR and caluclate CL
cave = (c_t + c_r)/2;
AR = b/cave;

c_L = A(1)*pi*AR;

% delta calculation
    del = 0;
for k = 2:N

    del = del + ((2*k -1) * (A(k)/A(1))^2);

end 

    %CD and E calulation based on delta
    c_Di = ((c_L)^2) / (pi*AR) * (1 + del);

    e = 1/(1+del);
end
