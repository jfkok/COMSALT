function [x_prev, v_xR, v_zR, alpha_R] = bounce_trajectory (x_prev, z, z_prev, v_xI, v_zI, d)

%bounce_trajectory calculates the rebound parameters for a given impacting particle

%these are global parameters
global mu_angle mu_alpha_R sigma_alpha_R pi;

vI = sqrt(v_xI^2+v_zI^2);  %the impact speed
angle_imp = (180/pi)*atan(-v_zI/v_xI); %rebound angle in degrees
P = rand(1); %a random number
angle_reb = log(1/(1-P))*mu_angle; %determining angle_reb from an exponential distribution
if (angle_reb<1)  %in case theta is unrealistically small (this chance is ~0.1%)
    angle_reb = 1;
end %if, in case theta is unrealistically small
if (angle_reb>179)  %in case theta is unrealistically large (this chance is ~0.1%)
    angle_reb = 179;
end %if, in case theta is unrealistically small
angle_reb = (angle_reb/180)*3.14159;  %converting theta to radiants
alpha_R = min(0.999,sqrt(max(0.01,Gaussian(mu_alpha_R^2,sigma_alpha_R^2)))); %limiting alpha_R^2 between 0.01 and 0.999

%The rebound parameters are then:
v_xR = alpha_R*vI*cos(angle_reb); %horizontal speed of rebounding grain
v_zR = alpha_R*vI*sin(angle_reb); %vertical speed of rebounding grain
del_t = -((z_prev-d/2)/v_zI); %the time difference between t(k,i-2)+h and the time when the particle actually struck the surface. This is used to correct t(k,i-1) in saltation_trajectory
x_prev = del_t*v_xI+x_prev; %this is the horizontal distance at which the particle struck the surface