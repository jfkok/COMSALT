function [v_t] = Terminal_velocity (d)

%This function calculates the particle's terminal velocity

%these are global parameters used by Terminal_velocity
global g mu rho rho_p

sigma = rho_p/rho;
v_t = 1; %guessing the initial terminal velocity
delta = 1; %the relative difference in v_t between subsequent iterations in the while loop below
while (delta>0.0001) %iterating until v_t converges
    Re = (d*rho/mu)*v_t;  %the Reynolds number
    Cd = ((32/Re)^(2/3)+1)^(3/2); %the drag coefficient
    v_t_old = v_t; %saving the old v_t for comparison with the new one
    v_t = 4*g*d*sigma/(3*Cd*v_t); %recalculating v_t
    delta = abs((v_t - v_t_old)/v_t); %the relative difference in v_t between subsequent iterations in the while loop below
end %while, iterating until v_t converges