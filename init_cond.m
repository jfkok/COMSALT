function [v_x, v_z, spin] = init_cond (pEj, d, u_fr_surf, min_lift)

%init_cond sets the parameters of a newly lifted surface particle

%these are global parameters
global g rho d_mean splash_fraction u_fr_thr_soil init_mu_angle spin_sigma spin_mu u_fr rho_p;
aer_entrain = 0; %keeps track of whether the particle was aerodynamically entrained

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The following determines the lift-off speed    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (max(pEj)==0 && min_lift || (splash_fraction < rand(1) && u_fr_thr_soil < u_fr_surf))  %in this case, aerodynamic lifting takes place in steady-state saltation because the fluid threshold is below the friction velocity at the surface OR this is the initial iteration, and the model is initiated by aerodynamic lifting of particles with lowest fluid threshold
    %This describes aerodynamic entrainment.
    aer_entrain = 1;
    shields_no = rho*u_fr^2/((rho_p-rho)*g*d);
    w0x_avg = 9.0*u_fr; %after Eq. 26 in Hu and Hui (1996)
    P = rand(1); %generating a random number
    w0x = log(1/(1-P))*w0x_avg; %assuming that the lift-off velocity is an exponential distribution with mean of w0_avg
    if (shields_no < 1.2) %after Eq. 28 in Hu and Hui (1996)
        w0z_avg = (3.2-4.5*log(shields_no))*u_fr;
    else
        w0z_avg = 3.1*u_fr;
    end
    P = rand(1); %generating a random number
    w0z = log(1/(1-P))*w0z_avg; %assuming that the lift-off velocity is an exponential distribution with mean of w0_avg    
    w0 = sqrt(w0x^2+w0z^2);
elseif (max(pEj)==0) %in this case, the model has just been initialized, and pEj has not yet been obtained. So particle lifting is done aerodynamically, assuming that the lift-off velocity is 2*sqrt(2*g*d_mean). This assumption does not affect the model results
    w0_avg = 2*sqrt(2*g*d_mean)/sin(3.14159*init_mu_angle/180); %It is assumed that, on average, an aerodynamically entrained particle has enough vertical velocity to reach an average height of ~4 diameteters when non-gravitational forces are neglected
    P = rand(1); %generating a random number
    w0 = log(1/(1-P))*w0_avg; %assuming that the lift-off velocity is an exponential distribution with mean of w0_avg
else %if there is no aerodynamic entrainment, which is almost always the case in steady-state saltation, w0 is randomly determined from the distribution of ejected velocities
    a = rand(1); %generating a random number
    w0 = real(pEj(round((a*(size(pEj,2)-1))+1))); %the speed of ejected surface particles is derived from the running list pEj of characteristics of splashed surface particles
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The following determines the lift-off angle    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (aer_entrain) %in this case, aerodynamic lifting takes place in steady-state saltation because the fluid threshold is below the friction velocity at the surface
    theta = atan(w0z/w0x);
else    
    P = rand(1); %generating a random number
    theta = log(1/(1-P))*init_mu_angle; %determining the lift-off angle from an exponential distribution with mean init_mu_angle
    if (theta>179)  %in case theta is unrealistically large
        theta = 179;
    end %if, in case theta is unrealistically small
    if (theta<1)  %in case theta is unrealistically large (this chance is 2%)
        theta = 1;
    end %if, in case theta is unrealistically small
    theta = (theta/180)*3.14159;  %converting theta to radiants
end

%Finally, the lift-off speed of the particle can be set:
v_x = cos(theta)*w0; %the horizontal lift-off speed
v_z = sin(theta)*w0; %the vertical lift-off speed
spin = Gaussian(spin_mu, spin_sigma); %the initial particle spin is assumed to follow a Gaussian distribution