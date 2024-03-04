function [n, u, u_acc, tau_p, tau_p_acc, tau_p_surface, u_fr_surf_new] = wind_profile (d, n_relative, particle_shear_stress, particle_shear_stress_acc, total_surface_shear_stress, initial_iteration, n_power, n_old, u_fr_surf)

%wind_profile calculates the wind profile based on the shear velocity and
%the particle shear stress as a function of height

%these are global parameters used by wind_profile
global j_end j_acc_end u_fr u_fr_thr u_fr_thr_soil delta_z rho z0 delta_z_acc Owen pi rho_p splash_fraction no_fine_grid Mars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Initializing various parameters    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = (pi/6)*rho_p*d.^3; %the particle mass in kg
z = delta_z*(1:1:j_end); %the vertical coordinate on the coarse grid (don't confuse with the particle position z as used in saltation_trajectory and saltation_model!!!)
z_acc = delta_z_acc*(1:1:j_acc_end); %the vertical coordinate on the fine grid 
n_negative = false; %n_negative is set to true when n < 0 (this can happen when the winds are much too low - having n < 0 will reduce the running average of n used in saltation_model and increase the wind profile towards steady state)
if (total_surface_shear_stress==0) %the surface shear stress can not be zero
    error('total_surface_shear_stress==0')
end %if, the surface shear stress can not be zero
if (delta_z<z0) %in this case delta_z is set smaller than z0, which inhibits the calculation of the wind profile
    error('delta_z set too small/z0 set too large');
end %if, in this case delta_z is set smaller than z0, which inhibits the calculation of the wind profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Calculating the particle concentration and corresponding particle shear stress    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (initial_iteration || Owen) %Determining the particle concentration. Using Owen's hypothesis for the initial iteration and for subsequent iterations only as specified by the user. According to Owen's hypothesis, the surface shear stress should approximately equal the impact threshold
    n_fact_increase=rho*(max((~Owen)*0.05*u_fr^2,u_fr^2-u_fr_thr^2))/total_surface_shear_stress; %this is the multiplication factor for n_relative that will yield n (i.e., n_relative*n_fact_increase = n)
elseif (u_fr_thr_soil < u_fr_surf) %in this case, aerodynamic lifting occurs and the surface friction speed should equal the fluid threshold
    if (splash_fraction == 1) %if splash_fraction == 1, u_fr_surf might be larger than u_fr_thr_soil due to model fluctuations
        n_fact_increase=max(max(1,rho*(u_fr^2-u_fr_thr_soil^2)/total_surface_shear_stress),2^n_power*sum(m.*n_old)/sum(m.*n_relative)); %since splash_fraction == 1, u_fr_surf might be larger than u_fr_thr_soil due to model fluctuations. Thus, the number of lifted particles is set to the maximum of (1) the number of particles required to set the surface friction velocity to the threshold and (2) the number of particles as obtained from the splashing process; n_fact_increase is the multiplication factor for n_relative that will yield n (i.e., n_relative*n_fact_increase = n)
    else %in this case, aerodynamic entrainment is important and Owen's hypothesis applies
        n_fact_increase=max(1,rho*(u_fr^2-u_fr_thr_soil^2)/total_surface_shear_stress); %lifting enough particles for the surface friction speed to equal the fluid threshold; n_fact_increase is the multiplication factor for n_relative that will yield n (i.e., n_relative*n_fact_increase = n)
    end %if, if splash_fraction == 1, u_fr_surf might be larger than u_fr_thr_soil due to model fluctuations
else %in this case, splash lifting dominates and the model varies the surface friction speed to match the condition that, in steady-state, the number of splashed particles equals the number of particles lost to the soil bed.  Since the surface friction speed is below the fluid threshold, all particles are splashed, with no aerodynamic entrainment.
    n_power = max(min(0.4,n_power),-0.4); %limiting n_power between -0.4 and 0.4 to damp model oscillations
    if (total_surface_shear_stress < 0 || mean(particle_shear_stress)<0 || mean(particle_shear_stress_acc)<0) %these are all conditions that shouldn't occur and might mean something is wrong
        n_fact_increase = 0; %this is the multiplication factor for n_relative that will yield n (i.e., n_relative*n_fact_increase = n)
        n_negative = true; %n_negative is set to true when n < 0 (this can happen when the winds are much too low - having n < 0 will reduce the running average of n used in saltation_model and increase the wind profile towards steady state)
    else %this is the condition that should normally occur, and the particle concentration is determined from the splashing process
        n_fact_increase = sqrt(1+2*n_power)*sum(m.*n_old)/sum(m.*n_relative);  %this is the multiplication factor for n_relative that will yield n (i.e., n_relative*n_fact_increase = n)    end %if
    end %if
end %if, Determining the particle concentration. 
n = max(1,n_fact_increase*n_relative); %obtaining n from n_fact_increase. Note that this is handed back to saltation_model, where this is n absorbed into the running average from which the particle concentration used by saltation_model is derived
u_fr_surf_new=sign(u_fr^2-n_fact_increase*total_surface_shear_stress/rho)*sqrt(abs(u_fr^2-n_fact_increase*total_surface_shear_stress/rho)); %calculating the surface friction speed that would result from the n calculated in the previous line
if (u_fr_surf_new < -u_fr) %bounding n by the surface friction velocity that it produces. If it's too negative, it can cause model instability
    n_fact_increase = (rho/total_surface_shear_stress)*(u_fr^2+u_fr^2); %this is the multiplication factor for n_relative that will yield n (i.e., n_relative*n_fact_increase = n)
    u_fr_surf_new=sign(u_fr^2-n_fact_increase*total_surface_shear_stress/rho)*sqrt(abs(u_fr^2-n_fact_increase*total_surface_shear_stress/rho));
end %if, bounding n by the surface friction velocity that it produces. If it's too negative, it can cause model instability
tau_p = n_fact_increase*particle_shear_stress;  %tau_p is the total particle shear stress on the coarse grid
tau_p_acc = n_fact_increase*particle_shear_stress_acc; %tau_p is the total particle shear stress on the fine grid
tau_p_surface = n_fact_increase*total_surface_shear_stress; %tau_p is the total particle shear stress at the surface
tau_air = rho*u_fr^2-tau_p; %the fluid shear stress on the coarse grid. The fluid shear stress and the particle shear stress add up to the total shear stress (or momentum flux), which equals rho*u_fr^2
tau_air_acc = rho*u_fr^2-tau_p_acc; %the fluid shear stress on the fine grid. The fluid shear stress and the particle shear stress add up to the total shear stress (or momentum flux), which equals rho*u_fr^2
dudz = sign(tau_air).*sqrt(abs(tau_air)./(rho*0.4^2*z.^2)); %the slope of the wind profile, used to calculate the wind profile on the coarse grid. See description in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 (also at http://arxiv.org/abs/0901.0270 )
dudz_acc = sign(tau_air_acc).*sqrt(abs(tau_air_acc)./(rho*0.4^2*z_acc.^2)); %the slope of the wind profile, used to calculate the wind profile on the fine grid. See description in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 (also at http://arxiv.org/abs/0901.0270 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The below code subdivides the bottom two delta_z_acc grid points into many (~2^m-1) %%% 
%%% further grid points, in order to make the numerical integration of the wind profile %%%
%%% near the surface, which is very sensitive to grid size, more accurate               %%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=6; %this sets the smallest increment of z_near_surf to delta_z_acc/2^m-1
delta_z_wind = delta_z_acc; %the increment of z_near_surf
z_near_surf_temp(1) = 2*delta_z_acc; %the first grid point on z_near_surf is 2*delta_z_acc; later this is inverted so that this is the last (highest) grid point
p=1; %keeps track of how many times delta_z_wind has been halved
while (p<m)
    delta_z_wind = delta_z_wind/2;
    p=p+1;
    z_near_surf_temp(p) = z_near_surf_temp(p-1)-delta_z_wind; %subsequent grid points on z_near_surf_temp have sequentially halved grid box sizes
end
remainder_p = round(((z_near_surf_temp(p) - z0)/delta_z_wind)-0.5); %calculating how many increments of delta_z_wind are needed to fill the gap between z_near_surf_temp(p) and z0
if (remainder_p<2) %in this case there would not be enough values of remainder_p for proper numerical integration. This shouldn't occur and can be prevented by increasing delta_z_acc or m
    error('not enough z_near_surf values! Enlarge delta_z_acc or m');
end %if, in this case there would not be enough values of remainder_p for proper numerical integration. This shouldn't occur and can be prevented by increasing delta_z_acc or m
for p=m+1:1:m+remainder_p %filling the gap between z_near_surf_temp(p) and z0
    z_near_surf_temp(p)=z_near_surf_temp(p-1)-delta_z_wind;
end %for, z_near_surf_temp(p) and z0
z_near_surf(1:1:size(z_near_surf_temp,2)) = z_near_surf_temp(size(z_near_surf_temp,2):-1:1); %inverting z_near_surf_temp so that z_near_surf(1) is closest to z0 and z_near_surf(m+remainder_p) = 2*delta_z_acc
for p=1:1:size(z_near_surf,2) %matching z_near_surf to a value of tau_p_acc through linear interpolation
    j_near_surf_top = round(z_near_surf(p)/delta_z_acc+0.4999999); %calculating nearest grid point on the fine grid above z_near_surf(p)
    j_near_surf_bottom = j_near_surf_top-1; %calculating nearest grid point on the fine grid below z_near_surf(p)
    if (j_near_surf_bottom == 0) %calculating tau_p for the nearest grid point on the fine grid below z_near_surf(p)
        tau_p_bottom = tau_p_surface; %if j_near_surf_bottom == 0, then the nearest grid point is the surface
    else %in this case, the nearest grid point is not the surface
        tau_p_bottom = tau_p_acc(j_near_surf_bottom); %setting tau_p_bottom equal to the tau_p at the nearest grid point below z_near_surf(p)
    end %if, calculating tau_p for the nearest grid point on the fine grid below z_near_surf(p)
    if (j_near_surf_top > size(tau_p_acc,2)) %just checking that j_near_surf_top is not above tau_p_acc, which should never happen
        warning('something wrong with calculating very-near-surface wind profile')
    else
        tau_p_top = tau_p_acc(j_near_surf_top); %setting tau_p_bottom equal to the tau_p at the nearest grid point above z_near_surf(p)
    end %if, just checking that j_near_surf_top is not above tau_p_acc, which should never happen
    tau_p_near_surf(p) = ((j_near_surf_top*delta_z_acc-z_near_surf(p))/delta_z_acc)*tau_p_bottom+((z_near_surf(p)-j_near_surf_bottom*delta_z_acc)/delta_z_acc)*tau_p_top; %calculating the particle shear stress at z_near_surf(p) through linear interpolation between tau_p_bottom and tau_p_top
end %for, matching z_near_surf to a value of tau_p_acc through linear interpolation
tau_air_near_surf = rho*u_fr^2-tau_p_near_surf; %calculating the fluid shear stress on z_near_surf; The fluid shear stress and the particle shear stress add up to the total shear stress (or momentum flux), which equals rho*u_fr^2
dudz_near_surf = sign(tau_air_near_surf).*sqrt(abs(tau_air_near_surf)./(rho*0.4^2*z_near_surf.^2)); %calculating the wind profile slope. See description in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 (also at http://arxiv.org/abs/0901.0270 )
k1 = abs(z_near_surf-delta_z_acc)-min(abs(z_near_surf-delta_z_acc)); %this function is used to match to u_acc(1) and dudz_acc(1)
k2 = abs(z_near_surf-2*delta_z_acc)-min(abs(z_near_surf-2*delta_z_acc)); %this function is used to match to u_acc(2) and dudz_acc(2)
for p=1:1:size(z_near_surf,2) %retrieving u_near_surf from dudz_near_surf
    if (imag(dudz_near_surf(p))~=0) %checking for imaginary components, which should never occur
        warning('imaginary component of dudz_near_surf!')
        dudz_near_surf(p) = u_fr_thr/(0.4*z_near_surf(p)); %eliminating the imaginary component to prevent the model from crashing
    end %if, checking for imaginary components, which should never occur
    if (p==1) %retrieving u_near_surf(1) from dudz_near_surf(1)
        u_near_surf(1)=(z_near_surf(1)-z0)*dudz_near_surf(1); %calculating u_near_surf from dudz_near_surf
    elseif (p==2) %retrieving u_near_surf(2) from dudz_near_surf(2)
    	u_near_surf(2)=u_near_surf(1)+0.5*(z_near_surf(2)-z_near_surf(1))*(dudz_near_surf(2)+dudz_near_surf(1)); %calculating u_near_surf from dudz_near_surf
    elseif (p==3) %retrieving u_near_surf(3) from dudz_near_surf(3)
    	u_near_surf(3)=u_near_surf(2)+((z_near_surf(3)-z_near_surf(2))/12)*(-dudz_near_surf(1)+8*dudz_near_surf(2)+5*dudz_near_surf(3)); %calculating u_near_surf from dudz_near_surf
    elseif ((z_near_surf(p)-z_near_surf(p-1)) == delta_z_wind) %in this case p>3 and the difference between successive z_near_surf values is not constant such the Adams-Moulton scheme can be used
        u_near_surf(p)=u_near_surf(p-1)+((z_near_surf(p)-z_near_surf(p-1))/24)*(dudz_near_surf(p-3)-5*dudz_near_surf(p-2)+19*dudz_near_surf(p-1)+9*dudz_near_surf(p)); %calculating u_near_surf from dudz_near_surf
    else %in this case the difference between successive z_near_surf values is not constant
        del_p = round(log((z_near_surf(p)-z_near_surf(p-1))/delta_z_wind)/log(2)); %this corrects for the variation in grid box size such that z_near_surf(p-2-del_p) - z_near_surf(p-1) = z_near_surf(p) - z_near_surf(p-1)
        u_near_surf(p)=u_near_surf(p-1)+((z_near_surf(p)-z_near_surf(p-1))/12)*(-dudz_near_surf(p-2-del_p)+8*dudz_near_surf(p-1)+5*dudz_near_surf(p)); %calculating u_near_surf from dudz_near_surf
    end %if-else, retrieving u_near_surf from dudz_near_surf
    if (k1(p)==0) %matching u_acc(1) and dudz_acc(1) to the corresponding values on z_near_surf
        u_acc(1) = u_near_surf(p); %matching u_acc(1)
        dudz_acc(1) = dudz_near_surf(p); %matching dudz_acc(1)
    end %if, matching u_acc(1) and dudz_acc(1) to the corresponding values on z_near_surf
    if (k2(p)==0) %if, matching u_acc(2) and dudz_acc(2) to the corresponding values on z_near_surf
        u_acc(2) = u_near_surf(p); %matching u_acc(2)
        dudz_acc(2) = dudz_near_surf(p); %matching dudz_acc(2)
    end %if, matching u_acc(1) and dudz_acc(1) to the corresponding values on z_near_surf
end %for, retrieving u_near_surf from dudz_near_surf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Calculating u_acc and u         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j_acc=3:j_acc_end %applying the Adam-Bashforth-Moulton Predictor-Corrector Method to obtain u_acc from dudz_acc now that dudz_acc(1) and dudz_acc(2) are known
    if (abs(imag(dudz_acc(j_acc)))>0) %checking for imaginary components, which should never occur
    	warning('imaginary component of dudz_acc wind profile!')
        dudz_acc(j_acc)=u_fr_thr/(0.4*z_acc(j_acc)); %eliminating the imaginary component to prevent the model from crashing
    end %if, checking for imaginary components, which should never occur
    if (j_acc==3) %slightly different procedure for j_acc==3, where fewer data points are available
        u_acc(j_acc) = u_acc(j_acc-1)+((z_acc(j_acc)-z_acc(j_acc-1))/12)*(-dudz_acc(j_acc-2)+8*dudz_acc(j_acc-1)+5*dudz_acc(j_acc)); %calculating u_acc from dudz_acc
    else %normal Adams-Moulton scheme
        u_acc(j_acc) = u_acc(j_acc-1)+((z_acc(j_acc)-z_acc(j_acc-1))/24)*(dudz_acc(j_acc-3)-5*dudz_acc(j_acc-2)+19*dudz_acc(j_acc-1)+9*dudz_acc(j_acc)); %calculating u_acc from dudz_acc
    end %if, slightly different procedure for j_acc==3, where fewer data points are available
    if (rem(j_acc*delta_z_acc,delta_z)==0) %matching c(1:no_fine_grid) and dudz(1:no_fine_grid) to u_acc and dudz_acc
        u(round(j_acc*delta_z_acc/delta_z)) = u_acc(j_acc); %matching u to u_acc
        dudz(round(j_acc*delta_z_acc/delta_z)) = dudz_acc(j_acc); %matching dudz to dudz_acc
    end %if, matching c(1:no_fine_grid) and dudz(1:no_fine_grid) to u_acc and dudz_acc
end %for, applying the Adam-Bashforth-Moulton Predictor-Corrector Method to obtain u_acc from dudz_acc now that dudz_acc(1) and dudz_acc(2) are known
for j=no_fine_grid+1:j_end %applying the Adam-Bashforth-Moulton Predictor-Corrector Method	
    if (abs(imag(dudz(j)))>0) %checking for imaginary components, which should never occur
    	warning('imaginary component of dudz wind profile!')
        dudz(j)=u_fr_thr/(0.4*z(j)); %eliminating the imaginary component to prevent the model from crashing
    end %if, checking for imaginary components, which should never occur
    u(j) = u(j-1)+((z(j)-z(j-1))/24)*(dudz(j-3)-5*dudz(j-2)+19*dudz(j-1)+9*dudz(j)); %calculating u from dudz
end %for, applying the Adam-Bashforth-Moulton Predictor-Corrector Method	
if (n_negative == true) %this is to reduce the instability caused by getting n < 0 (this can happen when the winds are much too low - having n < 0 will reduce the running average of n used in saltation_model and increase the wind profile towards steady state)
    n = n_old/2;
end %if, this is to reduce the instability caused by getting n < 0 (this can happen when the winds are too low)