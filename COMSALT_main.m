function [avg_z50, avg_z50_mass, avg_vI, avg_thetaI, avg_vx, avg_z0s, avg_ufr_surf, avg_log_ejecta, avg_mass, avg_massflux, avg_length, avg_delVx, total_iterations] = COMSALT_main (t_max_temp, h_temp, delta_z_temp, u_fr_temp, psd_flag_temp, u_fr_thr_temp, filetext_temp, switches_temp, spare_temp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          Defining global parameters                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global u_fr u_fr_thr u_fr_thr_soil z0 delta_z turbulence delta_z_acc no_fine_grid j_acc_end elevation cohesion;
global rho_p pi S P d_thr a_splash Saffman j_end part_no h t_max d_mean major_iteration splash_fraction Owen spare g rho mu Mars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Processing the input of the function call into some of the global parameters    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_max = t_max_temp  %the total time of simulated orbits of saltating particles in a given particle for a single iteration of saltation_trajectory (should be set to >= 10 seconds for Earth ambient conditions)
h = h_temp   %the time step in seconds with which the saltating particles' motion is modeled (should be set to <= 0.075 seconds for Earth ambient conditions)
u_fr = u_fr_temp  %the wind friction or shear velocity in m/s (square root of the ratio of wind shear stress to air density)
delta_z = delta_z_temp %the vertical grid size in meters
switches = switches_temp; %this contains switches that set Owen's 2nd hypothesis, turbulence, etc.
psd_flag = psd_flag_temp; %this flag sets the particle size distribution
u_fr_thr = u_fr_thr_temp; %the impact threshold friction velocity required to sustain saltation. This needs to be supplied to the model, but can be calculated for a given soil size distribution using the function findImpThr
filetext = filetext_temp; %this is the text that's appended to the filename of output files
spare = spare_temp; %this 'spare' variable can be used to diagnoze model problems

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Setting 'switches'        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(switches,2)~=2) %the function call should contain 2 switches
    error('too few or too many switches!');
end %if
Owen = switches(1); %The parameter 'Owen' sets whether the model should employ Owen's second hypothesis, i.e., that the surface shear stress remains at the threshold for particle entrainment. This switch should be set to zero unless the model operates in conditions where Owen's hypothesis is probably valid (on Venus, for example)
turbulence = switches(2); %This parameters sets whether the effects of turbulence on the saltation trajectories is simulated. See description in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 (also at http://arxiv.org/abs/0901.0270 )
Saffman = false; %This switch sets whether the Saffman force is calculated, which is computationally intensive and has little effect on trajectories. See description in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 (also at http://arxiv.org/abs/0901.0270 ) 
Mars = false; %If this switch is set to true, it will set gravity and thermodynamic conditions to those at the Martian surface
cohesion = false; %If this switch is set to true, the model will account for the effects of cohesion on the splashing of soil particles, following Kok, GRL, 2010
plot_wind_profile = false; %this sets whether the wind profile is plotted after every major iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Initializing numerical parameters    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_major_iterations = 4;  %the number of iterations of the wind profile, saltation trajectories, mass flux, and particle concentration that are used to determine convergence of the model. Two of these iterations are used to initiate the model.
total_minor_iterations = 25;  %the number of iterations with which the wind profile and saltation trajectories are determined. The loop over wind-profile iterations is contained within the loop over total_major_iterations, such that the total number of iterations equals (total_major_iterations+2)*(total_minor_iterations-carry_minor_iterations) + carry_minor_iterations
carry_minor_iterations = 8; %the number of wind iterations that's carried over to the next iteration. This is done to provide improved stability between iterations
initial_minor_iterations = 1;  %The number of initial wind iterations that are used to initialize the saltation trajectories (should normally be set to 1)
tolerance = 0.15; %the tolerance within which essential parameters (z50 and mass flux) must be simulated for the model to converge. Dependent on conditions and desired model accuracy, this should normally be set to 0.02 to 0.20
hor_interval = 0.01; %the horizontal interval used in determining the horizontal mass flux (important for comparison with measurements by Namikas, Sedimentology, 2003)
delta_z_acc = delta_z/25; %delta_z_acc the spacing of the fine grid that is used to calculate the wind profile close to the surface, where higher resolution is needed
no_fine_grid = 5; %the size of the fine grid in units of the coarser grid spacing
if (no_fine_grid < 3 || round(no_fine_grid) ~= no_fine_grid) %no_fine_grid must be an integer and must equal at least 3
    error('Increase no_fine_grid to at least 3 and make it an integer)')
end %no_fine_grid must be an integer and must equal at least 3
j_end = no_fine_grid + 1; %the total number of vertical layers on the coarse grid. j_end is automatically increased if particles exceed j_end*delta_z
j_acc_end = no_fine_grid*round(delta_z/delta_z_acc); %the total number of vertical layers on the fine grid.

%////////////////////////////////////////////////////////////////////////////////////
%   This portion of the code runs initializing functions related to the
%   particle   //
%   size distribution, soil characteristics, and other physical parameters         //
%////////////////////////////////////////////////////////////////////////////////////
d_test = 0; %use this function to test particle particle sizes within the 'default' soil particle size distribution of Namikas, Sedimentology, 2003
if (psd_flag==1) %this sets the particle size distribution, d_temp = 1 is a flag for the size distribution described in Namikas, Sedimentology, 2003
    [part_no, d, d_mean, mass_fraction] = PSD_Namikas (d_test); %this sets the soil size distribution
    d_thr = d_mean; %the particle diameter in m that determines the threshold friction velocity
    elevation = 0; %elevation in m
    z0 = d_mean/30 %the aerodynamic surface roughness (Nikuradse, 1933)
elseif (psd_flag==2) %this sets the particle size distribution, d_temp = 2 is a flag for the size distribution described in Bagnold, Proceedings of the Royal Society of London Series A - Mathematical and Physical Sciences, 1938
    [part_no, d, d_mean, mass_fraction] = PSD_Bagnold; %this sets the soil size distribution
    d_thr = d_mean; %the particle diameter in m that determines the threshold friction velocity
    elevation = 0; %elevation in m
    z0 = d_mean/30 %the aerodynamic surface roughness (Nikuradse, 1933)
elseif (psd_flag==3) %this sets the particle size distribution, d_temp = 3 is a flag for a size distribution at the crest of a dune described in Lancaster, Journal of Sedimentary Petrology, 1986
    [part_no, d, d_mean, mass_fraction] = PSD_Dune; %this sets the soil size distribution
    elevation = 0; %elevation in m
    d_thr = d_mean; %the particle diameter in m that determines the threshold friction velocity
    z0 = d_mean/30 %the aerodynamic surface roughness (Nikuradse, 1933)
elseif (psd_flag==4) %this sets the particle size distribution, d_temp = 4 is a flag for size distribution described in Williams, Sedimentology, 1964
    [part_no, d, d_mean, mass_fraction] = PSD_Williams; %this sets the soil size distribution
    elevation = 0; %elevation in m
    d_thr = d_mean; %the particle diameter in m that determines the threshold friction velocity    
    z0 = d_mean/30 %the aerodynamic surface roughness (Nikuradse, 1933)
else %in this case the user set a monodisperse size distribution, with particle size equal to psd_flag
    mass_fraction(1) = 1; %the volume fraction of the soil occupied by a particular particle bin
    part_no = 1;
    d = psd_flag; d_thr = psd_flag; d_mean = psd_flag; %the particle diameter in m that determines the threshold friction velocity    
    elevation = 0; %elevation in m
    z0 = d_mean/30 %the aerodynamic surface roughness (Nikuradse, 1933)
end %if, this sets the particle size distribution
if (delta_z<z0) %for the model to function, z0 must be smaller than delta_z
    error('delta_z set too small/z0 set too large');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Loading all the physical parameters        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_parameters(d); %this function contains many of the physical parameters that are global variables, such as particle density and parameters for the turbulence and splash function parameterizations (see description in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 (also at http://arxiv.org/abs/0901.0270 ) 
m = (1/6)*pi*rho_p*d.^3; %particle mass in kg
if (Mars == true) %this sets the conditions for a simulation for Martian conditions
    rho_p = 3000; %particle density in kg/m3, after Claudin and Andreotti, Earth and Planetary Science Letters, 2006
    g = 3.71; %gravitational acceleration in m/s2
    T = 227; %Mars surface temperature in K
    m_molar = 0.0435; %molar weight of a mol of air
    P0 = 627;  %Sea-level pressure in Pa
    scale_height = 8.31*T/(m_molar*g); %the scale height of the atmosphere in m
    P = P0*exp(-elevation/scale_height);  %the pressure in Pascal
    rho = P/(8.31*T)*m_molar; %air density in kg/m3 from the ideal gas l    
    mu = 0.00751*(T/291.15)^(1.5)/(T+120); %the viscosity in kg/(m s) from the Sutherland relation
    u_fr_thr_soil = 0.10*sqrt(((rho_p-rho)/rho)*g*d_mean); %the fluid threshold
end %if, this sets the conditions for a simulation for Martian conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This portion of the code initializes various other parameters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convergence = false; %when the model converges, this will be set to true and the while (convergence == false) loop will terminate
particle_coll_freq = zeros(part_no, part_no); %the frequency of collisions between saltating particles that occur in mid-air
v_x_avg_down = zeros(part_no, j_end); v_x_avg_up = zeros(part_no, j_end); %initializing the matrix that holds the average horizontal speed for respectively downward and upward-moving particles
v_z_avg_down = zeros(part_no, j_end); v_z_avg_up = zeros(part_no, j_end); %initializing the matrix that holds the average vertical speed for respectively downward and upward-moving particles
density_down=zeros(part_no,j_end); density_up=zeros(part_no,j_end); %initializing the concentration matrix, describing the likelihood that any given particle is at a particular height, for respectively downward and upward-moving particles. Note that density_down and density_up sum to 1 for a particular particle bin.
ejecta_totals = zeros(part_no,1); %the total number of surface particles ejected by saltating particles in a particular particle bin
splash_fraction = 0; %the fraction of new saltating particles that is splashed - this is only important when the fluid threshold is equal to or below the impact threshold and only occurs in high-density fluids (water, Venus)
pEj = zeros(part_no,1); %this holds the speeds of the ejected particles of particular particle bins. This information is used in saltation_trajectory when a new saltating particle is generated.
weighting_factor_part = 10*t_max*ones(part_no,1)'; %initializing weighting_factor_part by setting it to a large number
major_iteration = 0; %this holds the major_iteration number.
z_plot = delta_z:delta_z:j_end*delta_z; z_plot_acc = delta_z_acc:delta_z_acc:j_acc_end*delta_z_acc; %the vertical position arrays for plotting various parameters, such as the wind profile
u = (min(u_fr,u_fr_thr)/0.4)*log(z_plot/z0); u_acc = (min(u_fr,u_fr_thr)/0.4)*log(z_plot_acc/z0); %this initializes the wind profile
n = zeros(1,part_no); %the particle concentration in #/m2
z_s = 0.20; %initial height of saltation layer (used to calculate Lagrangian time scale)
cum_minor_iteration = 0; %the total number of minor iterations over all major iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This portion of the code iteratively simulates the trajectories of saltating particles, the   %%%
% ejection of surface particles, the particle concentration, and the wind profile in saltation, %%%
% until convergence is achieved.                                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (convergence == false) %iterating until convergence is achieved
    major_iteration = major_iteration + 1 %keeping track of the major_iteration number
    if (major_iteration > 1)
        %the below carries all the variables from the final minor iteration of the previous major_iteration into the new iterations. This provides stability between subsequent iterations.
        start_position = start_position_carry; %start_position keeps track of the starting position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
        time_end_total = time_end_total_carry; %keeping track of the total time that trajectories of a particular particle bin are simulated. Used in calculating average quantities.
        surface_coll_freq_total = surface_coll_freq_total_carry; no_rebound_prob_total = no_rebound_prob_total_carry; %the frequency per second of surface collicions, and the average probability that an impacting particle will not rebound
        vI_total = vI_total_carry; v_xI_total = v_xI_total_carry; v_zI_total = v_xI_total_carry; vR_total = vR_total_carry; v_xR_total = v_xI_total_carry; v_zR_total = v_xI_total_carry; %parameters that hold the total, x, and z-direction speeds for impacting and rebounding particles
        vI_distr_total = vI_distr_total_carry; v_lift_total = v_lift_total_carry; %parameters that hold the average impacting speed with which surface particle are ejected, the average lifting speed of splashed particles
        thetaI_total = thetaI_total_carry; thetaR_total = thetaR_total_carry; %the average angle with the horizontal of impacting and rebounding particles
        v_x_avg_total = v_x_avg_total_carry; v_x_avg_acc_total = v_x_avg_acc_total_carry; %the average horizontal speeds of particles in the coarse and fine grids
        bounce_length_total = bounce_length_total_carry; height_total = height_total_carry; length_total = length_total_carry; alt_length_total = alt_length_total_carry; %parameters that hold approximate measures of the length of individual bounces, the average height and length of bounces per particle bin, and an alternative measure of the length of saltation trajectories, weighted by the impact speed
        z_append = z_append_carry; t_append = t_append_carry; v_x_append = v_x_append_carry; v_z_append = v_z_append_carry; %parameters that help reconstruct z, t, v_x, and v_z after a major iteration
        density_down_total = density_down_total_carry; density_up_total = density_up_total_carry; density_down_acc_total = density_down_acc_total_carry; density_up_acc_total = density_up_acc_total_carry; %parameters that help reconstruct the density of upward and downward moving particles after a major iteration        
        aerosol_emis_tot = aerosol_emis_tot_carry; no_bounces_total = no_bounces_total_carry; %aerosol_emis_tot is a measure of dust aerosol emission, calculated from the total kinetic energy impacting the soil surface, and no_bounces_total holds the total number of bounces of particles in a particular particle bin
        counter_down_total(:,:,1:carry_minor_iterations) = counter_down_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these counter variables are calculated in calc_avg_velocity_density and keep track of the total amount of time spent by upward and downward moving particles in particular particle bins. It's used to calculate the average particle velocity and density for each grid points after a major iteration
        counter_down_acc_total(:,:,1:carry_minor_iterations) = counter_down_acc_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these counter variables are calculated in calc_avg_velocity_density and keep track of the total amount of time spent by upward and downward moving particles in particular particle bins. It's used to calculate the average particle velocity and density for each grid points after a major iteration
        counter_up_total(:,:,1:carry_minor_iterations) = counter_up_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these counter variables are calculated in calc_avg_velocity_density and keep track of the total amount of time spent by upward and downward moving particles in particular particle bins. It's used to calculate the average particle velocity and density for each grid points after a major iteration
        counter_up_acc_total(:,:,1:carry_minor_iterations) = counter_up_acc_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these counter variables are calculated in calc_avg_velocity_density and keep track of the total amount of time spent by upward and downward moving particles in particular particle bins. It's used to calculate the average particle velocity and density for each grid points after a major iteration
        v_x_avg_up_total(:,:,1:carry_minor_iterations) = v_x_avg_up_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these parameters are used to calculate the total average particle speeds after a major iteration
        v_x_avg_up_acc_total(:,:,1:carry_minor_iterations) = v_x_avg_up_acc_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these parameters are used to calculate the total average particle speeds after a major iteration
        v_z_avg_up_total(:,:,1:carry_minor_iterations) = v_z_avg_up_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these parameters are used to calculate the total average particle speeds after a major iteration
        v_x_avg_down_total(:,:,1:carry_minor_iterations) = v_x_avg_down_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these parameters are used to calculate the total average particle speeds after a major iteration
        v_x_avg_down_acc_total(:,:,1:carry_minor_iterations) = v_x_avg_down_acc_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these parameters are used to calculate the total average particle speeds after a major iteration
        v_z_avg_down_total(:,:,1:carry_minor_iterations) = v_z_avg_down_total(:,:,total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); %these parameters are used to calculate the total average particle speeds after a major iteration
        weighting_factor_temp = weighting_factor; %creating temporary functions to carry over several variables
        tau_p_surface_total_temp = tau_p_surface_total; u_fr_surf_save_temp = u_fr_surf_save; n_power_save_temp = n_power_save; ejecta_orbit_ratio_temp = ejecta_orbit_ratio; %creating temporary functions to carry over several variables
        weighting_factor(1:carry_minor_iterations) = weighting_factor_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); weighting_factor(carry_minor_iterations+1:total_minor_iterations) = weighting_factor_temp(1:total_minor_iterations-carry_minor_iterations); %resetting weighting_factor; weighting_factor(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in weighting_factor
        tau_p_surface_total(1:carry_minor_iterations) = tau_p_surface_total_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); tau_p_surface(carry_minor_iterations+1:total_minor_iterations) = tau_p_surface_total_temp(1:total_minor_iterations-carry_minor_iterations); %resetting tau_p_surface_total; tau_p_surface_total(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in tau_p_surface_total
        u_fr_surf_save(1:carry_minor_iterations) = u_fr_surf_save_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); u_fr_surf_save(carry_minor_iterations+1:total_minor_iterations) = u_fr_surf_save_temp(1:total_minor_iterations-carry_minor_iterations); %resetting u_fr_surf_save; u_fr_surf_save(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in u_fr_surf_save
        n_power_save(1:carry_minor_iterations) = n_power_save_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); n_power_save(carry_minor_iterations+1:total_minor_iterations) = n_power_save_temp(1:total_minor_iterations-carry_minor_iterations); %resetting n_power_save; n_power_save(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in n_power_save
        ejecta_orbit_ratio(1:carry_minor_iterations) = ejecta_orbit_ratio_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations); ejecta_orbit_ratio(carry_minor_iterations+1:total_minor_iterations) = ejecta_orbit_ratio_temp(1:total_minor_iterations-carry_minor_iterations); %resetting ejecta_orbit_ratio; ejecta_orbit_ratio(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in ejecta_orbit_ratio
        impacts_temp = impacts; rebounds_ejecta_temp = rebounds_ejecta; lost_temp = lost; created_temp = created;
        impacts(1:carry_minor_iterations,:) = impacts_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations,:); impacts(carry_minor_iterations+1:total_minor_iterations,:) = impacts_temp(1:total_minor_iterations-carry_minor_iterations,:); %resetting impacts; impacts(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in impacts
        rebounds_ejecta(1:carry_minor_iterations,:) = rebounds_ejecta_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations,:); rebounds_ejecta(carry_minor_iterations+1:total_minor_iterations,:) = rebounds_ejecta_temp(1:total_minor_iterations-carry_minor_iterations,:); %resetting rebounds_ejecta; rebounds_ejecta(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in rebounds_ejecta
        lost(1:carry_minor_iterations,:) = lost_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations,:); lost(carry_minor_iterations+1:total_minor_iterations,:) = lost_temp(1:total_minor_iterations-carry_minor_iterations,:); %resetting lost; lost(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in lost
        created(1:carry_minor_iterations,:) = created_temp(total_minor_iterations-carry_minor_iterations+1:total_minor_iterations,:); created(carry_minor_iterations+1:total_minor_iterations,:) = created_temp(1:total_minor_iterations-carry_minor_iterations,:); %resetting created; created(1:carry_minor_iterations) holds the last carry_minor_iterations number of minor iterations, and so the total_minor_iteration-carry_minor_iterations that were done before are now in carry_minor_iterations+1:total_iterations. In this way, new minor iterations will be written over the oldest minor iterations stored in created
        initial_minor_iterations = 0; %setting initial_minor_iterations to 0 after the first major iteration
        minor_iteration = carry_minor_iterations; %resetting minor_iteration to carry_minor_iterations
    else %for the first iteration, there are no variables to carry over
        minor_iteration = 0; %setting minor_iteration to 0 for the first major iteration
        u_fr_surf = u_fr_thr; %setting the surface friction velocity to the impact threshold for the first major iteration
    end %if, the below carries all the variables from the final minor iterations of the previous iteration into the new iterations. This provides stability between subsequent iterations
    while (minor_iteration<total_minor_iterations) %performing a number of iterations while continuously adjusting the wind profile and (re)calculating particle trajectories
        minor_iteration = minor_iteration+1 %increasing the minor_iteration by 1
        cum_minor_iteration = cum_minor_iteration+1; %the total number of minor iterations over all major iterations        
        z = 0; %this resets the particle height for the while loop on the next line
        while (max(max(z')-d) < delta_z_acc) %the while loop prevents the exceedingly rare case that no particle leaves the lowest grid box from crashing the model. If this occurs repeatedly, the size of the fine grid box (delta_z_acc) is likely too large
            %the below function computes the saltating particle trajectories
            [t, x, z, v_x, v_z, time_end, i_end, surface_coll_freq, no_rebound_prob, height, length, alt_length, v_xI, v_zI, vI, thetaI, v_xR, v_zR, vR, thetaR, vI_distr, surface_shear_stress, bounce_length, no_bounces, n_relative, created(minor_iteration,1:part_no), lost(minor_iteration,1:part_no), impacts(minor_iteration,1:part_no), rebounds_ejecta(minor_iteration,1:part_no), ejecta_orbit_ratio(minor_iteration), weighting_factor(minor_iteration), weighting_factor_part(minor_iteration,:), pEj_temp, ejecta_totals_temp, v_lift, aerosol_emis, j_end_prev] = saltation_trajectory (d, u, u_acc, pEj, n, mass_fraction, particle_coll_freq, v_x_avg_down, v_x_avg_up, v_z_avg_down, v_z_avg_up, density_down, density_up, mean(weighting_factor_part,1), u_fr_surf(max(1,major_iteration-1)), z_s(max(1,major_iteration-1)));
        end %while, the while loop prevents the exceedingly rare case that no particle leaves the lowest grid box from crashing the model
        for k=1:1:part_no %adding pEj_temp to the distribution of splashed particle speeds pEj
            if (ejecta_totals_temp(k) > 0) %in this case there are particle ejecta to add to pEj
                pEj(k,ejecta_totals(k)+1:ejecta_totals(k)+ejecta_totals_temp(k)) = pEj_temp(k,1:ejecta_totals_temp(k)); %adding the result of the last calculation pEj_temp of particle trajectories by saltation_trajectory to pEj
                ejecta_totals(k) = ejecta_totals(k) + ejecta_totals_temp(k); %keeping track of the total number of splashed particles per particle bin
            end %if, in this case there are particle ejecta to add to pEj
        end %for, adding pEj_temp to the distribution of splashed particle speeds pEj        
        splash_fraction = min(1,sum(ejecta_orbit_ratio(max(1,minor_iteration-4):minor_iteration)*weighting_factor(max(1,minor_iteration-4):minor_iteration)')/sum(weighting_factor((max(1,minor_iteration-4):minor_iteration)))); %splash_fraction calculates which fraction of particle ejecta is from splash entrainment. This is only used when particle lifting is not dominated by splash entrainment, which can occur due to strong vertical fluid drag (on Venus or under water)
        
        %The following computes average velocities from the above calculated trajectories and computes the normalized density of  particles as a function of height.  The integral over height of the density is identically 1.
        [v_x_avg, v_x_avg_acc, counter_down, counter_up, counter_down_acc, counter_up_acc, v_x_avg_up_temp, v_x_avg_down_temp, v_z_avg_up_temp, v_z_avg_down_temp, density_up_temp, density_down_temp, flux_down, wind_total_mass_flux(minor_iteration), flux_down_acc, v_x_avg_down_acc_temp, v_x_avg_up_acc_temp, v_z_avg_down_acc_temp, v_z_avg_up_acc_temp, density_up_acc_temp, density_down_acc_temp] = calc_avg_velocity_density (d, x, z, t, v_x, v_z, time_end, i_end, n);
        
        %The following calculates the particle shear stress due to the absorption of wind momentum by the saltating particles
        [total_surface_shear_stress, particle_shear_stress, particle_shear_stress_acc] = calc_part_shear_stress (n_relative, d, flux_down, v_x_avg_up_temp, v_x_avg_down_temp, flux_down_acc, v_x_avg_down_acc_temp, v_x_avg_up_acc_temp);
        
        %This loop provides a running average of all important saltation variables, after a certain number of initial iterations (usually 1) have been performed to initialize the run
        if(minor_iteration>initial_minor_iterations) %the model must equilibriate first (during the initial_minor_iteration(s))
            counter_up_total(1:part_no,1:j_end,minor_iteration) = counter_up; %these counter variables are calculated in calc_avg_velocity_density and keep track of the total amount of time spent by upward and downward moving particles in particular particle bins. It's used to calculate the average particle velocity and density for each grid points after a major iteration
            counter_down_total(1:part_no,1:j_end,minor_iteration) = counter_down; 
            counter_up_acc_total(1:part_no,1:j_acc_end,minor_iteration) = counter_up_acc;
            counter_down_acc_total(1:part_no,1:j_acc_end,minor_iteration) = counter_down_acc;            
            v_x_avg_up_total(1:part_no,1:j_end,minor_iteration) = counter_up.*v_x_avg_up_temp; %these parameters are used to calculate the total average particle speeds after a major iteration
            v_x_avg_up_acc_total(1:part_no,1:j_acc_end,minor_iteration) = counter_up_acc.*v_x_avg_up_acc_temp;
            v_z_avg_up_total(1:part_no,1:j_end,minor_iteration) = counter_up.*v_z_avg_up_temp; 
            v_z_avg_up_acc_total(1:part_no,1:j_acc_end,minor_iteration) = counter_up_acc.*v_z_avg_up_acc_temp; 
            v_x_avg_down_total(1:part_no,1:j_end,minor_iteration) = counter_down.*v_x_avg_down_temp;
            v_x_avg_down_acc_total(1:part_no,1:j_acc_end,minor_iteration) = counter_down_acc.*v_x_avg_down_acc_temp;
            v_z_avg_down_total(1:part_no,1:j_end,minor_iteration) = counter_down.*v_z_avg_down_temp;            
            v_z_avg_down_acc_total(1:part_no,1:j_acc_end,minor_iteration) = counter_down_acc.*v_z_avg_down_acc_temp;                        
            if(minor_iteration==initial_minor_iterations+1 && major_iteration==1) %checking whether this is the first iteration to be taken into the running average             
                start_position = ones(part_no,1); %start_position keeps track of the starting position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
                no_bounces_total = no_bounces; %keeping track of the total number particle hops or bounces
                end_position = i_end;  %end_position keeps track of the ending position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
                time_end_total = time_end; %keeping track of the total time that trajectories of a particular particle bin are simulated. Used in calculating average quantities.
                surface_coll_freq_total = surface_coll_freq; %the frequency per second of surface collicions
                v_xI_total = v_xI*weighting_factor(minor_iteration); %this weighs the impacting and rebounding particle speed by the total number of saltation trajectories simulated in an iteration. The "total" variable is later divided by the sum of weight_factor to obtain the average particle speed
                v_zI_total = v_zI*weighting_factor(minor_iteration);
                v_xR_total = v_xR*weighting_factor(minor_iteration);
                v_zR_total = v_zR*weighting_factor(minor_iteration);                
                vI_distr_total = vI_distr; %the average impacting speed with which surface particle are ejected
                v_lift_total = v_lift; %the average speed with which splashed particles are ejected
                vI_total = vI; vR_total = vR; %the average impacting and rebounding particle speed
                thetaI_total = thetaI; thetaR_total = thetaR; %the average angles with the horizontal for impacting and rebounding particles
                v_x_avg_total = v_x_avg*time_end; %the average horizontal particle speed on the coarse grid
                v_x_avg_acc_total = v_x_avg_acc*time_end; %the average horizontal particle speed on the fine grid
                bounce_length_total = bounce_length; %this holds the horizontal length of individual saltating particle hops
                no_rebound_prob_total = no_rebound_prob.*time_end; %the average probability that an impacting particle will not rebound
                height_total = height'.*time_end; %this holds a measure of the average height of saltating particle trajectories
                length_total = length'.*time_end; %this holds a measure of the average length of saltating particle trajectories
                alt_length_total = alt_length'.*time_end; %this holds an alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed
                for k=1:1:part_no %cycling over all particle size bins
                    z_append(k,round(start_position(k)):1:round(end_position(k))) = z(k,1:1:round(end_position(k))); %saving z by appending it to z_append
                    t_append(k,round(start_position(k)):1:round(end_position(k))) = t(k,1:1:round(end_position(k))); %saving t by appending it to t_append
                    v_x_append(k,round(start_position(k)):1:round(end_position(k))) = v_x(k,1:1:round(end_position(k))); %saving v_x by appending it to v_x_append
                    v_z_append(k,round(start_position(k)):1:round(end_position(k))) = v_z(k,1:1:round(end_position(k))); %saving v_z by appending it to v_z_append
                    density_down_total(k,1:j_end) = density_down_temp(k,1:j_end)*time_end(k); density_down_acc_total(k,:) = density_down_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of downward moving particles on the coarse and fine grids
                    density_up_total(k,1:j_end) = density_up_temp(k,1:j_end)*time_end(k); density_up_acc_total(k,:) = density_up_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of upward moving particles on the coarse and fine grids
                end %for, cycling over all size particle size bins
                aerosol_emis_tot = aerosol_emis; %aerosol_emis is a measure of dust aerosol emission, calculated from the total kinetic energy impacting the soil surface, and no_bounces_total holds the total number of bounces of particles in a particular particle bin
            else %this is not the first iteration to be taken into the running average
                end_position = start_position + i_end - 1; %end_position keeps track of the ending position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
                time_end_total = time_end + time_end_total; %keeping track of the total time that trajectories of a particular particle bin are simulated. Used in calculating average quantities.
                surface_coll_freq_total = surface_coll_freq_total + surface_coll_freq; %the frequency per second of surface collicions
                v_xI_total = v_xI_total + v_xI*weighting_factor(minor_iteration); %this weighs the impacting and rebounding particle speed by the total number of saltation trajectories simulated in an iteration. The "total" variable is later divided by the sum of weight_factor to obtain the average particle speed
                v_zI_total = v_zI_total + v_zI*weighting_factor(minor_iteration);
                v_xR_total = v_xR_total + v_xR*weighting_factor(minor_iteration);
                v_zR_total = v_zR_total + v_zR*weighting_factor(minor_iteration);
                v_lift_total(size(v_lift_total,2)+1:1:size(v_lift_total,2)+size(v_lift,2)) = v_lift; %the average speed with which splashed particles are ejected
                vI_distr_total(size(vI_distr_total,2)+1:1:size(vI_distr_total,2)+size(vI_distr,2)) = vI_distr; %the average impacting speed with which surface particle are ejected
                v_x_avg_total = v_x_avg_total + v_x_avg'.*time_end; %the average horizontal particle speed on the coarse grid
                v_x_avg_acc_total = v_x_avg_acc_total + v_x_avg_acc'.*time_end; %the average horizontal particle speed on the fine grid
                no_rebound_prob_total = no_rebound_prob_total + no_rebound_prob.*time_end; %the average probability that an impacting particle will not rebound
                height_total = height'.*time_end + height_total; %this holds a measure of the average height of saltating particle trajectories
                length_total = length'.*time_end + length_total; %this holds a measure of the average length of saltating particle trajectories
                alt_length_total = alt_length'.*time_end + alt_length_total; %this holds an alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed
                for k=1:1:part_no  %cycling over all particle size bins                   
                    z_append(k,round(start_position(k)):1:round(end_position(k))) = z(k,1:1:i_end(k)); %saving z by appending it to z_append                
                    t_append(k,round(start_position(k)):1:round(end_position(k))) = t(k,1:1:i_end(k)) + t_append(k,round(start_position(k))-1); %saving t by appending it to t_append
                    v_x_append(k,round(start_position(k)):1:round(end_position(k)))  = v_x(k,1:1:i_end(k)); %saving v_x by appending it to v_x_append
                    v_z_append(k,round(start_position(k)):1:round(end_position(k)))  = v_z(k,1:1:i_end(k)); %saving v_z by appending it to v_z_append
                    density_down_total(k,1:j_end) = [density_down_total(k,1:j_end_prev),zeros(1,j_end-j_end_prev)] + density_down_temp(k,:)*time_end(k); %calculating the average density of downdard moving particles on the coarse grid
                    density_down_acc_total(k,1:j_acc_end) = density_down_acc_total(k,1:j_acc_end) + density_down_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of downdard moving particles on the fine grid
                    density_up_total(k,1:j_end) = [density_up_total(k,1:j_end_prev),zeros(1,j_end-j_end_prev)] + density_up_temp(k,:)*time_end(k); %calculating the average density of upward moving particles on the coarse grid
                    density_up_acc_total(k,1:j_acc_end) = density_up_acc_total(k,1:j_acc_end) + density_up_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of upward moving particles on the fine grid
                    bounce_length_total(k,no_bounces_total(k)+1:1:no_bounces_total(k)+no_bounces(k)) = bounce_length(k,1:no_bounces(k)); %this holds the horizontal length of individual saltating particle hops
                    vI_total(k,no_bounces_total(k)+1:no_bounces_total(k)+no_bounces(k)) = vI(k,1:no_bounces(k)); vR_total(k,no_bounces_total(k)+1:no_bounces_total(k)+no_bounces(k)) = vR(k,1:no_bounces(k)); %the average impacting and rebounding particle speed
                    thetaI_total(k,no_bounces_total(k)+1:no_bounces_total(k)+no_bounces(k)) = thetaI(k,1:no_bounces(k)); thetaR_total(k,no_bounces_total(k)+1:no_bounces_total(k)+no_bounces(k)) = thetaR(k,1:no_bounces(k)); %the average angles with the horizontal for impacting and rebounding particles
                end  %for, cycling over all particle size bins                         
                no_bounces_total = no_bounces_total + no_bounces; %keeping track of the total number of bounces or hops of saltating particles
                aerosol_emis_tot = aerosol_emis_tot + aerosol_emis; %aerosol_emis is a measure of dust aerosol emission, calculated from the total kinetic energy impacting the soil surface, and no_bounces_total holds the total number of bounces of particles in a particular particle bin
            end %if, checking whether this is the first iteration to be taken into the running average             
            start_position = end_position+1; %start_position keeps track of the starting position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
            
            %the below saves the last few minor_iteratations and carries them over to the next iteration to provide stability between iterations
            if (total_minor_iterations - minor_iteration < carry_minor_iterations) %the iteration is one of the final minor_iterations, which are carried over to the next iteration
                if(minor_iteration==total_minor_iterations - carry_minor_iterations+1) %the first minor_iteration to be carried over into the running average for the next iteration
                    start_position_carry = ones(part_no,1); %start_position keeps track of the starting position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
                    no_bounces_total_carry = no_bounces; %keeping track of the total number of bounces or hops of saltating particles
                    end_position_carry = i_end; %end_position keeps track of the ending position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
                    time_end_total_carry = time_end; %keeping track of the total time that trajectories of a particular particle bin are simulated. Used in calculating average quantities.
                    surface_coll_freq_total_carry = surface_coll_freq; %the frequency per second of surface collicions
                    v_xI_total_carry = v_xI*weighting_factor(minor_iteration); %this weighs the impacting and rebounding particle speed by the total number of saltation trajectories simulated in an iteration. The "total" variable is later divided by the sum of weight_factor to obtain the average particle speed
                    v_zI_total_carry = v_zI*weighting_factor(minor_iteration);
                    v_xR_total_carry = v_xR*weighting_factor(minor_iteration);
                    v_zR_total_carry = v_zR*weighting_factor(minor_iteration);                
                    vI_distr_total_carry = vI_distr; %the average impacting speed with which surface particle are ejected
                    v_lift_total_carry = v_lift; %the average speed with which splashed particles are ejected
                    vI_total_carry = vI; vR_total_carry = vR; %the average impacting and rebounding particle speed
                    thetaI_total_carry = thetaI; thetaR_total_carry = thetaR; %the average angles with the horizontal for impacting and rebounding particles
                    v_x_avg_total_carry = v_x_avg*time_end; %the average horizontal particle speed on the coarse grid
                    v_x_avg_acc_total_carry = v_x_avg_acc*time_end; %the average horizontal particle speed on the fine grid
                    bounce_length_total_carry = bounce_length; %this holds the horizontal length of individual saltating particle hops
                    no_rebound_prob_total_carry = no_rebound_prob.*time_end; %the average probability that an impacting particle will not rebound
                    height_total_carry = height'.*time_end; %this holds a measure of the average height of saltating particle trajectories
                    length_total_carry = length'.*time_end; %this holds a measure of the average length of saltating particle trajectories
                    alt_length_total_carry = alt_length'.*time_end; %this holds an alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed
                    for k=1:1:part_no %cycling over all particle size bins
                        z_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k))) = z(k,1:1:round(end_position_carry(k))); %saving z by appending it to z_append_carry
                        t_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k))) = t(k,1:1:round(end_position_carry(k))); %saving t by appending it to t_append_carry
                        v_x_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k))) = v_x(k,1:1:round(end_position_carry(k))); %saving v_x by appending it to v_x_append_carry
                        v_z_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k))) = v_z(k,1:1:round(end_position_carry(k))); %saving v_z by appending it to v_z_append_carry
                        density_down_total_carry(k,1:j_end) = density_down_temp(k,1:j_end)*time_end(k); density_down_acc_total_carry(k,1:j_acc_end) = density_down_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of downdard moving particles on the coarse and fine grids
                        density_up_total_carry(k,1:j_end) = density_up_temp(k,1:j_end)*time_end(k); density_up_acc_total_carry(k,1:j_acc_end) = density_up_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of upward moving particles on the coarse and fine grids
                    end %for, cycling over all particle size bins
                    aerosol_emis_tot_carry = aerosol_emis; %aerosol_emis is a measure of dust aerosol emission, calculated from the total kinetic energy impacting the soil surface, and no_bounces_total holds the total number of bounces of particles in a particular particle bin
                else %this is not the first iteration to be taken into the running average
                    end_position_carry = start_position_carry + i_end - 1; %start_position keeps track of the starting position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
                    time_end_total_carry = time_end + time_end_total_carry; %keeping track of the total time that trajectories of a particular particle bin are simulated. Used in calculating average quantities.
                    surface_coll_freq_total_carry = surface_coll_freq_total_carry + surface_coll_freq; %the frequency per second of surface collicions
                    v_xI_total_carry = v_xI_total_carry + v_xI*weighting_factor(minor_iteration); %this weighs the impacting and rebounding particle speed by the total number of saltation trajectories simulated in an iteration. The "total" variable is later divided by the sum of weight_factor to obtain the average particle speed
                    v_zI_total_carry = v_zI_total_carry + v_zI*weighting_factor(minor_iteration);
                    v_xR_total_carry = v_xR_total_carry + v_xR*weighting_factor(minor_iteration);
                    v_zR_total_carry = v_zR_total_carry + v_zR*weighting_factor(minor_iteration);
                    v_lift_total_carry(size(v_lift_total_carry,2)+1:1:size(v_lift_total_carry,2)+size(v_lift,2)) = v_lift; %the average speed with which splashed particles are ejected
                    vI_distr_total_carry(size(vI_distr_total_carry,2)+1:1:size(vI_distr_total_carry,2)+size(vI_distr,2)) = vI_distr; %the average impacting speed with which surface particle are ejected
                    v_x_avg_total_carry = v_x_avg_total_carry + v_x_avg'.*time_end; %the average horizontal particle speed on the coarse grid
                    v_x_avg_acc_total_carry = v_x_avg_acc_total_carry + v_x_avg_acc'.*time_end; %the average horizontal particle speed on the fine grid
                    no_rebound_prob_total_carry = no_rebound_prob_total_carry + no_rebound_prob.*time_end; %the average probability that an impacting particle will not rebound
                    height_total_carry = height'.*time_end + height_total_carry; %this holds a measure of the average height of saltating particle trajectories
                    length_total_carry = length'.*time_end + length_total_carry; %this holds a measure of the average length of saltating particle trajectories
                    alt_length_total_carry = alt_length'.*time_end + alt_length_total_carry; %this holds an alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed
                    for k=1:1:part_no  %cycling over all particle size bins                   
                        z_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k))) = z(k,1:1:i_end(k)); %saving z by appending it to z_append_carry
                        t_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k))) = t(k,1:1:i_end(k)) + t_append_carry(k,round(start_position_carry(k))-1); %saving t by appending it to t_append_carry
                        v_x_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k)))  = v_x(k,1:1:i_end(k)); %saving v_x by appending it to v_x_append_carry
                        v_z_append_carry(k,round(start_position_carry(k)):1:round(end_position_carry(k)))  = v_z(k,1:1:i_end(k)); %saving v_z by appending it to v_z_append_carry
                        density_down_total_carry(k,1:j_end) = [density_down_total_carry(k,1:j_end_prev),zeros(1,j_end-j_end_prev)] + density_down_temp(k,:)*time_end(k); %calculating the average density of downdard moving particles on the coarse grid 
                        density_down_acc_total_carry(k,1:j_acc_end) = density_down_acc_total_carry(k,1:j_acc_end) + density_down_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of downdard moving particles on the fine grid
                        density_up_total_carry(k,1:j_end) = [density_up_total_carry(k,1:j_end_prev),zeros(1,j_end-j_end_prev)] + density_up_temp(k,:)*time_end(k); %calculating the average density of upward moving particles on the coarse grid
                        density_up_acc_total_carry(k,1:j_acc_end) = density_up_acc_total_carry(k,1:j_acc_end) + density_up_acc_temp(k,1:j_acc_end)*time_end(k); %calculating the average density of upward moving particles on the fine grid
                        bounce_length_total_carry(k,no_bounces_total_carry(k)+1:1:no_bounces_total_carry(k)+no_bounces(k)) = bounce_length(k,1:no_bounces(k)); %this holds the horizontal length of individual saltating particle hops
                        vI_total_carry(k,no_bounces_total(k)+1:no_bounces_total(k)+no_bounces(k)) = vI(k,1:no_bounces(k)); vR_total_carry(k,no_bounces_total(k)+1:no_bounces_total(k)+no_bounces(k)) = vR(k,1:no_bounces(k)); %the average impacting and rebounding particle speed
                        thetaI_total_carry(k,no_bounces_total_carry(k)+1:no_bounces_total_carry(k)+no_bounces(k)) = thetaI(k,1:no_bounces(k)); thetaR_total_carry(k,no_bounces_total_carry(k)+1:no_bounces_total_carry(k)+no_bounces(k)) = thetaR(k,1:no_bounces(k)); %the average angles with the horizontal for impacting and rebounding particles
                    end  %for, cycling over all particle size bins                              
                    no_bounces_total_carry = no_bounces_total_carry + no_bounces; %keeping track of the total number of bounces or hops of saltating particles
                    aerosol_emis_tot_carry = aerosol_emis_tot_carry + aerosol_emis; %aerosol_emis is a measure of dust aerosol emission, calculated from the total kinetic energy impacting the soil surface, and no_bounces_total holds the total number of bounces of particles in a particular particle bin
                end %if, the first minor_iteration to be carried over into the running average for the next iteration
                start_position_carry = end_position_carry+1; %start_position keeps track of the starting position in units of model times steps in different parameters, such as z, v_x, v_z, etc.
            end %if, the iteration is one of the final minor_iterations, which are carried over to the next iteration            
        end %if, the model must equilibriate a bit first (during the initial_minor_iterations)     

        log_ejecta = log(sum(sum(created))/sum(sum(lost))) %log_ejecta = log(ejecta/#particles that did not rebound) and thus denotes the imbalance between ejected particles and particles settling back to the soil surface. The model strives for a dynamic equilibrium around log_ejecta = 0
        log_ejecta_save(cum_minor_iteration) = log_ejecta; %saving log_ejecta for debugging        
        if (minor_iteration > initial_minor_iterations) %setting n_power, which is the power to which n is raised in wind_profile to account for the deviation of log(ejecta_orbit_ratio) from 0
            n_power = (log_ejecta+log(sum(created(minor_iteration,1:part_no))/sum(lost(minor_iteration,1:part_no))))/2; %this is tuned to produce fast convergence
            u_fr_surf_save(cum_minor_iteration) = 0.4*(u_acc(2)-u_acc(1))/log(z_plot_acc(2)/z_plot_acc(1)); %%the friction velocity right at the surface (this is normally less than the friction velocity directly above the saltation layer, because of the absorption of wind momentum flux by saltating particles
        else %in this case, it's the initial wind iteration, such that n_power should be zero
            n_power = 0;
            u_fr_surf_save(cum_minor_iteration) = u_fr_thr;
        end %if, setting n_power, which is the power to which n is raised in wind_profile to account for the deviation of log(ejecta_orbit_ratio) from 0
        n_power_save(cum_minor_iteration) = n_power; %saving n_power for debugging
        
        %wind_profile calculates the wind profile, which in the absence of saltation follows the "law of the wall," but is modified by absorption of wind momentum by the saltating particles
        [n_temp, u_temp, u_acc_temp, tau_p, tau_p_acc, tau_p_surface_temp, u_fr_surf_temp(minor_iteration)] = wind_profile (d, n_relative, particle_shear_stress, particle_shear_stress_acc, total_surface_shear_stress, minor_iteration <= initial_minor_iterations, n_power, n, u_fr_surf_save(cum_minor_iteration));
        
        if (j_end > j_end_prev) %if j_end > j_end_prev, then the wind profile and other variables need to be extended to past j_end
            z_plot = delta_z:delta_z:j_end*delta_z; %extending the plotting variable
            v_x_avg_down(1:part_no,j_end_prev+1:j_end) = zeros(part_no, j_end-j_end_prev); v_x_avg_up(1:part_no,j_end_prev+1:j_end) = zeros(part_no, j_end-j_end_prev); %re-initializing the matrix that holds the average horizontal speed (because j_end might have increased)
            v_z_avg_down(1:part_no,j_end_prev+1:j_end) = zeros(part_no, j_end-j_end_prev); v_z_avg_up(1:part_no,j_end_prev+1:j_end) = zeros(part_no, j_end-j_end_prev); %re-initializing the matrix that holds the average horizontal speed (because j_end might have increased)
            density_down(1:part_no,j_end_prev+1:j_end) = zeros(part_no, j_end-j_end_prev); density_up(1:part_no,j_end_prev+1:j_end) = zeros(part_no, j_end-j_end_prev); %initializing the concentration matrix, describing the likelihood that the particle is at a particular height.
            u(j_end_prev+1:j_end) = u(j_end_prev) + (u_fr/0.4)*(log(delta_z*(j_end_prev+1:j_end)/(j_end_prev*delta_z)));
        end %if, if j_end > j_end_prev, then the wind profile needs and other variables need to be extended to past j_end
        tau_p_total(minor_iteration, 1:size(u_temp,2)) = tau_p; %storing the particle shear stress on the coarse grid (tau_p), for later calculation of the average particle shear stress
        tau_p_acc_total(minor_iteration, 1:size(tau_p_acc,2)) = tau_p_acc; %storing the particle shear stress on the fine grid (tau_p_acc), for later calculation of the average particle shear stress on the fine grid
        tau_p_surface_total(cum_minor_iteration) = tau_p_surface_temp; %the particle shear at the surface, as calculated by the individual minor iterations
        tau_p_surface = mean(tau_p_acc_total(:,1)); %the particle shear stress at the surface
        n_save(cum_minor_iteration,1:part_no) = n; %saving the particle concentration for debugging
        if (major_iteration==1) %this calculates the particle concentration and wind profile, which differs depending on whether this is the initial major iteration
            n = ((minor_iteration-1)*n+n_temp)/minor_iteration; %calculating the particle concentration
            u = ((minor_iteration)*u+u_temp)/(minor_iteration+1); %calculating the wind profile on the coarse grid
            u_acc = ((minor_iteration)*u_acc+u_acc_temp)/(minor_iteration+1); %calculating the wind profile on the fine grid
        else %in this case it's not the initial major_iteration
            n = ((total_minor_iterations-1)*n+n_temp)/total_minor_iterations; %calculating the particle concentration
            u = ((total_minor_iterations-1)*u+u_temp)/total_minor_iterations; %calculating the wind profile on the coarse grid
            u_acc = ((total_minor_iterations-1)*u_acc+u_acc_temp)/total_minor_iterations; %calculating the wind profile on the fine grid
        end %if, this calculates the particle concentration and wind profile, which differs depending on whether this is the initial major iteration
    end %while, performing a number of iterations while continuously adjusting the wind profile, particle trajectories, and particle concentration
    
    if (plot_wind_profile) %plotting the wind profile
        if (major_iteration==1)
            figure(1); clf; %clearing the figure in case this is the first iteration
        end %if
        figure(1); hold on; plot([0,u_acc],[z0,z_plot_acc],'r'); plot(u,z_plot(1:j_end),'g'); %plotting the wind profile on both the fine (in red) and coarse (in green) grids
    end %if, plotting the wind profile

    %the following calculates various parameters from the iterative loop above    
    j_end_save (major_iteration) = j_end; %saving j_end, which will differ among iterations because j_end is increased when particles trajectories go above delta_z*j_end
    z = z_append; t = t_append; v_x = v_x_append; v_z = v_z_append; %obtaining z, t, v_x, and v_z
    time_end = time_end_total; %keeping track of the total time that trajectories of a particular particle bin are simulated. Used in calculating average quantities.
    salt_mass_fraction = n.*m/sum(n.*m); %the fraction of the total mass of saltating particles contained in each particle bin
    surface_coll_freq = surface_coll_freq_total./time_end; %the frequency per second of surface collicions
    surface_coll_freq_time(major_iteration) = salt_mass_fraction*surface_coll_freq; %saving the frequency per second of surface collicions for output
    if (size(pEj,2)>total_minor_iterations) %this stores the particle ejection matrix pEj, but limits its size for computational efficiency by randomly selecting a subset
        pEj_temp = zeros(part_no,1);
        for k=1:part_no
            for p=1:1:min(2*total_minor_iterations,ejecta_totals(k))
                a=rand(1);
                pEj_temp(k,p) = pEj(k,round(a*(ejecta_totals(k)-1)+1));
            end
        end
        pEj = pEj_temp;
        ejecta_totals = min(2*total_minor_iterations,ejecta_totals);
    end %if, this stores the particle ejection matrix pEj, but limits its size for computational efficiency by randomly selecting a subset
    salt_mass_fraction = n.*m/sum(n.*m); %the fraction of the total mass of saltating particles contained in each particle bin
    avg_vI_distr(major_iteration) = sum(vI_distr_total)/size(vI_distr_total,2); %the average impacting speed with which surface particle are ejected
    avg_v_lift(major_iteration) = sum(v_lift_total)/size(v_lift_total,2); %saving the average speed with which splashed particles are ejected
    v_xI_time(major_iteration) = (v_xI_total/sum(weighting_factor))'*salt_mass_fraction'; %calculating the average impacting and rebounding particle speeds by dividing by the sum of the weight functions, and scaling by the mass fraction of the PSD
    v_zI_time(major_iteration) = (v_zI_total/sum(weighting_factor))'*salt_mass_fraction';
    v_xR_time(major_iteration) = (v_xR_total/sum(weighting_factor))'*salt_mass_fraction';   
    v_zR_time(major_iteration) = (v_zR_total/sum(weighting_factor))'*salt_mass_fraction';
    delVx_time(major_iteration) = v_xI_time(major_iteration) - v_xR_time(major_iteration); %saving the average particle speed gained from acceleration by the wind between ejection and impact
    v_x_avg_time(major_iteration) = (v_x_avg_total./time_end)'*salt_mass_fraction'; %saving the average horizontal particle speed on the coarse grid
    v_x_avg_acc_time(major_iteration) = (v_x_avg_acc_total./time_end)'*salt_mass_fraction'; %saving the average horizontal particle speed on the fine grid
    v_x_avg_down = sum(v_x_avg_down_total,3)./(sum(counter_down_total,3)+10^(-50)); %saving the average horizontal and vertical particle speeds of upward and downward moving particles on the fine and coarse grids. The very small number is in the denominator to prevent division by zero - it does not affect the results
    v_x_avg_down_acc = sum(v_x_avg_down_acc_total,3)./(sum(counter_down_acc_total,3)+10^(-50)); 
    v_z_avg_down = sum(v_z_avg_down_total,3)./(sum(counter_down_total,3)+10^(-50));
    v_z_avg_down_acc = sum(v_z_avg_down_acc_total,3)./(sum(counter_down_acc_total,3)+10^(-50)); 
    v_x_avg_up = sum(v_x_avg_up_total,3)./(sum(counter_up_total,3)+10^(-50));
    v_x_avg_up_acc = sum(v_x_avg_up_acc_total,3)./(sum(counter_up_acc_total,3)+10^(-50));
    v_z_avg_up_acc = sum(v_z_avg_up_acc_total,3)./(sum(counter_up_acc_total,3)+10^(-50));    
    v_z_avg_up = sum(v_z_avg_up_total,3)./(sum(counter_up_total,3)+10^(-50));
    log_ejecta_orbit_ratio(major_iteration) = log_ejecta %saving the average log_ejecta_orbit_ratio. This is an important indicator of steady-state saltation and is used to assess model convergence below
    vertflux_ratio_save(major_iteration) = sum(sum(rebounds_ejecta))/sum(sum(impacts)) %saving the ratio of particles impacting the soil bed to outgoing particles (rebounding + ejected) and should be ~1 in steady state.
    for k=1:1:part_no
        density_down(k,1:j_end) = density_down_total(k,1:j_end)/time_end(k); %reconstructing the average density of downward-moving particles on the coarse grid
        density_down_acc(k,1:j_acc_end) = density_down_acc_total(k,1:j_acc_end)/time_end(k); %reconstructing the average density of downward-moving particles on the fine grid
        density_up(k,1:j_end) = density_up_total(k,1:j_end)/time_end(k); %reconstructing the average density of upward-moving particles on the coarse grid
        density_up_acc(k,1:j_acc_end) = density_up_acc_total(k,1:j_acc_end)/time_end(k); %reconstructing the average density of upward-moving particles on the fine grid
        v_x_avg(k,major_iteration) = sum(((v_x(k,2:end_position(k))+v_x(k,1:end_position(k)-1))/2).*(t(k,2:end_position(k))-t(k,1:end_position(k)-1)))/time_end_total(k); %calculating the average particle speed
        vI_time(major_iteration,k) = mean(vI_total(k,1:no_bounces_total(k))); vR_time(major_iteration,k) = mean(vR_total(k,1:no_bounces_total(k))); %saving the average impacting and rebounding particle speed
        vR_time(major_iteration,k) = mean(vR_total(k,1:no_bounces_total(k))); %saving the average impacting and rebounding particle speed        
        thetaI_time(major_iteration,k) = mean(thetaI_total(k,1:no_bounces_total(k))); thetaR_time(major_iteration,k) = mean(thetaR_total(k,1:no_bounces_total(k))); %saving the average angles with the horizontal for impacting and rebounding particles
    end            
    density_down_time(major_iteration,1:part_no,1:j_end) = density_down; density_up_time(major_iteration,1:part_no,1:j_end) = density_up; %saving the particle density on the coarse grid
    density_down_acc_time(major_iteration,1:part_no,1:j_acc_end) = density_down_acc; density_up_acc_time(major_iteration,1:part_no,1:j_acc_end) = density_up_acc; %saving the particle density on the coarse grid
    v_x_avg_time(major_iteration) = v_x_avg(:,major_iteration)'*salt_mass_fraction'; %calculating the average particle speed, averaged over all particle bins
    no_rebound_prob = no_rebound_prob_total./time_end; %the average probability that an impacting particle will not rebound
    no_rebound_prob_time(major_iteration) = salt_mass_fraction*no_rebound_prob; %saving the average probability that an impacting particle will not rebound for output
    height = height_total./time_end; %this holds a measure of the average height of saltating particle trajectories
    length = length_total./time_end; %this holds a measure of the average length of saltating particle trajectories
    alt_length = alt_length_total./time_end; %this holds an alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed
    z0s_time(major_iteration) = z_plot(j_end)*exp(-0.4*u(j_end)/u_fr); %the average aerodynamic roughness length of the saltation layer. Can be compared to measurements
    tau_p = sum(tau_p_total,1)/size(tau_p_total,1); %the average particle shear stress on the coarse grid
    tau_p_acc = sum(tau_p_acc_total,1)/size(tau_p_acc_total,1); %the average particle shear stress on the fine grid
    tau_p_time(major_iteration,1:1:size(tau_p,2)) = tau_p(1:1:size(tau_p,2)); %saving the average particle shear stress on the coarse grid
    tau_p_acc_time(major_iteration,1:1:size(tau_p_acc,2)) = tau_p_acc(1:1:size(tau_p_acc,2)); %saving the average particle shear stress on the fine grid
    wind_time(major_iteration,1:1:size(u,2)) = u(1:1:size(u,2)); %saving the wind profile on the coarse grid
    wind_acc_time(major_iteration,1:1:size(u_acc,2)) = u_acc(1:1:size(u_acc,2)); %saving the wind profile on the fine grid
    u_fr_surf(major_iteration) = 0.4*(u_acc(2)-u_acc(1))/log(z_plot_acc(2)/z_plot_acc(1))  %saving the friction velocity right at the surface (this is normally less than the friction velocity directly above the saltation layer, because of the absorption of wind momentum flux by saltating particles
    total_mass(major_iteration) = (pi/6)*n*rho_p*d'.^3; %saving the total mass of saltating particles per square meter of the soil surface
    height_time(major_iteration) = n*height/sum(n); %saving this measure of the average height of saltating particle trajectories
    length_time(major_iteration) = n*length/sum(n); %saving this measure of the average length of saltating particle trajectories
    alt_length_time(major_iteration) = n*alt_length/sum(n); %saving this alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed
    n_time(major_iteration,1:part_no) = n; %saving the particle concentration (in # per square meter)
    tau_p_surface_time(major_iteration) = tau_p_surface; %saving the particle shear stress at the surface
    aerosol_emis_time(major_iteration) = aerosol_emis_tot/(total_minor_iterations-initial_minor_iterations); %saving aerosol_emis, which is a measure of dust aerosol emission, calculated from the total kinetic energy impacting the soil surface, and no_bounces_total holds the total number of bounces of particles in a particular particle bin    
    
    %the following calculates the horizontal and vertical mass flux
    [z50(major_iteration), z50mass(major_iteration), z50_bin_temp, L50(major_iteration), total_mass_flux(major_iteration), vert_mf_temp, vert_mf_acc_temp, hor_mf_temp, z_s(major_iteration)] = calc_mass_flux(bounce_length_total, no_bounces_total, hor_interval, time_end, density_down, density_up, v_x_avg_down, v_x_avg_up, density_down_acc, density_up_acc, v_x_avg_down_acc, v_x_avg_up_acc, n, d, max(max(z)));
    vert_mf(major_iteration,1:part_no,1:j_end) = vert_mf_temp(1:part_no,1:j_end); %saving the vertical particle mass flux on the coarse grid
    vert_mf_acc(major_iteration,1:part_no,1:j_acc_end) = vert_mf_acc_temp(1:part_no,1:j_acc_end); %saving the vertical particle mass flux on the fine grid
    hor_mf(major_iteration,1:size(hor_mf_temp,1),1:size(hor_mf_temp,2)) = hor_mf_temp(:,:); %saving the horizontal particle mass flux
    z50_bin(major_iteration,1:part_no) = z50_bin_temp; %z50_bin is the height below which half of the mass flux of a *a particular particle bin* occurs. z50 is the height below which half of the total particle mass flux occurs
    z50, z50_bin, total_mass_flux
    
    up_fraction(1:j_end) = salt_mass_fraction*min(1,density_up(:,1:j_end)./max(10^(-10),density_up(:,1:j_end)+density_down(:,1:j_end))); up_fraction_acc(1:j_acc_end)= up_fraction(no_fine_grid) - 2*((j_acc_end-(1:j_acc_end))/j_acc_end)*(up_fraction(no_fine_grid) - sum(density_up_acc)/sum(density_up_acc+density_down_acc)); %calculating the fraction of particles that are moving upwards. This is generally less than 50 %, because upward-moving particles have higher speed than downward-moving particles because of fluid drag
    for j=1:1:j_end %this calculates the average horizontal particle speed as a function of height on the coarse grid
        v_x_total(major_iteration,j) = salt_mass_fraction*(v_x_avg_down(:,j)*(1-up_fraction(j))+v_x_avg_up(:,j)*up_fraction(j));
    end %for, this calculates the average horizontal particle speed as a function of height on the coarse grid
    for j_acc=1:1:j_acc_end %this calculates the average horizontal particle speed as a function of height on the fine grid
        v_x_acc_total(major_iteration,j_acc) = salt_mass_fraction*(v_x_avg_down_acc(:,j_acc)*(1-up_fraction_acc(j_acc))+v_x_avg_up_acc(:,j_acc)*up_fraction_acc(j_acc));
    end %for, this calculates the average horizontal particle speed as a function of height on the coarse grid
    del = floor(min(d/2)/delta_z); del_acc = floor(min(d)/delta_z_acc); %del is used in the next line to correct v_x_total when min(d) is larger than a grid box
    v_x_acc_total(major_iteration,1:del_acc) = v_x_acc_total(major_iteration,1+del_acc); v_x_total(major_iteration,1:del) = v_x_total(major_iteration,1+del); %correcting v_x_total in the grid boxes that are below min(d), where changes in particle speed between grid boxes is meaningless
  
    if (rem(major_iteration,total_major_iterations)==0 && major_iteration > 2*total_major_iterations) %to prevent the model from looping indefinitely, the tolerance is doubled once major_iteration == (3,4,5,...)*total_major_iterations
        tolerance = 2*tolerance;
    end
    if (major_iteration > total_major_iterations+1) %checking if the model converged
        convergence = true;
        z50_avg = sum(z50(major_iteration+1-total_major_iterations:major_iteration))/total_major_iterations %the mean of z50 of the past total_major_iterations major iterations
        mass_flux_avg = sum(total_mass_flux(major_iteration+1-total_major_iterations:major_iteration)+0.001)/total_major_iterations; %the mean of the mass flux of the past total_major_iterations major iterations. The fraction 0.001 is added to have it converge at low wind speeds, where the absence of Owen's hypothesis makes the mass flux very sensitively dependent on u_fr_surf
        fit_z50 = polyfit(-(total_major_iterations-1)/2:(total_major_iterations-1)/2,z50(major_iteration+1-total_major_iterations:major_iteration),1); %fitting a straight line to z50
        fit_MF = polyfit(-(total_major_iterations-1)/2:(total_major_iterations-1)/2,total_mass_flux(major_iteration+1-total_major_iterations:major_iteration)+0.001,1); %fitting a straight line to the mass flux. The fraction 0.001 is added to have it converge at low wind speeds, where the absence of Owen's hypothesis makes the mass flux very sensitively dependent on u_fr_surf
        %Below are a number of criteria that need to be satisfied for the model results to represent convergence
        if ((Owen==false && (splash_fraction == 1 || min(u_fr_surf(major_iteration+1-total_major_iterations:major_iteration)-u_fr_thr_soil) < 0) && abs(sum(log_ejecta_orbit_ratio(major_iteration+1-total_major_iterations:major_iteration))/total_major_iterations) > tolerance/2)) %checking whether saltation is in steady-state
            convergence = false;
            'saltation not in steady state'
        end %if, checking whether saltation is in steady-state
        if (std(z50(major_iteration+1-total_major_iterations:major_iteration))/z50_avg > 0.5*tolerance) %checking whether the z50 standard deviation is acceptable
            convergence = false;
            'z50 standard deviation exceeds tolerance'
        end %if checking whether the z50 standard deviation is acceptable
        if (abs(fit_z50(1)/fit_z50(2)) > tolerance/total_major_iterations) %checking whether the trend in z50 is acceptable
            convergence = false;
            'trend in z50 exceeds tolerance'
        end %if, checking whether the trend in z50 is acceptable
        if (std(total_mass_flux(major_iteration+1-total_major_iterations:major_iteration))/mass_flux_avg > tolerance) %checking whether the mass flux standard deviation is acceptable
            convergence = false;
            'mass flux standard deviation exceeds tolerance'
        end %if, checking whether the mass flux standard deviation is acceptable
        if (abs(total_major_iterations*fit_MF(1))/fit_MF(2) > tolerance) %checking whether the trend in the mass flux is acceptable
            convergence = false;
            'trend in mass flux exceeds tolerance'
        end %if, checking whether the trend in the mass flux is acceptable
        if (convergence == false && major_iteration > 10*total_major_iterations) %to prevent from looping indefinitely, the model will terminate if convergence has not been achieved after 10*total_major_iterations major iterations
            error('model is not converging');
        end %if, to prevent from looping indefinitely, the model will terminate if convergence has not been achieved after 10*total_major_iterations major iterations
    end %if, checking if the model converged   
end %for, iterating until convergence is achieved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Declaration of files to which the data will be written:   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fidWind = fopen(strcat('Wind',filetext,'.txt'), 'wt'); %file holding the wind profile
fidVertMassFlux = fopen(strcat('VMF',filetext,'.txt'), 'wt'); %file holding the vertical mass flux profile
fidHorMassFlux = fopen(strcat('HMF',filetext,'.txt'), 'wt'); %file holding the horizontal mass flux profile
fidParameters = fopen(strcat('Par',filetext,'.txt'),'wt'); %file holding important parameters
fidImpReb = fopen(strcat('ImpReb',filetext,'.txt'),'wt'); %file holding impact and rebound angles and speeds
fidPSD = fopen(strcat('PSD',filetext,'.txt'), 'wt'); %file holding the size distribution of saltating particles
fidVx = fopen(strcat('Vx',filetext,'.txt'), 'wt'); %file holding the horizontal particle speed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%         The following writes data to files         %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculating average parameters that will be written to files
total_iterations = major_iteration; %the total number of iterations before convergence was achieved
avg_z50 = mean(z50(total_iterations-total_major_iterations+1:1:total_iterations)); %the height z50 below which half of the mass flux takes place, averaged over the last several major iterations
avg_z50_mass = mean(z50mass(total_iterations-total_major_iterations+1:1:total_iterations)); %the height z50_mass below which half of the mass of saltating particles is contained, averaged over the last several major iterations
z50_bin_avg = sum(z50_bin(total_iterations-total_major_iterations+1:1:total_iterations,:))/total_major_iterations; %the height z50 for the individual particle bins, averaged over the last several major iterations
avg_vI = mean(vI_time(total_iterations-total_major_iterations+1:1:total_iterations)); %the impacting particle speed, averaged over the last several major iterations
avg_thetaI = mean(thetaI_time(total_iterations-total_major_iterations+1:1:total_iterations)); %the angle with the horizontal for impacting particles, averaged over the last several major iterations
avg_vx = mean(v_x_avg_time(total_iterations-total_major_iterations+1:1:total_iterations)); %the mean horizontal particle speed, averaged over the last several major iterations
avg_aer_emis = mean(aerosol_emis_time(total_iterations-total_major_iterations+1:1:total_iterations)); %aerosol_emis, which is a measure of dust aerosol emission, averaged over the last several major iterations
avg_z0s = mean(z0s_time(total_iterations-total_major_iterations+1:1:total_iterations)); %the mean aerodynamic roughness length, averaged over the last several major iterations
avg_log_ejecta = mean(log_ejecta_orbit_ratio(total_iterations-total_major_iterations+1:1:total_iterations)); %the mean log_ejecta, which must be around zero for the model to be in steady-state, averaged over the last several major iterations
avg_mass = mean(total_mass(total_iterations-total_major_iterations+1:1:total_iterations)); %the load of saltating particles in kilograms per square meter, averaged over the last several major iterations
avg_massflux = mean(total_mass_flux(total_iterations-total_major_iterations+1:1:total_iterations)); %the particle mass flux, averaged over the last several major iterations
avg_vertflux_ratio = mean(vertflux_ratio_save(total_iterations-total_major_iterations+1:1:total_iterations)); %the ratio of impacting to ejected+rebounding particles, averaged over the last several major iterations
avg_length = mean(length_time(total_iterations-total_major_iterations+1:1:total_iterations)); %a measure of the average length of saltating trajectories, averaged over the last several major iterations
avg_alt_length = mean(alt_length_time(total_iterations-total_major_iterations+1:1:total_iterations)); %an alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed, averaged over the last several major iterations
avg_delVx = mean(delVx_time(total_iterations-total_major_iterations+1:1:total_iterations)); %the average particle speed gained from acceleration by the wind between ejection and impact, averaged over the last several major iterations

%The below calculates the average vertical profiles of several quantities
for j=1:j_end
    j_nonzero(j) = size(find(v_x_total(total_iterations-total_major_iterations+1:1:total_iterations,j)>0),1);
    if (j_nonzero(j)==0) %in this case, there's no data above j
        break;
    end
    v_x_total_avg(j) = sum(v_x_total(total_iterations-total_major_iterations+1:1:total_iterations,j))/j_nonzero(j); %the average horizontal particle speed on the coarse grid    
    avg_density(j) = sum(density_down_time(total_iterations-total_major_iterations+1:1:total_iterations,1,j))/(j_nonzero(j)*delta_z) + sum(density_up_time(total_iterations-total_major_iterations+1:1:total_iterations,1,j))/(j_nonzero(j)*delta_z);
    for k=1:1:part_no
        vert_mf_total(k,j) = sum(vert_mf(total_iterations-total_major_iterations+1:1:total_iterations,k,j))/j_nonzero(j); %the vertical mass flux per particle bin on the coarse grid, averaged over the last total_major_iterations major iterations
    end    
end
for j_acc=1:j_acc_end
    j_acc_nonzero(j_acc) = size(find(v_x_acc_total(total_iterations-total_major_iterations+1:1:total_iterations,j_acc)>0),1);
    if (j_acc_nonzero(j_acc)==0) %in this case, there's no data above j
        break;
    end
    v_x_acc_total_avg(j_acc) = sum(v_x_acc_total(total_iterations-total_major_iterations+1:1:total_iterations,j_acc))/j_acc_nonzero(j_acc); %the average horizontal particle speed on the fine grid    
    avg_density_acc(j_acc) = sum(density_down_acc_time(total_iterations-total_major_iterations+1:1:total_iterations,1,j_acc))/(j_acc_nonzero(j_acc)*delta_z_acc) + sum(density_up_acc_time(total_iterations-total_major_iterations+1:1:total_iterations,1,j_acc))/(j_acc_nonzero(j_acc)*delta_z_acc);
    for k=1:1:part_no
        vert_mf_acc_total(k,j_acc) = sum(vert_mf_acc(total_iterations-total_major_iterations+1:1:total_iterations,k,j_acc))/j_acc_nonzero(j_acc); %the vertical mass flux per particle bin on the fine grid, averaged over the last total_major_iterations major iterations        
    end
end

%saving the wind profile
for p=total_iterations-total_major_iterations+1:1:total_iterations %this accounts for differences in j_end for the last total_major_iterations major iterations, which affect the wind profile saved for each iteration
    wind_time(p,j_end_save(p)+1:j_end) = (wind_time(p,j_end_save(p)+1:j_end) + (u_fr/0.4)*(log(delta_z*(j_end_save(p)+1:j_end)/(j_end_save(p)*delta_z))));
    if (j_end_save(p)>j_end_save(p-1))
        for s=1:p-1
            wind_time(s,j_end_save(p-1)+1:j_end_save(p)) = wind_time(p,j_end_save(p-1)+1:j_end_save(p));
            wind_time(s,j_end_save(p-1)+1:j_end_save(p)) = wind_time(s,j_end_save(p-1)) + (u_fr/0.4)*(log(delta_z*(j_end_save(p-1)+1:j_end_save(p))/(j_end_save(p-1)*delta_z)));
        end
    end
end %for, this accounts for differences in j_end for the last total_major_iterations major iterations, which affect the wind profile saved for each major_iteration
wind_time_avg = sum(wind_time(total_iterations-total_major_iterations+1:1:total_iterations,:))/total_major_iterations; %the wind profile on the coarse grid averaged over the last total_major_iterations major iterations
wind_acc_time_avg = sum(wind_acc_time(total_iterations-total_major_iterations+1:1:total_iterations,:))/total_major_iterations; %the wind profile on the fine grid averaged over the last total_major_iterations major iterations
dUdz(1) = (-wind_acc_time_avg(5)+8*wind_acc_time_avg(4)-8*wind_acc_time_avg(2)+wind_acc_time_avg(1))/(12*delta_z_acc); dUdz(2) = dUdz(1); dUdz(3) = dUdz(1); %calculating the dUdz to obtain the surface friction speed and the grain-borne shear stress. The latter is assumed constant in the first 3 layers on the fine grid in order to obtain an accurate numerical differentiation of the wind profile.
avg_ufr_surf = 0.4*z_plot_acc(3)*dUdz(3);  %the mean surface friction speed
fprintf(fidWind, 'j   z   Wind     tau_p \n');
fprintf(fidWind, '%3.0f  %1.3e  %3.5f  %1.5f  \n', 0, z0, 0, rho*(u_fr^2-avg_ufr_surf^2));
for j_acc=1:max(3,floor(d/(2*delta_z_acc)+1))
    fprintf(fidWind, '%3.0f  %3.5f  %3.5f  %1.5f \n', j_acc, z_plot_acc(j_acc), wind_acc_time_avg(j_acc), rho*(u_fr^2-avg_ufr_surf^2));
end
for j_acc=(max(3,floor(d/(2*delta_z_acc)+1))+1):1:size(u_acc,2)-1
    if (j_acc<size(u_acc,2)-2)
        dUdz(j_acc) = (-wind_acc_time_avg(j_acc+2)+8*wind_acc_time_avg(j_acc+1)-8*wind_acc_time_avg(j_acc-1)+wind_acc_time_avg(j_acc-2))/(12*delta_z_acc);
    else
        dUdz(j_acc) = (wind_acc_time_avg(j_acc+1)-wind_acc_time_avg(j_acc))/(delta_z_acc);
    end
    tau_p_acc_time_avg(j_acc) = rho*u_fr^2-rho*(0.4*z_plot_acc(j_acc)*dUdz(j_acc))^2;
    fprintf(fidWind, '%3.0f  %3.5f  %3.5f  %1.5f \n', j_acc, z_plot_acc(j_acc), wind_acc_time_avg(j_acc), tau_p_acc_time_avg(j_acc));
end
for j=no_fine_grid:1:j_end-2
    dUdz(j) = (-wind_time_avg(j+2)+8*wind_time_avg(j+1)-8*wind_time_avg(j-1)+wind_time_avg(j-2))/(12*delta_z);
    tau_p_time_avg(j) = rho*u_fr^2-rho*(0.4*z_plot(j)*dUdz(j))^2;
    fprintf(fidWind, '%3.0f  %3.5f  %3.5f  %3.5f \n', j, z_plot(j), wind_time_avg(j), tau_p_time_avg(j));
end

%saving the vertical and horizontal mass flux
for k=1:1:part_no
    hor_mf_total(k,:) = sum(hor_mf(total_iterations-total_major_iterations+1:1:total_iterations,k,:))/total_major_iterations; %the horizontal mass flux per particle bin, averaged over the last total_major_iterations major iterations
end
vert_mf_total(part_no+1,:) = sum(vert_mf_total,1); %the column part_no+1 holds the total vertical mass flux on the coarse grid, summed over all particle bins
vert_mf_total_norm = vert_mf_total(part_no+1,:)/mass_flux_avg; %the vertical mass flux on the coarse grid, normalized by the total mass flux
vert_mf_acc_total(part_no+1,:) = sum(vert_mf_acc_total,1); %the column part_no+1 holds the total vertical mass flux on the coarse grid, summed over all particle bins
vert_mf_acc_total_norm = vert_mf_acc_total(part_no+1,:)/mass_flux_avg; %the vertical mass flux on the fine grid, normalized by the total mass flux
hor_mf_total(part_no+1,:) = sum(hor_mf_total,1); %the column part_no+1 holds the total horizontal mass flux, summed over all particle bins
hor_mf_total_norm = hor_mf_total(part_no+1,:)/mass_flux_avg; %the horizontal mass flux, normalized by the total mass flux
hor_mf_total_norm = (1-hor_mf_total_norm(1)/sum(hor_mf_total_norm))*hor_mf_total_norm; %correcting the normalized horizontal mass flux for the large contribution of the lowest bin, where particles are in creep rather than saltation
hor_mf_dist(1,1:size(hor_mf,3)) = 0.5*hor_interval:hor_interval:(size(hor_mf,3)-0.5)*hor_interval; %the horizontal distance coordinate corresponding to hor_mf_total and hor_mf_total_norm
hmf001z = 0;
for j=1:j_end 
    if (hor_mf_total_norm(j)<0.01) %finding hmf001z
        hmf001z = hor_mf_dist(j); %the horizontal distance where the normalized horizontal mass flux is smaller than 0.01. This gives an indication of the length of saltation trajectories and the decay of the horizontal saltation mass flux with distance
        break;
    end %if, finding hmf001z
end %for
fprintf(fidVertMassFlux, 'z        VertMFtot VermMFnorm ');
for k=1:1:part_no %writing the headers of VMF.txt
    fprintf(fidVertMassFlux, 'd%3.0f', d(k)*10^6);
    fprintf(fidVertMassFlux, 'um      ');
end %for
fprintf(fidVertMassFlux, '\n');
for j=1:size(vert_mf_total,2) %writing the vertical mass flux profile to VMF.txt
    fprintf(fidVertMassFlux, '%3.5f  %2.8f  %4.6f', (j-0.5)*delta_z, vert_mf_total(part_no+1,j), vert_mf_total_norm(j));
    for k=1:1:part_no
        fprintf(fidVertMassFlux, '  %2.8f', vert_mf_total(k,j));
    end %for
    fprintf(fidVertMassFlux, '\n');
end %for
fprintf(fidHorMassFlux, 'z        HorMFtot HorMFnorm  ');
for k=1:1:part_no %writing the headers of HMF.txt
    fprintf(fidHorMassFlux, 'd%3.0f', d(k)*10^6); 
    fprintf(fidHorMassFlux, 'um      ');
end %for
fprintf(fidHorMassFlux, '\n');
for j=1:size(hor_mf_dist,2) %writing the vertical mass flux profile to VMF.txt
    fprintf(fidHorMassFlux, '%3.5f  %2.8f  %4.6f', hor_mf_dist(j), hor_mf_total(part_no+1,j), hor_mf_total_norm(j));
    for k=1:1:part_no
        fprintf(fidHorMassFlux, '  %2.8f', hor_mf_total(k,j));
    end
    fprintf(fidHorMassFlux, '\n');
end

%saving important parameters to Par.txt
fprintf(fidParameters, 'ufr=%2.3f  ufrthr=%2.3f  h=%1.4f z0=%1.6f  delta_z=%1.4f t_max=%2.2f  a_splash=%1.4f d=%1.6f \n', u_fr, u_fr_thr, h, z0, delta_z, t_max, a_splash, psd_flag);
fprintf(fidParameters, 'elevation=%4.0f MinorIter=%2.0f S%1.1f hmf001z=%2.4f\n', elevation, total_minor_iterations, S, hmf001z);
fprintf(fidParameters, 'Iteration z50     vI     Vx     z0s        ufr_surf log_ejecta  mass     massflux length z50_mass \n');
for k=1:1:total_iterations
    fprintf(fidParameters, '%2.0f        %1.5f %2.4f %2.4f %1.3e %1.4f   %+1.3f      %1.6f %1.6f %2.4f %1.4f   \n', k, z50(k), vI_time(k), v_x_avg_time(k), z0s_time(k), u_fr_surf(k), log_ejecta_orbit_ratio(k), total_mass(k), total_mass_flux(k), length_time(k), z50mass(k));
end
fprintf(fidParameters, 'Avg:      %1.5f %2.4f %2.4f %1.3e %1.4f   %+1.3f      %1.6f %1.6f %2.4f %1.4f   \n', avg_z50, avg_vI, avg_vx, avg_z0s, avg_ufr_surf, avg_log_ejecta, avg_mass, avg_massflux, avg_length, avg_z50_mass);

%saving impact and rebound speeds and angles to ImpReb.txt
fprintf(fidImpReb, 'Iteration vI    thetaI     vR    thetaR  \n');
for k=1:1:total_iterations
    fprintf(fidImpReb, '%2.0f        %2.4f  %2.4f  %2.4f  %2.4f \n', k, vI_time(k), thetaI_time(k), vR_time(k), thetaR_time(k));
end
fprintf(fidImpReb, 'Avg:      %2.4f  %2.4f  %2.4f  %2.4f \n', avg_vI, mean(thetaI_time(total_iterations-total_major_iterations+1:1:total_iterations)), mean(vR_time(total_iterations-total_major_iterations+1:1:total_iterations)), mean(thetaR_time(total_iterations-total_major_iterations+1:1:total_iterations)));

%saving the size distribution of saltating particles to PSD.txt
fprintf(fidPSD, 'D        MassFract \n');
for k=1:part_no
    fprintf(fidPSD, '%1.6f  %2.6f \n', d(k), (sum(vert_mf_acc_total(k,:)*delta_z_acc)+sum((vert_mf_total(k,no_fine_grid+1:size(vert_mf_total,2)))*delta_z))/(sum(vert_mf_acc_total(part_no+1,:)*delta_z_acc)+sum((vert_mf_total(part_no+1,no_fine_grid+1:size(vert_mf_total,2)))*delta_z)));
end

%saving the horizontal particle speed to Vx.txt
fprintf(fidVx, 'z         Vx \n');
for j=1:min(j_acc_end,size(v_x_acc_total_avg,2))
    fprintf(fidVx,'%2.6f  %2.6f\n', (j-1)*delta_z_acc, v_x_acc_total_avg(j));
end
for j=no_fine_grid+1:size(v_x_total_avg,2)
    fprintf(fidVx,'%2.4f  %2.4f\n', (j-1)*delta_z, v_x_total_avg(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%         The following closes the files to which data was written         %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fidHorMassFlux);
fclose(fidVertMassFlux);
fclose(fidParameters);
fclose(fidImpReb);
fclose(fidWind);
fclose(fidPSD);
fclose(fidVx);