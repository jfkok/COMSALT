function [t, x, z, v_x, v_z, t_end, i_end, surface_coll_freq, no_rebound_prob, height, length, alt_length, v_x_impact, v_z_impact, vI, thetaI, v_x_rebound, v_z_rebound, vR, thetaR, vI_distr, surface_shear_stress, bounce_length, no_bounces, n_relative, created, lost, impacts, rebounds_ejecta, ejecta_orbit_ratio, weighting_factor_avg, weighting_factor_part, pEj, ejecta_totals, v_lift, aerosol_emis_tot, j_end_prev] = saltation_trajectory (d, u, u_acc, pEj, n, mass_fraction, particle_coll_freq, v_x_avg_down, v_x_avg_up, v_z_avg_down, v_z_avg_up, density_down, density_up, mean_weighting_factor_part, u_fr_surf, z_s)

%saltation_trajectory calculates the trajectories of saltating particles,
%and hands all relevant parameters back to the calling function, saltation_model

%these are the global parameters used by saltation_trajectory
global h t_max rho mu gamma sigma_W_prop sigma_U_prop u_fr no_fine_grid j_end j_acc_end;
global delta_z delta_z_acc P_reb alpha_ej_fract turbulence g d_mean;
global rho_p pi part_no beta_time spin_sigma spin_mu a_splash beta min_lift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   This portion of the code initializes various variables used in saltation_trajectory   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=d/2;  %particle radius in m
m = (4/3)*pi*rho_p*r.^3; %particle mass in kg
tau_rot = rho_p*d.^2/(60*mu); %the rotational relaxation time
h_permanent = h;  %the maximum time step, used when no collisions occurred in the particle's recent history (i.e., 4 timesteps)
h_bounce = h/16;  %after a collision, there is more numerical error in calculating the particle's trajectory, and so the time step is decreased in an adaptive manner to preserve model accuracy
x=zeros(part_no,1); %the horizontal coordinate
z=zeros(part_no,1); %the vertical coordinate
v_x=zeros(part_no,1);  %the horizontal speed
v_z=zeros(part_no,1);  %the vertical speed
vI=zeros(part_no,1); %the impact velocity for a collision
thetaI=zeros(part_no,1); %%the impact angle for a collision
orbit=zeros(part_no,1);  %the number of orbits; each time a saltating particle does not rebound and is counted as an orbit
t_end=zeros(part_no,1);  %the time at which the saltation orbit is terminated (usually when t>t_max and the saltating particle does not rebound)
no_rebound_prob=zeros(part_no,1); %this holds the approximate probability per second that a saltating particle will settle back to the soil surface by failing to rebound upon impact
no_bounces=zeros(part_no,1); %the number of bounces the particle has undergone
no_rebounds=ones(part_no,1); %the number of rebounds the particle has made, which is closely tied to no_bounces
total_ejecta = zeros(1,part_no); %this holds the total number of particles splashed by impacting saltating particles.  This is used to determine the distribution of impact speeds, scaled by the number of ejecta that each impact produces.  This in turn is used to determine the ejection velocity of a new saltating particle.
bounce_time = zeros(part_no,1); %this keeps track of the times at which particles collide with the surface
surface_shear_stress = zeros(part_no,1);  %this keeps track of the stress exerted by particles at the surface
ejecta_matrix = zeros(part_no,part_no); %this holds the number of ejected particles, as a function of the ejecting particle type, summed over all saltation orbits
vI_distr = null(1); %holds the impact speed with which particles are ejected 
v_lift = null(1); %holds the speed with which particles leave the surface
pEj_new = null(1); %holds the ejection speed of splashed particles
bounce_length = null(1); %holds the hop lengths of the saltating particles
height = zeros(1,part_no);  %a measure of the average height of saltation trajectories
aerosol_emis = zeros(part_no,1); %this is a crude quantitative measure of dust aerosols ejected from the soil, which is scaled proportionally to the energy of impacting particles
j_end_prev = j_end; %to keep track of whether j_end should be increased when particles exceed delta_z*j_end_prev
i_h_change = zeros(part_no,1); %this signals the last i at which h was doubled, which occurs if the particle was suspended relatively high above the ground, where the time step can be much reduced to achieve a similar accuracy
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       This portion of the code calculates the actual saltation trajectories      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:part_no %cycling over all the saltating particles in the system    
    %the following sets the initial conditions of the particle
    [v_x(k,1), v_z(k,1), spin(k,1)] = init_cond(pEj(k,1:max(1,size(find(pEj(k,:)>0),2))), d(k), u_fr_surf, min_lift==k); %the initial conditions of the ejected particle are taken from the running pEj list of the speed of splashed surface particles
    v_lift(sum(orbit)+1) = sqrt(v_x(k,1)^2+v_z(k,1)^2); %the speed with which saltating particles are lifted
    z(k,1) = d(k)/2; %the center of the particle is assumed to be a distance d/2 from the surface    
    v_x_rebound(k,1) = v_x(k,1); v_z_rebound(k,1) = v_z(k,1); vR(k,1) = sqrt(v_x_rebound(k,1)^2+v_z_rebound(k,1)^2); %the rebound speed
    thetaR(k,1) = (180/pi)*atan(v_z_rebound(k,1)/v_x_rebound(k,1));  %the rebound angle                    
    if (thetaR(k,no_rebounds(k)) < 0) thetaR(k,no_rebounds(k)) = thetaR(k,no_rebounds(k)) +180; end %mapping the rebound angle between 0 and 180 degrees from horizontal
    
    %the folliwing calculates the forces on the particle
    i=1; %initializing i, which indicates the model time step for this particle trajectory simulaation
    j = round(z(k,i)/delta_z+0.5);  %j marks the vertical position of the particle on the coarse grid
    j_acc = round(z(k,i)/delta_z_acc+0.5); %j_acc marks the vertical position of the particle on the fine grid, close to the surface
    [U(k,1),dU(k,1)] = calc_windspeed (u, u_acc, d(k), z(k,1), j_acc, j_end_prev, 1); %this calculates the wind speed felt by the particle
    V_rel_x(1) = v_x(k,1)-U(k,1); V_rel_z(1) = v_z(k,1); V_rel_mag(1) = sqrt(V_rel_x(1)^2+V_rel_z(1)^2); %the magnitude of the vector difference between the wind speed and the particle speed
    Re(k,1) = (d(k)*rho/mu)*sqrt(v_z(k,1)^2+(v_x(k,1)-U(k,1))^2); %the Reynolds number
    Cd(k,1) = ((32/Re(k,1))^(2/3)+1)^(3/2); %calculating Cd for natural sand after Cheng (1997)    
    v_t(k) = Terminal_velocity(d(k)); %calculating the particle's terminal velocity
    F_z(k,1) = Fz(Cd(k,1), d(k), V_rel_x(1), V_rel_z(1), V_rel_mag(1), dU(k,1), m(k), spin(k,1)); a_z(k,1) = F_z(k,1)/m(k); %the vertical force on the particle and the resulting vertical acceleration
    F_x(k,1) = Fx(Cd(k,1),d(k),V_rel_x(1),V_rel_z(1),V_rel_mag(1),spin(k,1),dU(k,1)); a_x(k,1) = F_x(k,1)/m(k); %the horizontal force on the particle and the resulting horizontal acceleration

    %initializing various parameters
    t(k,1) = 0; %setting time to zero at the beginning of this particle's trajectory simulation
    last_bounce_time = 0; last_bounce_i(k) = 1; i_end_orbit(k) = 1; %various parameters tat keep track of the time and the time step of the last collision with the surface
    done = false; adapt_h=true; suspend = false; %'done' signals when the model is done with modeling particle motion; 'adapt_h' signals that h is set smaller because the particle just left the surface; suspend signals when a particle has been suspended
    h = h_bounce; %the time step is set much smaller when the particle just leaves the surface, because higher accuracy is needed here
    g_eff = g+6*beta/(pi*rho_p*d(k)^2); %this 'effective' gravitational constant includes the effects of both soil cohesion and gravity and is used to account for the effect of soil cohesion on splashing (see Kok, GRL, 2010)
    
    if (turbulence==true) %setting the initial turbulent wind speeds; the turbulence parameterization is desribed in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 and is also available at http://arxiv.org/abs/0901.0270
        sigma_W = sigma_W_prop*dU(k,1)*0.4*z(k,1); %the turbulence level depends on the wind shear
        sigma_U = sigma_U_prop*dU(k,1)*0.4*z(k,1); %the turbulence level depends on the wind shear
        W_turb(k,1) = sigma_W*Gaussian(0,1); %the vertical turbulent wind speed
        U_turb(k,1) = sigma_U*Gaussian(0,1); %the horizontal turbulent wind speed
    else %turbulence is set to zero when the turbulence 'switch' is off
        W_turb(k,1) = 0;
        U_turb(k,1) = 0;
    end        
    
    % The following while looop integrates particle motion forward in time until a 
    % specified number of seconds is reached and at least one particle orbit is finished
    while (done == false) 
        i=i+1; %increasing the model time step
        h_prev = h; %keeping track of h in the previous time step
        t(k,i)=t(k,i-1) + h; %stepping forward in time - h varies between h_permanent and h_bounce depending how much time resolution is needed for an accurate simulation. This is affected by particle collisions in the particle's recent history that reduce simulation accuracy and thus require a smaller time step 
        v_z(k,i) = v_z(k,i-1) + h*a_z(k,i-1);  %these values are used as predictors - they are later corrected with more accurate calculations
        v_x(k,i) = v_x(k,i-1) + h*a_x(k,i-1); %these values are used as predictors - they are later corrected with more accurate calculations
        z(k,i) = z(k,i-1) + (h/2)*(v_z(k,i-1)+v_z(k,i)); %these values are used as predictors - they are later corrected with more accurate calculations
        spin(k,i) = spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(2*tau_rot(k)))*dU(k,i-1); %these values are used as predictors - they are later corrected with more accurate calculations
        if (i>2) %calculating the particle response time.  In order to properly model the effect of turbulence on the particle trajectories, we must have h < tau_r
            tau_r(k,i) = 4*d(k)*rho_p/(3*Cd(k,i-1)*rho*V_rel_mag(i-1));
        else
            tau_r(k,i) = rho_p*d(k)^2/(18*mu); %setting tau_r to the Reynolds regime response time when v_rel is not yet determined when i==2
        end %if, calculating the particle response time.  In order to properly model the effect of turbulence on the particle trajectories, we must have h < tau_r
        if (z(k,i) > 10*v_t(k)^2/g) %In order to determine whether the particle is suspended, which is done on the next line, the mean of a large array needs to be calculated. This is computationally expensive, and thus is only done if z is above a typical length scale, defined as v_t^2/g
            while (h < 0.025*z(k,i)/mean(abs(v_z(k,max(last_bounce_i(k),i-100):i))) && 2*h < tau_r(k,i) && turbulence == true && adapt_h == false) %in this case the particle is suspended and a larger time step does not degrade the accuracy and should be used for computational efficiency
                h = 2*h; %doubling the model time step
                i_h_change(k) = i; %this signals to the numerical integration of v and z that h has changed and so a different numerical scheme has to be used
                t(k,i)=t(k,i-1) + h;  %stepping forward in time - h varies between h_permanent and h_bounce depending how much time resolution is needed for an accurate simulation. This is affected by particle collisions in the particle's recent history that reduce simulation accuracy and thus require a smaller time step 
                v_z(k,i) = v_z(k,i-1) + h*a_z(k,i-1);  %these values are used as predictors - they are later corrected with more accurate calculations
                v_x(k,i) = v_x(k,i-1) + h*a_x(k,i-1);
                z(k,i) = z(k,i-1) + (h/2)*(v_z(k,i-1)+v_z(k,i));            
            end %while, in this case the particle is suspended and a larger time step is needed
        end %%if, in order to determine whether the particle is suspended, which is done on the next line, the mean of a large array needs to be calculated. This is computationally expensive, and thus is only done if z is above a typical length scale, defined as v_t^2/g
        while ((z(k,i-1)==d(k)/2 && v_z(k,i) < 0) || h > tau_r(k,i)) %in this case the model time step is too small for the motion to be resolved and the model time step is halved and the predictor values corrected
            if (h > h_permanent) %in this case h is reduced after it was increased due to particle suspension in the while loop above
                i_h_change(k) = i; %this signals to the numerical integration of v and z that h has changed and so a different numerical scheme has to be used
            end
            h = h/2; %halving h for better time resolution
            t(k,i)=t(k,i-1) + h; %stepping forward in time - h varies between h_permanent and h_bounce depending how much time resolution is needed for an accurate simulation. This is affected by particle collisions in the particle's recent history that reduce simulation accuracy and thus require a smaller time step 
            v_z(k,i) = v_z(k,i-1) + h*a_z(k,i-1);  %these values are used as predictors - they are later corrected with more accurate calculations
            v_x(k,i) = v_x(k,i-1) + h*a_x(k,i-1);
            z(k,i) = z(k,i-1) + (h/2)*(v_z(k,i-1)+v_z(k,i));
            last_bounce_time = t(k,i);  %this signals that the Adam-Bashforth-Moulton Predictor-Corrector Method, which is used to integrate the differential equations for the particle motion and spin, can not be used up to this point
        end %in this case the model time step is too small for the motion to be resolved and the model time step is halved
        if (h<h_permanent/(2^30)) %to prevent the model from stalling when h becomes unrealistically small
            error('h is too small')
        end %if, to prevent the program from stalling when h becomes unrealistically small
        j = round(z(k,i)/delta_z+0.5);  %j marks the vertical position of the particle on the coarse grid
        j_acc = round(z(k,i)/delta_z_acc+0.5); %j_acc marks the vertical position of the particle on the fine grid, close to the surface
        if (j_acc<1 || j<1) %in this case the particle will collide with the surface, but its vertical position has not been corrected yet
            j=1; j_acc=1;            
        end %if, in this case the particle will collide with the surface, but its vertical position has not been corrected yet
        bounce = false; %this signals when a bounce has taken place in this time step

        %calculating U(k,i) and dU(k,i)
        if (j_acc<no_fine_grid*round(delta_z/delta_z_acc)) %in this case the particle is in the fine grid
            [U(k,i),dU(k,i)] = calc_windspeed (u, u_acc, d(k), z(k,i), j_acc, j_end_prev, 1);
        else %in this case the particle is above the fine grid and is treated by the coarse grid
            [U(k,i),dU(k,i)] = calc_windspeed (u, u_acc, d(k), z(k,i), j, j_end_prev, 0);
        end %if, in this case the particle is in the fine grid

        
        % The following calculates the vertical wind speed due to turbulent fluctuations, as desribed in Kok and Renno, Journal 
        % of Geophysical Research-Atmospheres, 2009 and is also available at http://arxiv.org/abs/0901.0270
        sigma_W = sigma_W_prop*dU(k,i)*0.4*z(k,i); %the turbulence level depends on the wind shear
        sigma_U = sigma_U_prop*dU(k,i)*0.4*z(k,i);        
        V_rel_x(i) = U(k,i)+U_turb(k,i-1)-v_x(k,i); %the relative velocity in the x-direction between the wind and particle speeds
        V_rel_z(i) = W_turb(k,i-1)-v_z(k,i); %the relative velocity in the z-direction between the wind and particle speeds
        V_rel_mag(i) = sqrt(V_rel_x(i)^2+V_rel_z(i)^2); %the magnitude of the vector difference between the wind speed and the particle speed                            
        if (turbulence==true) %modeling the turbulent fluctuation along the particle's trajectory, as desribed in Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 and also available at http://arxiv.org/abs/0901.0270
            if (z(k,i)>d(k)/2) %to prevent division by zero when the particle is at the surface, where Tl = 0
                if (z(k,i) < z_s/4) %calculating Tl after Leuning et al., 2000
                    theta_Tl = 0.98; a_Tl = 0.850; b_Tl = 0.41; d_Tl = -1; x_Tl = 4*z(k,i)/z_s;
                    Tl(i) = (z_s/u_fr)*((a_Tl*x_Tl+b_Tl)+d_Tl*sqrt((a_Tl*x_Tl+b_Tl)^2-4*theta_Tl*a_Tl*b_Tl*x_Tl))/2*theta_Tl;
                elseif (z(k,i) < 2*z_s)
                    theta_Tl = 0.98; a_Tl = 0.256; b_Tl = 0.40; d_Tl = 1; x_Tl = z(k,i)/z_s - 0.8;
                    Tl(i) = (z_s/u_fr)*((a_Tl*x_Tl+b_Tl)+d_Tl*sqrt((a_Tl*x_Tl+b_Tl)^2-4*theta_Tl*a_Tl*b_Tl*x_Tl))/2*theta_Tl;
                else
                    Tl(i) = 0.4*z(k,i)/(sigma_W^2/u_fr);
                end %if, calculating Tl after Leuning et al., 2000
                Tl_par(k,i)  = max(0,Tl(i)/(sqrt(1+(beta_time*V_rel_mag(i)/sigma_W)^2))); %the reduction of the Lagrangian time scale for the direction parallel to gravity following Sawford and Guest (1991), p.156
                Tl_perp(k,i) = max(0,Tl(i)/(sqrt(1+(2*beta_time*V_rel_mag(i)/sigma_U)^2))); %the reduction of the Lagrangian time scale for the direction perpendicular to gravity following Sawford and Guest (1991), p.156
                W_turb(k,i) = W_turb(k,i-1)*exp(-h/Tl_par(k,i))+sigma_W*Gaussian(0,1)*sqrt(2)*(1-exp(-sqrt(h/Tl_par(k,i)))); %the turbulent vertical wind speed, following Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 and also available at http://arxiv.org/abs/0901.0270
                U_turb(k,i) = U_turb(k,i-1)*exp(-h/Tl_perp(k,i))+sigma_U*Gaussian(0,1)*sqrt(2)*(1-exp(-sqrt(h/Tl_perp(k,i)))); %the turbulent horizontal wind speed, following Kok and Renno, Journal of Geophysical Research-Atmospheres, 2009 and also available at http://arxiv.org/abs/0901.0270
            else %in this case the particle is at the surface and the turbulence intensity is zero
                W_turb(k,i) = sigma_W*Gaussian(0,1);           
                U_turb(k,i) = sigma_U*Gaussian(0,1);            
            end %if to prevent division by zero when the particle is at the surface, where Tl = 0
        else
            W_turb(k,i)=0;
            U_turb(k,i)=0;
        end
        
        Re(k,i) = (d(k)*rho/mu)*sqrt((v_z(k,i)-W_turb(k,i))^2+(v_x(k,i)-U(k,i)-U_turb(k,i))^2); %the Reynolds number
        Cd(k,i) = ((32/Re(k,i))^(2/3)+1)^(3/2); %calculating Cd for natural sand after Cheng (1997)
        
        % The following chunk of code accounts for the impact of saltating particles on the soil surface, and calculates 
        % rebound parameters and the splashing of surface particles, following the treatment described in Kok and Renno,
        % Journal of Geophysical Research-Atmospheres, 2009 and also available at http://arxiv.org/abs/0901.0270
        if (z(k,i-1) < d(k)/2 && z(k,i-2) >= d(k)/2) %if this is true, the particle impacted the soil surface during the previous time step
            if (no_bounces(k)~=0) %finding the time at which the particle last bounced, in order to determine whether enough data points exist to use the Adam-Bashforth-Moulton Predictor-Corrector Method
                last_bounce_time = bounce_time(k,no_bounces(k));
            end %if, finding the time at which the particle last bounced, in order to determine whether enough data points exist to use the Adam-Bashforth-Moulton Predictor-Corrector Method
            no_bounces(k) = no_bounces(k)+1; %keeping track of the total number of particle bounces or hops for this particular particle bin
            del_t = h_prev*(z(k,i-2)-d(k)/2)/(z(k,i-2)-z(k,i-1)); %this is to calculate the time at which the particle struck the surface (at z = d/2)
            v_x_impact(k,no_bounces(k)) = v_x(k,i-2) + (del_t/2)*(a_x(k,i-1)+a_x(k,i-2)); %indexing the horizontal impact speed
            v_z_impact(k,no_bounces(k)) = v_z(k,i-2) + (del_t/2)*(a_z(k,i-1)+a_z(k,i-2)); %indexing the vertical impact speed
            vI(k,no_bounces(k)) = sqrt(v_x_impact(k,no_bounces(k))^2+v_z_impact(k,no_bounces(k))^2);  %indexing the total impact speed
            thetaI(k,no_bounces(k)) = (180/pi)*atan(-v_z_impact(k,no_bounces(k))/v_x_impact(k,no_bounces(k)));  %indexing the impact angle
            bounce_length(k,no_bounces(k)) = x(k,i-1)-x(k,last_bounce_i(k)); %indexing the length of the bounce
            surface_shear_stress(k) = surface_shear_stress(k) + m(k)*(v_x_impact(k,no_bounces(k))-v_x_rebound(k,no_rebounds(k))); %this keeps track of the stress exerted by particles at the surface
            bounce = true; %this signals that a bounce has taken place in this time step
            adapt_h = true; %makes the time step adaptive; that is, the time step is reduced upon impacting the surface and is increased by factors of 2 after that, unti h = h_permanent again
            a = rand(1); %a random number between 0 and 1.  
            b = P_reb*(1-exp(-gamma*vI(k,no_bounces(k)))); %the chance that a particle will rebound (see Kok and Renno, JGR-Atm., 2009). If a>b, the particle won't rebound.
            no_rebound_prob(k) = no_rebound_prob(k) + (1-b); %this holds the approximate probability per second that a saltating particle will settle back to the soil surface by failing to rebound upon impact
            
            %the following accounts for the rebound of the particle and/or the creation of a new saltating particle with new initial conditions for the case when the particle fails to rebound from the surface
            if (a<b) %if a<b, then the particle will rebound
                if (t(k,i-1)-last_bounce_time>=3*h) %in this case enough data points after a rebound exist to use the Adam-Bashforth-Moulton Predictor-Corrector Method
                    z(k,i-1) = z(k,i-2) + (h/24)*(v_z(k,i-4)-5*v_z(k,i-3)+19*v_z(k,i-2)+9*v_z(k,i-1));                
                else %in this case, not enough data points after rebound exist for the ABM method
                    z(k,i-1) = z(k,i-2) + (h/2)*(v_z(k,i-1)+v_z(k,i-2));                   
                end %if, in this case enough data points after a rebound exist to use the Adam-Bashforth-Moulton Predictor-Corrector Method
                [x(k,i-1), v_x(k,i-1), v_z(k,i-1), alpha_R(k,no_bounces(k))] = bounce_trajectory (x(k,i-2), z(k,i-1), z(k,i-2), v_x_impact(k,no_bounces(k)), v_z_impact(k,no_bounces(k)), d(k)); %bounce_trajectory calculates the rebound parameters for a given impacting particle
                z(k,i-1) = d(k)/2; %the lowest point of the particle is assumed to be when z = d/2, that is, the particle lies on the surface.
                spin(k,i-1) = Gaussian(spin_mu, spin_sigma); %reinitializing the particle spin, which is assumed uncorrelated to the particle spin at the time of impact. This is probably not accurate, but there is too little experimental data available for a more rigorous description
                no_rebounds(k) = no_rebounds(k) + 1; %the number of rebounds the particle has made, which is closely tied to no_bounces
                v_x_rebound (k,no_rebounds(k)) = v_x(k,i-1); %indexing the horizontal rebound speed
                v_z_rebound (k,no_rebounds(k)) = v_z(k,i-1); %indexing the vertical rebound speed
                vR(k,no_rebounds(k)) = sqrt(v_x_rebound(k,no_rebounds(k))^2+v_z_rebound(k,no_rebounds(k))^2); %indexing the total rebound speed
                thetaR(k,no_rebounds(k)) = (180/pi)*atan(v_z_rebound(k,no_rebounds(k))/v_x_rebound(k,no_rebounds(k)));  %indexing the rebound angle
                if (thetaR(k,no_rebounds(k)) < 0) thetaR(k,no_rebounds(k)) = thetaR(k,no_rebounds(k)) +180; end %to account for the 180-degree ambiguity in atan (arctangent)
                h = h_bounce; %setting the time step to h_bounce. h_bounce < h_permanent, as defined above, in order to provide the needed accuracy immediately after a collision, where fewer past time steps are available for an accurate implicit solver of the differential equations describing the particle's position and motion
                j_acc = 1; %the particle's position on the fine grid, which is 1 at the surface
                [U(k,i-1),dU(k,i-1)] = calc_windspeed (u, u_acc, d(k), z(k,i-1), j_acc, j_end_prev, 1);  %calculating U(k,i) and dU(k,i)
                V_rel_x(i-1) = v_x(k,i-1)-U(k,i-1); %the difference in the x-direction between the particle and wind speeds
                V_rel_z(i-1) = v_z(k,i-1); %the difference in the z-direction between the particle and wind speeds
                V_rel_mag(i-1) = sqrt(V_rel_x(i-1)^2+V_rel_z(i-1)^2); %the magnitude of the vector difference between the wind speed and the particle speed                    
                Re(k,i-1) = (d(k)*rho/mu)*sqrt(v_z(k,i-1)^2+(v_x(k,i-1)-U(k,i-1))^2); %the Reynolds number
                Cd(k,i-1) = ((32/Re(k,i-1))^(2/3)+1)^(3/2); %calculating Cd after Cheng (1997)               
                [F_z(k,i-1)] = Fz(Cd(k,i-1), d(k), V_rel_x(i-1), V_rel_z(i-1), V_rel_mag(i-1), dU(k,i-1), m(k), spin(k,i-1)); a_z(k,i-1) = F_z(k,i-1)/m(k); %the vertical force on the particle and the resulting vertical acceleration                
                F_x(k,i-1) = Fx(Cd(k,i-1),d(k),V_rel_x(i-1),V_rel_z(i-1),V_rel_mag(i-1),spin(k,i-1),dU(k,i-1)); a_x(k,i-1) = F_x(k,i-1)/m(k); %the horizontal force on the particle and the resulting horizontal acceleration
                t(k,i-1) = t(k,i-2)+del_t; %recalculating the time at i-1, defined as the time when the particle rebounds from the surface, because the particle struck the surface between i-2 and i-1
                bounce_time(k,no_bounces(k)) = t(k,i-1); %indexing the times when the particle impacts the surface
                last_bounce_i(k) = i-1; %keeping track of the model time step of the last collision with the surface. Used to determine how many past data points to use in marching the model forward in time with the differential solver
                vR_aer = vR(k,no_rebounds(k)); %used in predicting aerosol emission
                i = i-1; %rewinding the time step by 1, because the particle impacted the surface after i-2, and i-1 represents the time at which the impact takes place and is the current time step
            else %in this case the particle does not rebound, and a new particle with new initial conditions is generated from the indexed list of splashed particles from previous iterations
                orbit(k)=orbit(k)+1; %keeping track of the total number of particle trajectories that have been simulated or, equivalently, the number of times that a saltating particle was lost to the soil surface because it failed to rebound
                [v_x(k,i-1), v_z(k,i-1), spin(k,i-1)] = init_cond(pEj(k,1:max(1,size(find(pEj(k,:)>0),2))),d(k),u_fr_surf, min_lift==k); %the initial conditions of the ejected particle are taken from a running list of the speed of splashed surface particles
                v_lift(sum(orbit)+1) = sqrt(v_x(k,i-1)^2+v_z(k,i-1)^2); %indexing the speed with which new particles are lifted from the surface
                if (t(k,i-1) < t_max  && orbit(k) < 10*t_max && orbit(k) < 2*mean_weighting_factor_part(k)) %%checking if t_max, the time that saltation trajectories for each particle bin are simulated, has been reached, in which case the while(done==false) loop is terminated, and limiting the number of simulated orbits to twice the average of previous simulations in order to prevent large swings in the model
                    no_rebounds(k) = no_rebounds(k) + 1; %the number of rebounds the particle has made, which is closely tied to no_bounces
                    v_x_rebound (k,no_rebounds(k)) = v_x(k,i-1); %indexing the horizontal rebound speed
                    v_z_rebound (k,no_rebounds(k)) = v_z(k,i-1); %indexing the vertical rebound speed
                    vR(k,no_rebounds(k)) = sqrt(v_x_rebound(k,no_rebounds(k))^2+v_z_rebound(k,no_rebounds(k))^2); %indexing the total rebound speed
                    thetaR(k,no_rebounds(k)) = (180/pi)*atan(v_z_rebound(k,no_rebounds(k))/v_x_rebound(k,no_rebounds(k)));  %indexing the rebound angle                    
                    if (thetaR(k,no_rebounds(k)) < 0) thetaR(k,no_rebounds(k)) = thetaR(k,no_rebounds(k)) +180; end   %to account for the 180-degree ambiguity in atan (arctangent)
                else %in this case t_max has been reached and done is set to true so that that the while(done==false) loop is terminated
                    done = true;
                end %if, checking if the max. orbits has been reached                
                z(k,i-1) = d(k)/2; %the lowest point of the particle is assumed to be when z = d/2, that is, the particle lies on the surface.
                t(k,i-1) = t(k,i-2)+del_t; %recalculating the time at i-1, defined as the time when the particle rebounds from the surface, because the particle struck the surface between i-2 and i-1
                x(k,i-1) = x(k,i-2)+del_t*v_x(k,i-2); %also recalculating x at i-1
                h = h_bounce; %setting the time step to h_bounce. h_bounce < h_permanent, as defined above, in order to provide the needed accuracy immediately after a collision, where fewer past time steps are available for an accurate implicit solver of the differential equations describing the particle's position and motion
                j_acc = 1; %the particle's position on the fine grid, which is 1 at the surface
                [U(k,i-1),dU(k,i-1)] = calc_windspeed (u, u_acc, d(k), z(k,i-1), j_acc, j_end_prev, 1); %calculating U(k,i) and dU(k,i)
                V_rel_x(i-1) = v_x(k,i-1)-U(k,i-1); %the difference in the x-direction between the particle and wind speeds
                V_rel_z(i-1) = v_z(k,i-1); %the difference in the z-direction between the particle and wind speeds
                V_rel_mag(i-1) = sqrt(V_rel_x(i-1)^2+V_rel_z(i-1)^2); %the magnitude of the vector difference between the wind speed and the particle speed                    
                Re(k,i-1) = (d(k)*rho/mu)*sqrt(v_z(k,i-1)^2+(v_x(k,i-1)-U(k,i-1))^2); %the Reynolds number
                Cd(k,i-1) = ((32/Re(k,i-1))^(2/3)+1)^(3/2); %calculating Cd after Cheng (1997)                 
                [F_z(k,i-1)] = Fz(Cd(k,i-1), d(k), V_rel_x(i-1), V_rel_z(i-1), V_rel_mag(i-1), dU(k,i-1), m(k), spin(k,i-1)); a_z(k,i-1) = F_z(k,i-1)/m(k); %the vertical force on the particle and the resulting vertical acceleration
                F_x(k,i-1) = Fx(Cd(k,i-1),d(k),V_rel_x(i-1),V_rel_z(i-1),V_rel_mag(i-1),spin(k,i-1),dU(k,i-1)); a_x(k,i-1) = F_x(k,i-1)/m(k); %the horizontal force on the particle and the resulting horizontal acceleration
                bounce_time(k,no_bounces(k)) = t(k,i-1); %indexing the times when the particle impacts the surface
                last_bounce_i(k) = i-1; %keeping track of the model time step of the last collision with the surface. Used to determine how many past data points to use in marching the model forward in time with the differential solver
                alpha_R(k,no_bounces(k)) = 0; vR_aer = 0; %used in predicting aerosol emission
                i=i-1; %rewinding the time step by 1, because the particle impacted the surface after i-2, and i-1 represents the time at which the impact takes place and is the current time step
                i_end_orbit(k) = i; %i_end_orbit denotes the last i at which an orbit was terminated
            end %if, if a<b, then the particle will rebound            
            
            %The following accounts for the ejection of surface particles, after Kok and Renno, JGR-Atm., 2009, 
            alpha_ej(k,no_bounces(k)) = sqrt(g/g_eff)*(1-alpha_R(k,no_bounces(k)))/alpha_ej_fract; %sets the fraction of the impacting particle momentum that's spend on ejecting surface particles
            C_vej_eff = (alpha_ej(k,no_bounces(k))/a_splash)*sqrt(g_eff*d_mean); %sets the speed of ejected particles as the impacting particle speed goes to infinite (see Kok and Renno, JGR, 2008)
            ejecta_number = 0; ejecta_char = 0; %resetting ejecta_number, which indexes the total number of ejected particles, and ejecta_char, which holds the particle bin and speed of the ejected particles
            aerosol_emis(k) = max(0,0.5*m(k)*(vI(k,no_bounces(k))^2-vR_aer^2)); %this is a crude quantitative measure of dust aerosols ejected from the soil, which is scaled proportionally to the energy of impacting particles
            no_ejecta_coll = (a_splash/sqrt(d_mean*g_eff))*(d(k)./d)*vI(k,no_bounces(k)).*mass_fraction; %calculating the total number of ejected surface particles
            if (alpha_R(k,no_bounces(k))==0 && sum(no_ejecta_coll) < 1) %in this case the collision would be considered a rebound so no_ejecta_coll is set to 0
                no_ejecta_coll = 10^(-20)*ones(part_no,1); %the small number prevents division by zero
            end %if, in this case the collision would be considered a rebound so N(i) is set to 0
            rel_ejecta_coll = no_ejecta_coll./sum(no_ejecta_coll); %the relative fraction of ejected particles that is from each particle bin. So sum(rel_ejecta_coll) = 1
            if (rem(sum(no_ejecta_coll),1)>rand(1)) %calculating total number of ejected particles
                sum_ejecta_coll = ceil(sum(no_ejecta_coll)); %rounding sum_ejecta_coll up
            else
                sum_ejecta_coll = floor(sum(no_ejecta_coll)); %rounding sum_ejecta_coll down
            end %for, calculating total number of ejected particles            
            no_ejecta_coll=zeros(1,part_no); %resetting no_ejecta_coll
            for p=1:sum_ejecta_coll %cycling over total number of ejected particles and assigning them to a particle bin based on their relative probability of ejection
                a = rand(1); %generating a random number
                for s=1:part_no %cycling over the particle bins
                    if((s==1 && a<rel_ejecta_coll(s)) || (s>1 && a>sum(rel_ejecta_coll(1:s-1)) && a<sum(rel_ejecta_coll(1:s)))) %the ejected particle is assigned to a particular particle bin if the random number falls in the range of [rel_ejecta_coll(1:s-1),rel_ejecta_coll(1:s)]
                        no_ejecta_coll(s)=no_ejecta_coll(s)+1; %increasing the number of ejected particles in this bin by 1
                        ejecta_number = ejecta_number + 1; %updating the ejecta_number
                        ejecta_char(ejecta_number,1)=s; %the first row of ejecta_char holds the bin number of the ejected particle
                    end %if, the ejected particle is assigned to a particular particle bin if the random number falls in the range of [rel_ejecta_coll(1:s-1),rel_ejecta_coll(1:s)]
                end %for, cycling over the particle bins
            end %for, cycling over total number of ejected particles and assigning them to a particle bin based on their relative probability of ejection
            if (ejecta_number > 0) %in this case, there are ejecta
                for s=1:ejecta_number %looping over the ejected particles
                    P(s) = rand(1); %generating a random number
                    ejecta_char(s,2) = log(1/(1-P(s)))*C_vej_eff; %the second row of ejecta_char holds the ejected particle's speed; the ejected particle speed is generated from an exponential distribution (see Kok and Renno, JGR-Atm., 2009)
                end %for, looping over the ejected particles          
                if (part_no > 1) %this is because the column notation is slightly different for part_no > 1, which leads to a matrix multiplication error if the two cases aren't separated
                    tot_momentum = m(k)*alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)) + sum(m(ejecta_char(1:ejecta_number,1))'.*ejecta_char(1:ejecta_number,2)); %tot_momentum is used in calculating the speed of splashed particles
                    tot_energy = 0.5*m(k)*(alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)))^2 + 0.5*sum((m(ejecta_char(1:ejecta_number,1))'.*(ejecta_char(1:ejecta_number,2)).^2)); %tot_energy is used in calculating the speed of splashed particles
                else %for when part_no > 1
                    tot_momentum = m(k)*alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)) + sum(m(ejecta_char(1:ejecta_number,1)).*ejecta_char(1:ejecta_number,2)); %tot_momentum is used in calculating the speed of splashed particles
                    tot_energy = 0.5*m(k)*(alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)))^2 + 0.5*sum((m(ejecta_char(1:ejecta_number,1)).*(ejecta_char(1:ejecta_number,2)).^2)); %tot_energy is used in calculating the speed of splashed particles
                end %if, this is because the column notation is slightly different for part_no > 1, which leads to a matrix multiplication error if the two cases aren't separated                    
                counter = 0; %counts the iterations needed to have the speeds of the ejected particles converged. If this is excessively large, something is wrong and the model will terminate
                while (tot_momentum>m(k)*vI(k,no_bounces(k)) || tot_energy>0.5*m(k)*vI(k,no_bounces(k))^2) %This loop calculates the speed of ejected particles iteratively, until conservation of energy and momentum are satisfied
                    counter = counter + 1; %counts the iterations needed to have the speeds of the ejected particles converged. If this is excessively large, something is wrong and the model will terminate
                    if (counter > 1000)
                        error('something wrong with ejecta speed calculation') %in this case, counter is excessively large, indicating that something is wrong and so the model will terminate
                    end
                    for s=1:ejecta_number %looping over the ejected particles and determining their speeds from an exponential probability distribution - see Kok and Renno, JGR-Atm., 2009
                        P(s) = P(s)*rand(1); %randomly lowering the random number until convergence is achieved
                        ejecta_char(s,2) = log(1/(1-P(s)))*C_vej_eff; %the second row of ejecta_char holds the ejected particle's speed; the ejected particle speed is generated from an exponential distribution (see Kok and Renno, JGR-Atm., 2009)
                        if (part_no > 1) %this is because the column notation is slightly different for part_no > 1, which leads to a matrix multiplication error if the two cases aren't separated
                            tot_momentum = m(k)*alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)) + sum(m(ejecta_char(1:ejecta_number,1))'.*ejecta_char(1:ejecta_number,2)); %tot_momentum is used in calculating the speed of splashed particles
                            tot_energy = 0.5*m(k)*(alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)))^2 + 0.5*sum((m(ejecta_char(1:ejecta_number,1))'.*(ejecta_char(1:ejecta_number,2)).^2)); %tot_energy is used in calculating the speed of splashed particles
                        else %for when part_no > 1
                            tot_momentum = m(k)*alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)) + sum(m(ejecta_char(1:ejecta_number,1)).*ejecta_char(1:ejecta_number,2)); %tot_momentum is used in calculating the speed of splashed particles
                            tot_energy = 0.5*m(k)*(alpha_R(k,no_bounces(k))*vI(k,no_bounces(k)))^2 + 0.5*sum((m(ejecta_char(1:ejecta_number,1)).*(ejecta_char(1:ejecta_number,2)).^2)); %tot_energy is used in calculating the speed of splashed particles
                        end %if, this is because the column notation is slightly different for part_no > 1, which leads to a matrix multiplication error if the two cases aren't separated
                        if (tot_momentum<m(k)*vI(k,no_bounces(k)) && tot_energy<0.5*m(k)*vI(k,no_bounces(k))^2) %checking whether momentum and energy conservation are satisfied
                            break; %in this case, convergence is achieved and the while loop terminates
                        end %if, checking whether momentum and energy conservation are satisfied
                    end %for, looping over the ejected particles and determining their speeds from an exponential probability distribution - see Kok and Renno, JGR-Atm., 2009
                end %while, This loop calculates the speed of ejected particles iteratively, until conservation of energy and momentum are satisfied
                for s=1:ejecta_number %looping over all the ejected particles
                    total_ejecta(ejecta_char(s,1)) = total_ejecta(ejecta_char(s,1)) + 1; %keeping track of the total number of surface particles ejected from each particular particle bin
                    pEj_new(ejecta_char(s,1),total_ejecta(ejecta_char(s))) = max(0.0001,ejecta_char(s,2)); %adding the characteristics of the splashed particles to pEj_new, which will be added to pEj, the running list of characteristics of splashed particles; the ejection speed is limited to a minimum of 1e-4 m/s, in order to prevent miniscule ejection speeds that require h << h_permanent and slow down the code
                    pEj_ejector(ejecta_char(s,1),total_ejecta(ejecta_char(s))) = k; %keeping track of the impacting saltator for each ejecta. This is used later to calculate the relative probability of occurrence of each ejecta, based mainly on the relative occurrence of the impacting saltator
                    vI_distr(sum(total_ejecta)) = vI(k,no_bounces(k)); %indexing the impact speed with which particles are ejected 
                end %for, looping over all the ejected particles
                ejecta_matrix(k,:) = ejecta_matrix(k,:) + no_ejecta_coll; %this holds the number of ejected particles, as a function of the ejecting particle type, summed over all saltation orbits
            end %if, in this case, there are ejecta
        end %if, seeing if the particle will rebound
        
        %the following calculates the horizontal and vertical accelerations the particle experiences 
        if (bounce==false && i>1)  %if bounce==false, then the particle didn't bounce and the forces on it during the time step are calculated to determine its motion. The if i>1 conditions prevents i-1 = 0 If i==1, then there was only one orbit simulated, which was one with infinitesimal velocity.  
            V_rel_x(i) = v_x(k,i)-(U(k,i)+U_turb(k,i)); %the difference in the x-direction between the particle and wind speeds
            V_rel_z(i) = v_z(k,i)-W_turb(k,i); %the difference in the z-direction between the particle and wind speeds
            V_rel_mag(i) = sqrt(V_rel_x(i)^2+V_rel_z(i)^2); %the magnitude of the vector difference between the wind speed and the particle speed                            
            F_x(k,i) = Fx(Cd(k,i),d(k),V_rel_x(i),V_rel_z(i),V_rel_mag(i),spin(k,i),dU(k,i-1)); a_x(k,i) = F_x(k,i)/m(k); %the horizontal force on the particle and the resulting horizontal acceleration
            if (j_acc <= j_acc_end) %checking whether the particle is on the fine or the coarse grid
                [F_z(k,i)] = Fz(Cd(k,i), d(k), V_rel_x(i), V_rel_z(i), V_rel_mag(i), dU(k,i), m(k), spin(k,i)); %for E, j = 0 corresponds to the surface, so must use E(j+1)       
            else
                [F_z(k,i)] = Fz(Cd(k,i), d(k), V_rel_x(i), V_rel_z(i), V_rel_mag(i), dU(k,i), m(k), spin(k,i)); %for E, j = 0 corresponds to the surface, so must use E(j+1)                
            end %if, checking whether the particle is on the fine or the coarse grid
            a_z(k,i) = F_z(k,i)/m(k); %the vertical acceleration
        end %if, if bounce==false, then the particle didn't bounce and the forces on it during the time step are calculated to determine its motion. The if i>1 conditions prevents i-1 = 0 If i==1, then there was only one orbit simulated, which was one with infinitesimal velocity.  
        
        %this integrates the differential equations that describe the particle motion, such that the new position, speed, and particle spin are obtained 
        if (bounce==false)  %the new position, speed, and particle spin do not need to be calculated if the particle just bounced, in which case these parameters are set by the particle's initial conditions
            if (no_bounces(k)==0) %finding the time at which the particle last bounced, in order to determine whether enough data points exist to use the Adam-Bashforth-Moulton Predictor-Corrector Method
                last_bounce_time = 0;
            else
                last_bounce_time = bounce_time(k,no_bounces(k));
            end %if            
            %this portion of the code sets the coefficients used in the Adams-Moulton (AM) method used below. Because h is not constant, this is somewhat involved.            
            if(adapt_h==false) %computing v and x,z based on a constant model time step for the past several model time steps. If h changed for the past several time steps, then adapt_h == true and the numerical integration must necessarily be different
                if (i - i_h_change(k) < 2 && i_h_change(k)~=0) %this accounts for when the particle is (practically) suspended and the time step is increased. in this case the accuracy of the calculation is much less important
                    if (i-i_h_change(k) == 0) %h was changed in the present time step due to particle suspension, such that only data from the current and the previous time step can be used to calculate v and x, z
                        v_x(k,i) = v_x(k,i-1) + (h/2)*(a_x(k,i-1)+a_x(k,i)); %integrating a_x to obtain v_x
                        v_z(k,i) = v_z(k,i-1) + (h/2)*(a_z(k,i-1)+a_z(k,i)); %integrating a_z to obtain v_z
                        x(k,i) = x(k,i-1) + (h/2)*(v_z(k,i-1)+v_z(k,i)); %integrating v_x to obtain x
                        z(k,i) = z(k,i-1) + (h/2)*(v_z(k,i-1)+v_z(k,i)); %integrating v_z to obtain z
                        spin(k,i) = spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(2*tau_rot(k)))*dU(k,i-1); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                    else %in this case, h was changed in the previous time step due to particle suspension, such that only data from the current and the two previous time steps can be used to calculate v and x, z
                        v_x(k,i) = v_x(k,i-1) + (h/12)*(-a_x(k,i-2)+8*a_x(k,i-1)+5*a_x(k,i)); %integrating a_x to obtain v_x
                        v_z(k,i) = v_z(k,i-1) + (h/12)*(-a_z(k,i-2)+8*a_z(k,i-1)+5*a_z(k,i)); %integrating a_z to obtain v_z
                        x(k,i)   = x(k,i-1)   + (h/12)*(-v_x(k,i-2)+8*v_x(k,i-1)+5*v_x(k,i)); %integrating v_x to obtain x
                        z(k,i)   = z(k,i-1)   + (h/12)*(-v_z(k,i-2)+8*v_z(k,i-1)+5*v_z(k,i)); %integrating v_z to obtain z
                        spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(12*2*tau_rot(k)))*(-dU(k,i-2)+8*dU(k,i-1)+5*dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                    end %h was changed in the present time step due to particle suspension, such that only data from the current and the previous time step can be used to calculate v and x, z
                elseif (round((t(k,i)-last_bounce_time-h_bounce)/h_permanent) == 2) %in this case h was less than h_permanent 2 time steps ago
                    del_i = log(h_permanent/h_bounce)/log(2); %del_i calculates the number of time steps before i-2 for which t(k,i-1) - t(k,i-2-del_i) = h_permanent
                    v_x(k,i) = v_x(k,i-1) + (h/12)*(-a_x(k,i-2-del_i)+8*a_x(k,i-1)+5*a_x(k,i)); %integrating a_x to obtain v_x
                    v_z(k,i) = v_z(k,i-1) + (h/12)*(-a_z(k,i-2-del_i)+8*a_z(k,i-1)+5*a_z(k,i)); %integrating a_z to obtain v_z
                    x(k,i)   = x(k,i-1)   + (h/12)*(-v_x(k,i-2-del_i)+8*v_x(k,i-1)+5*v_x(k,i)); %integrating v_x to obtain x
                    z(k,i)   = z(k,i-1)   + (h/12)*(-v_z(k,i-2-del_i)+8*v_z(k,i-1)+5*v_z(k,i)); %integrating v_z to obtain z
                    spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(12*2*tau_rot(k)))*(-dU(k,i-2-del_i)+8*dU(k,i-1)+5*dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                elseif (round((t(k,i)-last_bounce_time-h_bounce)/h_permanent) == 3) %in this case h was less than h_permanent 3 time steps ago
                    del_i = log(h_permanent/h_bounce)/log(2); %del_i calculates the number of time steps before i-3 for which t(k,i-2) - t(k,i-3-del_i) = h_permanent
                    v_x(k,i) = v_x(k,i-1) + (h/24)*(a_x(k,i-3-del_i)-5*a_x(k,i-2)+19*a_x(k,i-1)+9*a_x(k,i)); %integrating a_x to obtain v_x
                    v_z(k,i) = v_z(k,i-1) + (h/24)*(a_z(k,i-3-del_i)-5*a_z(k,i-2)+19*a_z(k,i-1)+9*a_z(k,i)); %integrating a_z to obtain v_z
                    x(k,i)   = x(k,i-1)   + (h/24)*(v_x(k,i-3-del_i)-5*v_x(k,i-2)+19*v_x(k,i-1)+9*v_x(k,i)); %integrating v_x to obtain x
                    z(k,i)   = z(k,i-1)   + (h/24)*(v_z(k,i-3-del_i)-5*v_z(k,i-2)+19*v_z(k,i-1)+9*v_z(k,i)); %integrating v_z to obtain z
                    spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(24*2*tau_rot(k)))*(dU(k,i-3-del_i)-5*dU(k,i-2)+19*dU(k,i-1)+9*dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                else %in this case h == h_permanent for the past 3 time steps
                    v_x(k,i) = v_x(k,i-1) + (h/24)*(a_x(k,i-3)-5*a_x(k,i-2)+19*a_x(k,i-1)+9*a_x(k,i)); %integrating a_x to obtain v_x
                    v_z(k,i) = v_z(k,i-1) + (h/24)*(a_z(k,i-3)-5*a_z(k,i-2)+19*a_z(k,i-1)+9*a_z(k,i)); %integrating a_z to obtain v_z
                    x(k,i)   = x(k,i-1)   + (h/24)*(v_x(k,i-3)-5*v_x(k,i-2)+19*v_x(k,i-1)+9*v_x(k,i)); %integrating v_x to obtain x
                    z(k,i)   = z(k,i-1)   + (h/24)*(v_z(k,i-3)-5*v_z(k,i-2)+19*v_z(k,i-1)+9*v_z(k,i)); %integrating v_z to obtain z
                    spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(24*2*tau_rot(k)))*(dU(k,i-3)-5*dU(k,i-2)+19*dU(k,i-1)+9*dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                end %if, this accounts for when the particle is (practically) suspended and the time step is increased. in this case the accuracy of the calculation is much less important
            else %in this case, the model time step h is less than h_permanent, which changes the way the differential equations must be integrated
                if ((t(k,i)-last_bounce_time)/min(h,h_bounce) <= 3) %in this case, not enough time steps at the lower h have been performed for desired accuracy, and the model is not yet ready to increase h back to h_permanent; Adams-Moulton (AM) coefficients are assigned and used to compute v and x,z for an adaptive h
                    del_i = max(1,round((t(k,i)-last_bounce_time)/h_bounce)); %this calculates the number of time steps since the collision
                    a_x_AM(1:1:del_i) = a_x(k,i:-1:i+1-del_i); %determining the Adams-Moulton (AM) coefficients for integrating a_x
                    a_z_AM(1:1:del_i) = a_z(k,i:-1:i+1-del_i); %determining the Adams-Moulton (AM) coefficients for integrating a_z
                    v_x_AM(2:1:1+del_i) = v_x(k,i-1:-1:i-del_i); %determining the Adams-Moulton (AM) coefficients for integrating v_x
                    v_z_AM(2:1:1+del_i) = v_z(k,i-1:-1:i-del_i); %determining the Adams-Moulton (AM) coefficients for integrating v_z
                    dU_AM(2:1:1+del_i) = dU(k,i-1:-1:i-del_i); %determining the Adams-Moulton (AM) coefficients for integrating dU
                    if (del_i == 1) %computing v and x,z based on how many AM coefficients are available; in this case only 1
                        v_x(k,i) = v_x(k,i-1) + h*a_x_AM(1); %integrating a_x to obtain v_x
                        v_z(k,i) = v_z(k,i-1) + h*a_z_AM(1); %integrating a_z to obtain v_z
                        x(k,i)   = x(k,i-1) + (h/2)*(v_x_AM(2)+v_x(k,i)); %integrating v_x to obtain x
                        z(k,i)   = z(k,i-1) + (h/2)*(v_z_AM(2)+v_z(k,i)); %integrating v_z to obtain z
                        spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(2*2*tau_rot(k)))*(dU_AM(2)+dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                    elseif(del_i == 2) %computing v and x,z based on how many AM coefficients are available; in this case 2
                        v_x(k,i) = v_x(k,i-1) + (h/2)*(a_x_AM(2)+a_x_AM(1)); %integrating a_x to obtain v_x
                        v_z(k,i) = v_z(k,i-1) + (h/2)*(a_z_AM(2)+a_z_AM(1)); %integrating a_z to obtain v_z
                        x(k,i)   = x(k,i-1) + (h/12)*(-v_x_AM(3)+8*v_x_AM(2)+5*v_x(k,i)); %integrating v_x to obtain x
                        z(k,i)   = z(k,i-1) + (h/12)*(-v_z_AM(3)+8*v_z_AM(2)+5*v_z(k,i)); %integrating v_z to obtain z
                        spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(12*2*tau_rot(k)))*(-dU_AM(3)+8*dU_AM(2)+5*dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                    elseif(del_i == 3) %computing v and x,z based on how many AM coefficients are available; in this case 3
                        v_x(k,i) = v_x(k,i-1) + (h/12)*(-a_x_AM(3)+8*a_x_AM(2)+5*a_x_AM(1)); %integrating a_x to obtain v_x
                        v_z(k,i) = v_z(k,i-1) + (h/12)*(-a_z_AM(3)+8*a_z_AM(2)+5*a_z_AM(1)); %integrating a_z to obtain v_z
                        x(k,i)   = x(k,i-1)   + (h/24)*(v_x_AM(4)-5*v_x_AM(3)+19*v_x_AM(2)+9*v_x(k,i)); %integrating v_x to obtain x
                        z(k,i)   = z(k,i-1)   + (h/24)*(v_z_AM(4)-5*v_z_AM(3)+19*v_z_AM(2)+9*v_z(k,i)); %integrating v_z to obtain z
                        spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(24*2*tau_rot(k)))*(dU_AM(4)-5*dU_AM(3)+19*dU_AM(2)+9*dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                        if (round(h_permanent/h) > 1) 
                            h = 2*h; %if del_i == 3 and h has not yet been increased back to h_permanent, then h is doubled
                        end %if
                        a_x_AM = 0; a_z_AM = 0; v_x_AM = 0; v_z_AM = 0; %resetting the AM coefficients
                    else
                        warning('something wrong with assigning AM coefficients'); %in this case the AM coefficients could not be set, which should not happen
                    end %if, computing v and x,z based on how many AM coefficients are available
                else %in this case, enough time steps at the lower h have been performed for the desired model accuracy, and the model is ready to increase h back to h_permanent; Adams-Moulton (AM) coefficients are assigned and used to compute v and x,z for an adaptive h
                    del_i = round(log(h/h_bounce)/log(2)); %del_i calculates the number of time steps before i-1 for which t(k,i-1) - t(k,i-2-del_i) = h
                    if (del_i<0) %in this rare case h has been set smaller than h_bounce due to very strong E-fields
                        del_i=0;
                    end %if, in this rare case h has been set smaller than h_bounce due to very strong E-fields
                    a_x_AM(1:2) = a_x(k,i:-1:i-1); a_x_AM(3) = a_x(k,i-2-del_i); %determining the Adams-Moulton (AM) coefficients for integrating a_x
                    a_z_AM(1:2) = a_z(k,i:-1:i-1); a_z_AM(3) = a_z(k,i-2-del_i); %determining the Adams-Moulton (AM) coefficients for integrating a_z
                    v_x_AM(2) = v_x(k,i-1); v_x_AM(3) = v_x(k,i-2-del_i); %determining the Adams-Moulton (AM) coefficients for integrating v_x
                    v_z_AM(2) = v_z(k,i-1); v_z_AM(3) = v_z(k,i-2-del_i); %determining the Adams-Moulton (AM) coefficients for integrating v_z
                    dU_AM(2) = dU(k,i-1); dU_AM(3) = dU(k,i-2-del_i); %determining the Adams-Moulton (AM) coefficients for integrating dU
                    
                    v_x(k,i) = v_x(k,i-1) + (h/12)*(-a_x_AM(3)+8*a_x_AM(2)+5*a_x_AM(1)); %integrating a_x to obtain v_x
                    v_z(k,i) = v_z(k,i-1) + (h/12)*(-a_z_AM(3)+8*a_z_AM(2)+5*a_z_AM(1)); %integrating a_z to obtain v_z
                    x(k,i) = x(k,i-1) + (h/12)*(-v_x_AM(3)+8*v_x_AM(2)+5*v_x(k,i)); %integrating v_x to obtain x
                    z(k,i) = z(k,i-1) + (h/12)*(-v_z_AM(3)+8*v_z_AM(2)+5*v_z(k,i)); %integrating v_z to obtain z
                    spin(k,i)= spin(k,i-1)*exp(-h/tau_rot(k)) + (h/(12*2*tau_rot(k)))*(-dU_AM(3)+8*dU_AM(2)+5*dU(k,i)); %obtaining the particle spin; see equation (7) in Anderson and Hallet, 1986
                    if (round(h_permanent/h) > 1)
                        h = 2*h; %h is increased by a factor of two each loop, until h_permanent is reached
                    end
                    if (h == h_permanent) %in this case, h has reached its value when no surface collisions occurred recently
                        adapt_h = false; %setting adapt_h = false to signal that h is back to h_permanent.
                    end %if, in this case, h has reached its value when no surface collisions occurred recently
                end %if-else, assigning AM coefficients and computing v and x,z for adaptive h
            end %if-else, computing v and x,z based on a constant model time step for the past 4 model time steps. If h changed for the past 4 time steps, then adapt_h == true and the numerical integration will be different
        end %if, the new position, speed, and particle spin do not need to be calculated if the particle just bounced, in which case these parameters are set by the particle's initial conditions

        %the following checks if the particle is suspended, in which case the loop terminates for computational efficiency; the second condition shouldn't occur, but it keeps the model from going into an infinite loop
        if (z(k,i) > 10*v_t(k)^2/g || t(k,i) > 10*t_max) %In order to determine whether the particle is suspened, which is done on the next line, the mean of a large array needs to be calculated. This is computationally expensive, and thus is only done if z is above a typical length scale, defined as v_t^2/g
            if (z(k,i)/mean(abs(v_z(k,max(last_bounce_i(k),i-100):i))) > 0.25*t_max || t(k,i) > 25*t_max) %checking whether the particle is suspended
                if (last_bounce_i(k)==1||max(z(k,1:last_bounce_i(k)))<2*d(k)) %if the particle is suspended and underwent no saltation orbits, the initial conditions are reset and time is 'rewound' to zero
                    i = 1; j_acc = 1; done = false; suspend = false; t(k,1) = 0; %resetting some parameters
                    [v_x(k,1), v_z(k,1), spin(k,1)] = init_cond(pEj(k,1:max(1,size(find(pEj(k,:)>0),2))),d(k),u_fr_surf, min_lift==k); %the initial conditions of the ejected particle are taken from a running list of the speed of splashed surface particles
                    v_lift(sum(orbit)+1) = sqrt(v_x(k,1)^2+v_z(k,1)^2); %keeping track of the speed with which saltating particles are lifted
                    z(k,1) = d(k)/2; %the center of the particle is assumed to be a distance d/2 from the surface    
                    v_x_rebound(k,1) = v_x(k,1); v_z_rebound(k,1) = v_z(k,1); vR(k,1) = sqrt(v_x_rebound(k,1)^2+v_z_rebound(k,1)^2); %keeping track of the rebound speed
                    thetaR(k,1) = (180/pi)*atan(v_z_rebound(k,1)/v_x_rebound(k,1));  if (thetaR(k,no_rebounds(k)) < 0) thetaR(k,no_rebounds(k)) = thetaR(k,no_rebounds(k)) +180; end   %mapping the rebound angle between 0 and 180 degrees from horizontal
                    [U(k,1),dU(k,1)] = calc_windspeed (u, u_acc, d(k), z(k,1), j_acc, j_end_prev, 1); %this calculates the wind speed felt by the particle
                    V_rel_x(1) = v_x(k,1)-U(k,1); V_rel_z(1) = v_z(k,1); V_rel_mag(1) = sqrt(V_rel_x(1)^2+V_rel_z(1)^2); %the magnitude of the vector difference between the wind speed and the particle speed                    
                    Re(k,1) = (d(k)*rho/mu)*sqrt(v_z(k,1)^2+(v_x(k,1)-U(k,1))^2); Cd(k,1) = ((32/Re(k,1))^(2/3)+1)^(3/2); %calculating Cd for natural sand after Cheng (1997)    
                    F_z(k,1) = Fz(Cd(k,1), d(k), V_rel_x(1), V_rel_z(1), V_rel_mag(1), dU(k,1), m(k), spin(k,1)); a_z(k,1) = F_z(k,1)/m(k); %the vertical force on the particle and the resulting vertical acceleration
                    F_x(k,1) = Fx(Cd(k,1),d(k),V_rel_x(1),V_rel_z(1),V_rel_mag(1),spin(k,1),dU(k,1)); a_x(k,1) = F_x(k,1)/m(k); %the horizontal force on the particle and the resulting horizontal acceleration
                    last_bounce_time = 0; last_bounce_i(k) = 1; i_end_orbit(k) = 1; adapt_h=true; h = h_bounce; %various parameters tat keep track of the time and the time step of the last collision with the surface and other particles
                    if (turbulence==true) %setting the initial turbulent wind speeds
                        sigma_W = sigma_W_prop*dU(k,1)*0.4*z(k,1); sigma_U = sigma_U_prop*dU(k,1)*0.4*z(k,1);  %the turbulence level depends on the wind shear
                        W_turb(k,1) = sigma_W*Gaussian(0,1); U_turb(k,1) = sigma_U*Gaussian(0,1); %the vertical and horizontal turbulent wind speed
                    else %turbulence is set to zero when the turbulence 'switch' is off
                        W_turb(k,1) = 0;
                        U_turb(k,1) = 0;
                    end
                else %if the particle did undergo saltation orbits, the while loop modeling its motion is terminated for computational efficiency               
                    done = true; suspend = true; %suspend signals that i_end(k) and t_end(k) should be set to the last i and t when the particle was at the surface
                    if (i_end_orbit(k) ~= last_bounce_i(k)) %if the particle did not settle back to the soil surface at the last bounce, in which case i_end_orbit and last_bounce_i would be equal, then the current orbit is terminated
                        orbit(k) = orbit(k)+1; %%keeping track of the total number of particle trajectories that have been simulated or, equivalently, the number of times that a saltating particle was lost to the soil surface
                    end
                    if (t(k,i) > 25*t_max) %warnings so that the user knows there is a potential problem
                        warning('t(k,i) > 25*t_max; increase t_max!')
                    else
                        warning('particle suspended, or t_max is too small')
                    end
                    d(k) %outputting what the particle size of the suspended particle is
                end %if, if the particle is suspended and underwent no saltation orbits, the initial conditions are reset and time is 'rewound' to zero
            end %%checking whether the particle is suspended
        end %if, %In order to determine whether the particle is suspened, which is done on the next line, the mean of a large array needs to be calculated. This is computationally expensive, and thus is only done if z is above a typical length scale, defined as v_t^2/g
        if (abs(imag(v_z(k,i)))>0) %terminating the model when the particle speed is imaginary, which shouldn't happen
            error('v_z is imaginary!');            
        end %if        
    end %while, integrating particle motion forward in time until a specified number of seconds is reached and the present particle orbit is finished; the second condition shouldn't occur, but it keeps the model from going into an infinite loop       
    t_end(k) = t(k,i);  %the time at which the saltation orbit is terminated (since we have time start at t=0 and i=1, we have that t_end = i_end*h = t(k,i) + h.        
    if (suspend == true) %if the particle was suspended during the last orbit, then only the trajectories before it was suspended, and was thus still in saltation, are taken into account
        i_end(k) = last_bounce_i(k);
        t_end(k) = t(k,i_end(k));
    elseif (t(k,i)>t_max) %in this case the particle wasn't suspended, but t_max was exceeded
        i_end(k) = i-1;
        q_int(k,i) = 0;
    else %in this case t_max wasn't exceeded
        i_end(k) = i;
    end %if, setting i_end(k)
    if (t_end(k)==0) %something is wrong if t_end(k)==0, which crashes the model
        error('t_end is zero!');
    end
    no_rebound_prob(k)=no_rebound_prob(k)/t_end(k); %this holds the approximate probability per second that a saltating particle will settle back to the soil surface by failing to rebound upon impact
    surface_shear_stress(k)=surface_shear_stress(k)/t_end(k); %normalizing the surface shear stress exerted by particles at the surface by the simulation time
end %for, cycling over all the particles in the system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The below code processes the data from the calculation of the saltation trajectories above. %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the below is to strip off the last data point, which is erroneous
v_x_temp = v_x(1:part_no, 1:size(z,2)-1); v_x=v_x_temp; v_z_temp = v_z(1:part_no, 1:size(z,2)-1); v_z=v_z_temp;
t_temp = t(1:part_no, 1:size(z,2)-1); t=t_temp; x_temp = x(1:part_no, 1:size(z,2)-1); x=x_temp; z_temp = z(1:part_no, 1:size(z,2)-1); z=z_temp;

%calculating parameters in order to determine the particle concentration and guide the model towards steady-state
if (max(n)>0 && min(orbit) > 0 && max((t_end./orbit).*(ejecta_matrix'*n')) > 0) %calculating the number of splashed particles and the number of particles lost to the soil surface; n>0 is there to prevent division by zero
    created = (n'.*(sum(ejecta_matrix,2)./t_end))'; %number of particles created through splashing by each particle bin, per second
    lost = (n'.*orbit./t_end)'; %number of particles lost through failure to rebound by each particle bin, per second
    ejecta_orbit_ratio = sum(created)/sum(lost); %the ratio of the number of created particles to the number of lost particles. This is an importatnt that is used to adjust the particle concentration and guide the model to steady-state
    weighting_factor_part = orbit; %this parameter is averaged over the wind_profile_iterations and fed back into saltation_trajectory in order to limi the number of simulated orbits to twice the average of previous simulations in order to prevent large swings in the model
    weighting_factor_avg = (n.*m)*orbit/sum(n.*m); %this is used to weigh the results of the current calculation of the saltation trajectories by the total number of orbits simulated
    n_relative = ((t_end./orbit).*(ejecta_matrix'*n'))'/max((t_end./orbit).*(ejecta_matrix'*n')); %the relative size distribution is the total number of particles splashed of each bin, times the average lifetime of particles in each bin
    impacts = n.*no_bounces'./t_end'; %the number of impacts on the soil bed by saltating particles per m2 per second
    rebounds = impacts.*(no_bounces'-orbit')./no_bounces'; %the number of rebounding particles leaving the soil ber per m2 per second
    ejecta = impacts.*sum(ejecta_matrix,2)'./no_bounces'; %the number of particles splashed from the soil bed per m2 per second
    rebounds_ejecta = rebounds + ejecta; %the total number of particles leaving the surface per m2 per second
else %this occurs mostly for when n==0
    ejecta_orbit_ratio = 1; %the ratio of the number of created particles to the number of lost particles. This is an importatnt that is used to adjust the particle concentration and guide the model to steady-state
    weighting_factor_part = orbit; %this parameter is averaged over the wind_profile_iterations and fed back into saltation_trajectory in order to limi the number of simulated orbits to twice the average of previous simulations in order to prevent large swings in the model
    weighting_factor_avg = (mass_fraction./m)*orbit/sum(mass_fraction./m); %this is used to weigh the results of the current calculation of the saltation trajectories by the total number of orbits simulated
    n_relative = (t_end./orbit)'.*(mass_fraction./m)/max((t_end./orbit)'.*(mass_fraction./m)); %the relative size distribution is the total number of particles splashed of each bin, times the average lifetime of particles in each bin
    impacts = ones(part_no,1); rebounds_ejecta = ones(part_no,1); created = ones(part_no,1); lost = ones(part_no,1);
end %if, calculating the number of splashed particles and the number of particles lost to the soil surface; n>0 is there to prevent division by zero

surface_coll_freq = no_bounces; %the total number of surface collisions. This is normalized in saltation_model by the total simulation time of saltation trajectories in order to obtain the frequency (per second) of surface collisions
for k=1:1:part_no %finding the average length and height of each particle's saltation trajectory
    height(k)=sum(z(k,1:i_end(k)))/i_end(k); %a measure of the average height of saltation trajectories
    length(k)=x(k,i_end(k))/no_bounces(k); %a measure of the average length of saltation trajectories
    alt_length(k) = bounce_length(k,:)*vI(k,:)'/sum(vI(k,:)); %this holds an alternative measure of the average length of saltating particle trajectories that is weighted by the impact speed
    v_x_impact_temp(k) = sum(v_x_impact(k,:))/no_bounces(k); %this helps calculate the average impact speed
    v_z_impact_temp(k) = sum(v_z_impact(k,:))/no_bounces(k);
    j_end = max(j_end, round(max(z(k,1:i_end(k)))/delta_z+0.5));
end %finding the horizontal force exerted by the wind, per particle, to later calculate u
v_x_rebound = sum(v_x_rebound,2)./no_rebounds; %calculating the average rebound speed
v_z_rebound = sum(v_z_rebound,2)./no_rebounds;   
v_x_impact = v_x_impact_temp'; %calculating the average impact speed
v_z_impact = v_z_impact_temp';
h = h_permanent; %resetting h to h_permanent

aerosol_emis_tot = sum(aerosol_emis*n/t_end'); %aerosol_emis is a measure of dust aerosol emission, calculated from the total kinetic energy impacting the soil surface, and no_bounces_total holds the total number of bounces of particles in a particular particle bin
for k=1:part_no %cycling over all particle bins
    avg_vI(k) = mean(vI(k,1:no_bounces(k))); %calculating the average impact speed
    avg_thetaI(k) = mean(thetaI(k,1:no_bounces(k)).*vI(k,1:no_bounces(k)))/avg_vI(k); %calculating the average impact angle, weighted by the impact speed in order to reduce the dominance of creep impacts
end %for, cycling over all particle bins
avg_vI = mass_fraction*avg_vI'; %normalizing vI by the mass_fraction
avg_thetaI = mass_fraction*avg_thetaI'; %normalizing thetaI by the mass fraction

max_pEj = 20; %the maximum number of ejecta per particle bin saved per minor iteration
ejecta_totals = min(max_pEj,total_ejecta); %limiting the total number of ejecta per particle bin that will be indexed to 20 for computational efficiency

%the following code prepares pEj, which is a matrix holding ejecta speeds for each particle bin.
pEj_temp = zeros(part_no,1); %a temporary variable to calculate a reduced pEj_new
if (part_no>1) %need to fill pEj, which holds ejecta speeds, in a manner consistent with the occurrence of the ejecting saltator
    for k=1:part_no  %cycling over the particle bins
        clear rel_prob; %clearing rel_prob so it can be reused
        if (total_ejecta(k)>0) %checking that there are ejecta
            for p=1:total_ejecta(k) %cycling over all ejected particles in particle bin k
                k_eject(p) = pEj_ejector(k,p); %ejecting particle number
                if (max(n)>0) %setting the relative probability of occurrence of particle ejections, based on particle number and simulated time
                    rel_prob(p) = n(k_eject(p))/t_end(k_eject(p)); %calculating the relative probability of each ejection as the number of ejector particles divided by the total simulated time
                else %in this case, the particle numbers have not been properly initialized yet
                    rel_prob(p) = mass_fraction(k_eject(p))/d(k_eject(p))^4/t_end(k_eject(p)); %calculating the relative probability of each ejection as the number of ejector particles divided by the total simulated time; the max(1,n) is there in case n=0 for the first iteration
                end %if, setting the relative probability of occurrence of particle ejections, based on particle number and simulated time
            end %for, cycling over all ejected particles in particle bin k
            cum_prob = cumsum(rel_prob)/sum(rel_prob); %calculating the cumulative probability, normalized to 1
            for p=1:1:max_pEj %choosing max_pEj random splashed particles to pass on to pEj in saltation_model
                a = rand(1); %generating a random number
                pEj_temp(k,p) = pEj_new(k,find(cum_prob>a,1)); %using the random number to select an ejecta, while accounting for the relative probability of occurrence of each ejecta
            end %for, choosing 20 random splashed particles to pass on to pEj in saltation_model        
        end %if, checking that there are ejecta
    end %for, cycling over the particle bins
else %in this case there's only one particle bin so that all ejecta have equal probability of occurring
	if (total_ejecta(k)>0) %checking that there are ejecta
        for p=1:1:max_pEj %choosing max_pEj random splashed particles to pass on to pEj in saltation_model
            a = rand(1); %generating a random number
            pEj_temp(k,p) = pEj_new(k,round(a*(total_ejecta(k)-1)+1)); %using the random number to select an ejecta            
        end %for, choosing max_pEj random splashed particles to pass on to pEj in saltation_model
	end %if, checking that there are ejecta
end  %if, need to fill pEj, which holds ejecta speeds, in a manner consistent with the occurrence of the ejecting saltator
pEj = pEj_temp; %this is returned to saltation_model and incorporated in the global pEj