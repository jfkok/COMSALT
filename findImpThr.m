function [u_fr_thr,z50_avg] = findImpThr (t_max, psd_flag, u_fr_thr, spare)

%the function findImpThr determines the impact threshold for a particle
%particle size distribution. Numerical parameters used are specified below

h = 0.0075; %the model time step
delta_z = 0.001; %the box spacing on the coarse grid
Owen = 1; %Owen's hypothesis must be used to find the impact threshold (that is, the friction velocity at the surface must equal the impact threshold)
turbulence = 1; %switch for turbulence
counter = 0; %initializing the counter that keeps track of the number of successful determination of the impact threshold
counter_end = 4; %the minimum of the number of function calls to saltation_model to determine the impact threshold
tolerance = 0.005; %the minimum acceptable standard deviation of the mean of the obtained array of the impact threshold
u_fr_thr_save = 0; %initializing u_fr_thr_save, which records the impact thresholds
fidData = fopen('FindUfrThr.txt','wt'); %opening a data file
fprintf(fidData, 'Counter u_fr_thr z50\n'); %creating the header for the data file

while(counter<counter_end || err > tolerance) %looping until the impact_threshold is determined for a sufficient number of times and with sufficient precision
    try %to prevent errors in saltation_model from crashing the code
        [avg_z50, avg_z50_mass, avg_vI, avg_thetaI, avg_vx, avg_z0s, avg_ufr_surf, mean_log_ejecta, avg_mass, avg_massflux, avg_length, avg_delVx, total_iterations] = COMSALT_main (t_max, h, delta_z, 1.001*u_fr_thr, psd_flag, u_fr_thr, '', [Owen, turbulence], spare); %calling saltation_model
        u_fr_thr_old = u_fr_thr; %saving the threshold friction velocity used in the previous call to saltation_model
        if (abs(mean_log_ejecta)<0.15) %mean_log_ejecta is a measure of the convergence of the result from saltation_model to steady-state. If this is small, the result can be considered representative of steady-state saltation
            counter = counter + 1; %increasing the counter by 1
            u_fr_thr_save(counter) = u_fr_thr_old*2^(-0.1*mean_log_ejecta) %adjusting the impact threshold according to the imbalance in the ratio of the number of splashed particles to saltating particles that failed to rebound
            z50_save(counter) = avg_z50; %saving the average value of z50
            fprintf(fidData, '%3.0f %1.5f %1.5f\n', counter, u_fr_thr_save(counter), z50_save(counter)); %writing the result to a data file
            u_fr_thr = mean(u_fr_thr_save) %calculating the new impact threshold to be fed back into saltation_model 
            err=(std(u_fr_thr_save)/sqrt(counter))/mean(u_fr_thr_save) %calculating the standard deviation of the mean of u_fr_thr_save
        elseif (counter==0) %in this case no successful determination of the impact threshold has yet occurred, and the used impact threshold led to a large imbalance in the ratio of the number of splashed particles to saltating particles that failed to rebound, and u_fr_thr is adjusted to achieve a better approximation of steady state
            u_fr_thr = u_fr_thr_old*2^(-0.1*mean_log_ejecta) %adjusting u_fr_thr, the impact threshold that is used by saltation_model, based on mean_log_ejecta
        else %in this case a successful determination of the impact threshold has already occurred, and the used impact threshold led to a large imbalance in the ratio of the number of splashed particles to saltating particles that failed to rebound, and u_fr_thr is adjusted to achieve a better approximation of steady state
            u_fr_thr = (sum(u_fr_thr_save)+u_fr_thr_old*2^(-0.1*mean_log_ejecta))/(size(u_fr_thr_save,2)+1) %adjusting u_fr_thr, the impact threshold that is used by saltation_model, based on mean_log_ejecta
        end %if, mean_log_ejecta is a measure of the convergence of the result from saltation_model to steady-state. If this is small, the result can be considered representative of steady-state saltation
    catch 
        warning('error in saltation_model!');
    end %try, to prevent errors in saltation_model from crashing the code
end %while, looping until the impact_threshold is determined for a sufficient number of times and with sufficient precision
z50_avg=mean(z50_save); %calculating the average z50
u_fr_thr