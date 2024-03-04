function UfrLoop (t_max, h, delta_z, u_fr, d, ufr_thr, filetext, owen, turbulence, spare)

%UfrLoop runs saltation_model over a specified range of u_fr and archives
%the result in uRun***.txt, where *** is a user-specified string

fidUfrLoop = fopen(strcat('uRun',filetext,'.txt'),'wt'); %opening a text 
fprintf(fidUfrLoop, 'ufr:      z50    vI     vx     z0s        ufr_surf log_ejecta mass       massflux   length delVx   z50_mass thetaI no_iterations \n'); %creating a header for uRun***.txt

for i = 1:1:size(u_fr,2) %cycling over the array of u_fr
    try %to catch any errors that might occur in saltation_model
        [avg_z50(i), avg_z50_mass(i), avg_vI(i), avg_thetaI(i), avg_vx(i), avg_z0s(i), avg_ufr_surf(i), avg_log_ejecta(i), avg_mass(i), avg_massflux(i), avg_length(i), avg_delVx(i), no_iterations(i)] = COMSALT_main (t_max, h, delta_z, u_fr(i), d, ufr_thr, '', [owen, turbulence], spare);
        fprintf(fidUfrLoop, '%1.4f    %1.4f %2.4f %2.4f %1.3e %1.4f   %+1.3f     %1.3e %1.3e %2.3f  %1.4f  %1.4f   %2.2f  %2.0f \n', u_fr(i), avg_z50(i), avg_vI(i), avg_vx(i), avg_z0s(i), avg_ufr_surf(i), avg_log_ejecta(i), avg_mass(i), avg_massflux(i), avg_length(i), avg_delVx(i), avg_z50_mass(i), avg_thetaI(i), no_iterations(i)); %writing the data to uRun***.txt
    catch 
        warning('error occurred in saltation_model');
    end %try, to catch any errors that might occur in saltation_model
end
fclose(fidUfrLoop); %closing uRun***.txt
