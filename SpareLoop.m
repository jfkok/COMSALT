function SpareLoop (t_max, h, delta_z, u_fr, d, ufr_thr, filetext, owen, turbulence, spare)

fidSpareLoop = fopen(strcat('SpareRun',filetext,'.txt'),'wt');
fprintf(fidSpareLoop, 'spare:      z50    vI     vx     z0s        ufr_surf log_ejecta mass       massflux   length delVx   z50_mass thetaI no_iterations \n');

for i = 1:1:size(spare,2)
    try
        [avg_z50(i), avg_z50_mass(i), avg_vI(i), avg_thetaI(i), avg_vx(i), avg_z0s(i), avg_ufr_surf(i), avg_log_ejecta(i), avg_mass(i), avg_massflux(i), avg_length(i), avg_delVx(i), no_iterations(i)] = COMSALT_main (t_max, h, delta_z, u_fr, d, ufr_thr, '', [owen, turbulence], spare(i));
        fprintf(fidSpareLoop, '%1.4f    %1.4f %2.4f %2.4f %1.3e %1.4f   %+1.3f     %1.3e %1.3e %2.3f   %1.4f  %1.4f   %2.2f  %2.0f \n', spare(i), avg_z50(i), avg_vI(i), avg_vx(i), avg_z0s(i), avg_ufr_surf(i), avg_log_ejecta(i), avg_mass(i), avg_massflux(i), hmf001z(i), avg_vertflux_ratio(i), avg_length(i), avg_alt_length(i), avg_delVx(i), avg_z50_mass(i), z50_bin_avg(i), avg_thetaI(i), avg_aer_emis(i), avg_ufrthrred(i), avg_Esurf(i), no_iterations(i));
    catch
    end    
end
fclose(fidSpareLoop);
