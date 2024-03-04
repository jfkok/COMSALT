function [z50, z50mass, z50_bin, L50, total_mass_flux, vert_mf, vert_mf_acc, hor_mf, z_s] = calc_mass_flux(bounce_length_total, no_bounces, hor_interval, time_end, density_down, density_up, v_x_avg_down, v_x_avg_up, density_down_acc, density_up_acc, v_x_avg_down_acc, v_x_avg_up_acc, n, d, height)

%calc_mass flux takes particle concentration and density and horizontal
%speed profiles and calculates the horizontal and vertical mass flux profiles

%these are global parameters used by calc_mass_flux
global part_no rho_p pi j_end j_acc_end delta_z delta_z_acc no_fine_grid z_s_frac Mars

m = (1/6)*pi*rho_p*d.^3; %the particle mass in kg
hor_mf = zeros (part_no, round(max(max(bounce_length_total))/hor_interval+0.5)); %initializing hor_mf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  The following calculates the horizontal profile of the particle mass flux %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:part_no %cycling over all particle bins
    for p=1:1:no_bounces(k) %cycling over all bounces made by particles in this particle bin
        j_hor_cm = round(bounce_length_total(k,p)/hor_interval+0.5); %calculating the number of hor_intervals that the bounce length corresponds to
        if (bounce_length_total(k,p)>0) %to prevent division by zero
            hor_mf(k,1:1:j_hor_cm-1) = hor_mf(k,1:1:j_hor_cm-1) + n(k)*m(k)*hor_interval/bounce_length_total(k,p); %the horizontal profile of the mass flux is distributed evenly over all increments up until j_hor_cm - 1, to reflect the randomness in the distance between the lift-off position and the (hypothetical) start of the collector position
            hor_mf(k,j_hor_cm) = hor_mf(k,j_hor_cm) + n(k)*m(k)*((bounce_length_total(k,p)-(j_hor_cm-1)*hor_interval)/bounce_length_total(k,p)); %the horizontal profile of the mass flux contribution at j_hor_cm depends on how much of that grid box is traversed
        end  %if, to prevent division by zero
    end %for, cycling over all bounces made by particles in this particle bin
    hor_mf(k,:) = hor_mf(k,:)/time_end(k); %normalizing the horizontal mass flux by the total simulated time
end %for, cycling over all particle bins
hor_mf = hor_mf/hor_interval; %the total horizontal mass flux needs to be normalized by hor_interval
total_hor_mf = sum(hor_mf,1); %obtaining the total horizontal mass flux summed over all bins
sum_total_hor_mf = sum(total_hor_mf);
for j=1:size(hor_mf,2) %calculating L50, which is the slope of the exponential decay of the horizontal profile of the mass flux with distance. This is usually a poor fit, because of the large contribution of creep for short horizontal distances; see Namikas (2003)
    cum_perc_mass_flux(j) = 100*sum(total_hor_mf(1:j))/sum_total_hor_mf; %calculating the cumulative percentage of the mass flux up to this horizontal distance
    if (j>1) %finding L50; j must be larger than 1 to prevent j-1 = 0
        if (cum_perc_mass_flux(j)>50 && cum_perc_mass_flux(j-1)<50);
            L50 = hor_interval*(j-1+(50-cum_perc_mass_flux(j-1))/(cum_perc_mass_flux(j)-cum_perc_mass_flux(j-1))); %finding L50 by interpolating between successive values of cum_perc_mass_flux once 50 % is reached
            break; %once L50 is found, the for loop is terminated
        end
    elseif (cum_perc_mass_flux(j)>50) %in case L50 is reached for j = 1;
        L50 = hor_interval*0.5;
    end %finding L50
end %for, calculating L50, which is the slope of the exponential decay of the horizontal mass flux with distance. This is usually a poor fit, because of the large contribution of creep for short horizontal distances; see Namikas (2003)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   The following calculates the vertical profile of the particle mass flux and z50 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:part_no %cycling over all particle bins
    for j_acc = 1:j_acc_end %cycling from j_acc to j_acc_end
        vert_mf_acc(k,j_acc) = (density_up_acc(k,j_acc)*v_x_avg_up_acc(k,j_acc)+density_down_acc(k,j_acc)*v_x_avg_down_acc(k,j_acc))*m(k)*n(k)/delta_z_acc; %calculating the vertical profile of the mass flux on the fine grid
    end %for cycling from j_acc to j_acc_end
    for j=1:j_end %cycling from j to j_end
        vert_mf(k,j) = (density_up(k,j)*v_x_avg_up(k,j)+density_down(k,j)*v_x_avg_down(k,j))*m(k)*n(k)/delta_z; %calculating the vertical profile of the mass flux on the coarse grid
    end %for, cycling from j to j_end
end %for, cycling over all particle bins
total_vert_mf_acc = sum(vert_mf_acc,1); %obtaining the vertical mass flux on the fine grid by summing over all particle bins
total_vert_mf = sum(vert_mf,1); %obtaining the vertical mass flux on the coarse grid by summing over all particle bins
total_mass = n*m'; %the total mass in the saltation layer
total_density_down_acc = 100*(n.*m)*density_down_acc/total_mass; %the percentage of total mass moving downward on the fine grid
total_density_up_acc = 100*(n.*m)*density_up_acc/total_mass; %the percentage of total mass moving upward on the fine grid
total_density_down = 100*(n.*m)*density_down/total_mass; %the percentage of total mass moving downward on the coarse grid
total_density_up = 100*(n.*m)*density_up/total_mass; %the percentage of total mass moving upward on the coarse grid
if (Mars == false) %calculating the j until where the vertical mass flux should be integrated, to exclude flux from suspended particles
    j_mf_end = min(j_end,ceil(0.5/delta_z)); %to limit the vertically integrated mass flux to the saltation layer and exclude flux from suspended particles at large shear velocities. The 0.35 meters is based on the height of vertical collectors of Namikas (1999, 2003, 2006). This constraint does not apply for other planets such as Mars where the saltation layer is much deeper.
else
    j_mf_end = j_end; %for Mars, j_me_end = j_end
end %if, calculating the j until where the vertical mass flux should be integrated, to exclude flux from suspended particles
total_mass_flux_acc = sum(total_vert_mf_acc)*delta_z_acc; %calculating the total mass flux on the fine grid
total_mass_flux = total_mass_flux_acc + sum(total_vert_mf(no_fine_grid+1:j_mf_end))*delta_z; %calculating the total mass flux   
for j_acc=1:j_acc_end %for, finding z50, z_s, and z50mass if they are reached on the fine grid, which does not usually occur
    cum_perc_mass_flux_acc(j_acc) = 100*delta_z_acc*sum(total_vert_mf_acc(1:j_acc))/total_mass_flux; %the cumulative percentage horizontal mass flux up until this point on the fine grid
    cum_perc_mass_acc(j_acc) = sum(total_density_down_acc(1:j_acc))+sum(total_density_up_acc(1:j_acc)); %the cumulative percentage saltating mass up until this point on the fine grid
    if (j_acc>1) %to prevent j_acc-1 == 0
        if (cum_perc_mass_flux_acc(j_acc)>50 && cum_perc_mass_flux_acc(j_acc-1)<50); %in this case, cum_perc_mass_flux_acc just passed the 50% mark
            z50 = delta_z_acc*(j_acc-1+(50-cum_perc_mass_flux_acc(j_acc-1))/(cum_perc_mass_flux_acc(j_acc)-cum_perc_mass_flux_acc(j_acc-1))); %calculating z50, the height below which half of the total mass flux is contained
        end %if, z50
        if (cum_perc_mass_flux_acc(j_acc)>100*z_s_frac && cum_perc_mass_flux_acc(j_acc-1)<100*z_s_frac); %in this case, cum_perc_mass_flux_acc just passed the z_s_frac mark, which is defined in load_parameters
            z_s = delta_z_acc*(j_acc-1+(99-cum_perc_mass_flux_acc(j_acc-1))/(cum_perc_mass_flux_acc(j_acc)-cum_perc_mass_flux_acc(j_acc-1))); %calculating z_s
        end %if, z_s
        if (cum_perc_mass_acc(j_acc)>50 && cum_perc_mass_acc(j_acc-1)<50) %in this case, cum_perc_mass_acc just passed the 50% mark
            z50mass = delta_z_acc*(j_acc-1+(50-cum_perc_mass_acc(j_acc-1))/(cum_perc_mass_acc(j_acc)-cum_perc_mass_acc(j_acc-1))); %calculating z50mass, the height below which half of the total saltating mass is contained
        end %if, z50mass
    else %in this case j_acc==1
        if (cum_perc_mass_flux_acc(1)>50) %in case z50 is reached for j_acc = 1
            z50 = delta_z_acc*0.5;
        end %if
        if (cum_perc_mass_flux_acc(1)>100*z_s_frac) %in case z_s is reached for j_acc = 1
            z_s = delta_z_acc*0.5;
        end %if
        if (cum_perc_mass_acc(1)>50) %in case z50mass is reached for j_acc = 1
            z50mass = delta_z_acc*0.5;
        end %if
    end %if, to prevent j_acc-1 == 0
end %for, finding z50, z_s and z50mass if they are reached on the fine grid, which does not usually occur
cum_perc_mass_flux(no_fine_grid) = cum_perc_mass_flux_acc(j_acc_end); %setting cum_perc_mass_flux(no_fine_grid) equal to that at the top of the fine grid
cum_perc_mass(no_fine_grid) = cum_perc_mass_acc(j_acc_end); %setting cum_perc_mass(no_fine_grid) equal to that at the top of the fine grid
for j=no_fine_grid+1:j_end %for, finding z50, z_s, and z50mass
    cum_perc_mass_flux(j) = (100*total_mass_flux_acc+100*delta_z*sum(total_vert_mf(no_fine_grid+1:j)))/total_mass_flux; %the cumulative percentage horizontal mass flux up until this point on the coarse grid
    cum_perc_mass(j) = cum_perc_mass(no_fine_grid)+sum(total_density_down(no_fine_grid+1:j))+sum(total_density_up(no_fine_grid+1:j)); %the cumulative percentage saltating mass up until this point on the coarse grid
    if (cum_perc_mass_flux(j)>50 && cum_perc_mass_flux(j-1)<50); %in this case, cum_perc_mass_flux just passed the 50% mark
        z50 = delta_z*(j-1+(50-cum_perc_mass_flux(j-1))/(cum_perc_mass_flux(j)-cum_perc_mass_flux(j-1))); %calculating z50, the height below which half of the total mass flux is contained
    end %if
    if (cum_perc_mass_flux(j)>100*z_s_frac && cum_perc_mass_flux(j-1)<100*z_s_frac); %in this case, cum_perc_mass_flux just passed the z_s_frac mark, which is defined in load_parameters
        z_s = delta_z*(j-1+(99-cum_perc_mass_flux(j-1))/(cum_perc_mass_flux(j)-cum_perc_mass_flux(j-1))); %calculating z_s
    end %if
    if (cum_perc_mass(j)>50 && cum_perc_mass(j-1)<50) %in this case, cum_perc_mass just passed the 50% mark
        z50mass = delta_z*(j-1+(50-cum_perc_mass(j-1))/(cum_perc_mass(j)-cum_perc_mass(j-1))); %calculating z50mass, the height below which half of the total saltating mass is contained
    end %if
    if (cum_perc_mass_flux(j)>max(50,100*z_s_frac) && cum_perc_mass(j)>50) %breaking if z_s, z50, and z50mass have been determined
        break;
    end %if, breaking if z_s, z50, and z50mass have been determined
end %for, finding z50 and z50mass

for k=1:1:part_no %cycling over all particle bins to find z50_bin
    for j=1:j_end %for, finding z50
        if (sum(vert_mf(k,1:j_end))>0) %to prevent division by zero
            cum_perc_mass_flux(j) = 100*sum(vert_mf(k,1:j))/sum(vert_mf(k,1:j_end)); %the cumulative percentage horizontal mass flux for this particle bin up until this point on the coarse grid
            if (j>1) %to prevent j-1 == 0
                if (cum_perc_mass_flux(j)>50 && cum_perc_mass_flux(j-1)<50); %in this case, cum_perc_mass_flux just passed the 50 % mark for this particle bin
                    z50_bin(k) = delta_z*(j-1+(50-cum_perc_mass_flux(j-1))/(cum_perc_mass_flux(j)-cum_perc_mass_flux(j-1))); %calculating z50_bin
                    break; %breaking if z50_bin has been determined for this particle bin
                end
            elseif (cum_perc_mass_flux(j)>50) %in case z50_bin has been reached for j == 1
                z50_bin(k) = delta_z*0.5; %calculating z50_bin
                break; %breaking if z50_bin has been determined for this particle bin
            end %if, to prevent j-1 == 0
        else
            z50_bin(k) = delta_z*0.5; %calculating z50_bin
            break; %breaking if z50_bin has been determined for this particle bin
        end %to prevent division by zero
    end %for, finding z50
end %for, cycling over all particle bins to find z50_bin