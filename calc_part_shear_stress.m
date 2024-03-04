function [total_surface_shear_stress, all_particles_shear_stress, all_particles_shear_stress_acc] = calc_part_shear_stress (n_relative, d, flux_down, v_x_avg_up, v_x_avg_down, flux_down_acc, v_x_avg_down_acc, v_x_avg_up_acc)

%calc_part_shear_stress takes the output from calc_avg_velocity_density and
%calculates the particle shear stress on both the coarse and fine grids

%these are global parameters used by calc_avg_velocity_density
global pi rho_p part_no delta_z delta_z_acc no_fine_grid;

r=d/2; %the particle radius in m
m = (4/3)*pi*rho_p*r.^3; %the particle mass in kg
for k=1:1:part_no %cycling over all particle bins
    single_particle_shear_stress(k,1:1:size(flux_down,2)) = m(k)*flux_down(k,:)*n_relative(k).*(v_x_avg_down(k,:) - v_x_avg_up(k,:)); %calculating the shear stress exerted by particles from a particular particle bin on the coarse grid
    single_particle_shear_stress_acc(k,1:1:size(flux_down_acc,2)) = m(k)*flux_down_acc(k,:)*n_relative(k).*(v_x_avg_down_acc(k,:) - v_x_avg_up_acc(k,:)); %calculating the shear stress exerted by particles from a particular particle bin on the fine grid
end %for, cycling over all particle bins
all_particles_shear_stress = sum(single_particle_shear_stress,1); %calculating the particle shear stress due to all particle bins on the coarse grid
all_particles_shear_stress_acc = sum(single_particle_shear_stress_acc,1); %calculating the particle shear stress due to all particle bins on the fine grid
all_particles_shear_stress(1) = all_particles_shear_stress_acc(1); %the calculated shear stress is for the bottom of each grid box. Hence, all_particles_shear_stress(1) = all_particles_shear_stress_acc(1), and all_particles_shear_stress(2:no_fine_grid) = all_particles_shear_stress_acc(round(delta_z/delta_z_acc):round(delta_z/delta_z_acc):size(single_particle_shear_stress_acc,2)-round(delta_z/delta_z_acc))
all_particles_shear_stress(2:no_fine_grid) = all_particles_shear_stress_acc(round(delta_z/delta_z_acc):round(delta_z/delta_z_acc):size(single_particle_shear_stress_acc,2)-round(delta_z/delta_z_acc)); %the calculated shear stress is for the bottom of each grid box. Hence, all_particles_shear_stress(1) = all_particles_shear_stress_acc(1), and all_particles_shear_stress(2:no_fine_grid) = all_particles_shear_stress_acc(round(delta_z/delta_z_acc):round(delta_z/delta_z_acc):size(single_particle_shear_stress_acc,2)-round(delta_z/delta_z_acc))
total_surface_shear_stress = all_particles_shear_stress_acc(1); %the shear stress at the surface is equal to that at the lowest grid box on the fine grid
if (total_surface_shear_stress==0) %the surface shear stress can not be zero
    error('total_surface_shear_stress==0')
end %if, the surface shear stress can not be zero
