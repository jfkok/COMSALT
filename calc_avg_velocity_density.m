function [total_v_x_avg, total_v_x_avg_acc, t_down, t_up, t_down_acc, t_up_acc, v_x_avg_up, v_x_avg_down, v_z_avg_up, v_z_avg_down, density_up, density_down, flux_down, total_mass_flux, flux_down_acc, v_x_avg_down_acc, v_x_avg_up_acc, v_z_avg_down_acc, v_z_avg_up_acc, density_up_acc, density_down_acc] = calc_avg_velocity_density (d, x, z, t, v_x, v_z, time_end, i_end, n)

%calc_avg_velocity_density takes the output from saltation_trajectory and
%extracts information on the particle density, flux, and speed.

%these are global parameters used by calc_avg_velocity_density
global delta_z delta_z_acc no_fine_grid;
global part_no j_end rho_p j_acc_end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         The following initializes various parameters         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_x_avg = zeros(part_no, 1); %the average horizontal particle speed
v_x_avg_down = zeros(part_no, j_end); %the average horizontal particle speed of downward-moving particles on the coarse grid
v_x_avg_down_acc = zeros(part_no, j_acc_end); %the average horizontal particle speed of downward-moving particles on the fine grid
v_x_avg_up = zeros(part_no, j_end); %the average horizontal particle speed of upward-moving particles on the coarse grid
v_x_avg_up_acc = zeros(part_no, j_acc_end); %the average horizontal particle speed of upward-moving particles on the fine grid
v_z_avg_down = zeros(part_no, j_end); %the average vertical particle speed of downward-moving particles on the coarse grid
v_z_avg_down_acc = zeros(part_no, j_acc_end); %the average vertical particle speed of downward-moving particles on the fine grid
v_z_avg_up = zeros(part_no, j_end); %the average vertical particle speed of upward-moving particles on the coarse grid
v_z_avg_up_acc = zeros(part_no, j_acc_end); %the average vertical particle speed of upward-moving particles on the fine grid
flux_counter_down = zeros(part_no, j_end); %this counts the flux of downward moving particle passing a given grid box on the coarse grid
flux_counter_down_acc = zeros(part_no, j_acc_end); %this counts the flux of downward moving particle passing a given grid box on the fine grid
flux_counter_up = zeros(part_no, j_end); %this counts the flux of upward moving particle passing a given grid box on the coarse grid
flux_counter_up_acc = zeros(part_no, j_acc_end); %this counts the flux of upward moving particle passing a given grid box on the fine grid
t_down = zeros(part_no, j_end); %keeps track of the total time spent by downward moving particles in a given grid box on the coarse grid
t_down_acc = zeros(part_no, j_acc_end); %keeps track of the total time spent by downward moving particles in a given grid box on the fine grid
t_up = zeros(part_no, j_end); %keeps track of the total time spent by upward moving particles in a given grid box on the coarse grid
t_up_acc = zeros(part_no, j_acc_end); %keeps track of the total time spent by upward moving particles in a given grid box on the fine grid
density_down=zeros(part_no,j_end); %initializing the concentration matrix, describing the likelihood that the particle is at a particular height.
density_up=zeros(part_no,j_end); %initializing the concentration matrix, describing the likelihood that the particle is at a particular height.
density_down_acc=zeros(part_no,j_acc_end); %initializing the concentration matrix, describing the likelihood that the particle is at a particular height.
density_up_acc=zeros(part_no,j_acc_end); %initializing the concentration matrix, describing the likelihood that the particle is at a particular height.
r=d/2; %the particle radius in m
m = (4/3)*pi*rho_p*r.^3; %the particle mass in kg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following for loops cycles over the particle's trajectories and indexes their speed and the time spent in each particle bin %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1:part_no %cycling over all particle bins
    j(1:1:size(z,2)) = round(z(k,:)/delta_z+0.5); %the position of the particle on the coarse grid
    j_acc(1:1:size(z,2)) = round(z(k,:)/delta_z_acc+0.5); %the position of the particle on the fine grid
    speed_z = v_z(k,1); %the particle speed in the z-direction
    speed_x = v_x(k,1); %the particle speed in the x-direction
    h_prev = 0; %h_prev = 0 signals no collision with the surface occurred during this time step, while h_prev ~= 0 signals a collision with the surface, and communicates the remainder step size to the next loop    
    for i=2:i_end(k)
        if (k==2 && i==138)
            1;
        end
        speed_change = false;  
        speed_z_prev = speed_z; %saving the previous particle speed in the z-direction
        speed_x_prev = speed_x; %saving the previous particle speed in the z-direction
        if (t(k,i-1)~=t(k,i)) %data points are skipped in the rare instance when t(k,i-1) == t(k,i)
            speed_z = (z(k,i) - z(k,i-1))/(t(k,i) - t(k,i-1)); %the particle speed in the z-direction
            speed_x = (x(k,i) - x(k,i-1))/(t(k,i) - t(k,i-1)); %the particle speed in the x-direction
            box_acc=0; %the list of grid points on the fine scale that the particle traverses in this time step
                if (speed_z>=0)   %seeing if the particle moved up or down (a particle can move up even when v_z(i) < 0, if it just reached the top of its orbit
                    box = j(i-1):1:j(i); %the list of grid points on the coarse scale that the particle traverses in this time step
                    z_down = z(k,i-1); %the lowest point of the trajectory between the previous and the current time step
                    z_up = z(k,i); %the highest point of the trajectory between the previous and the current time step
                    h_step = t(k,i)-t(k,i-1); %the model time step
                	collision = false; %keeping track of whether a collision occurred in this time step - in this case it didn't
                    if (z_down < j_acc_end*delta_z_acc) %seeing if the particle's trajectory between the previous and the current time step is at least partially on the fine grid
                        box_acc = j_acc(i-1):1:min(j_acc(i),j_acc_end); %the list of grid points on the fine scale that the particle traverses in this time step
                    end %if, seeing if the particle's trajectory between the previous and the current time step is at least partially on the fine grid
                else %particle is moving down
                    box = j(i-1):-1:j(i); %the list of grid points on the coarse scale that the particle traverses in this time step
                    z_down = z(k,i); %the lowest point of the trajectory between the previous and the current time step
                    z_up = z(k,i-1); %the highest point of the trajectory between the previous and the current time step
                    h_step = t(k,i)-t(k,i-1); %the model time step
                    if (speed_z_prev > 0) %in this case the particle reached the top of the trajectory
                        speed_change = true;
                    end %if, %in this case the particle reached the top of the trajectory                
                    collision = false; %keeping track of whether a collision occurred in this time step - in this case it didn't
                    if (z_down < j_acc_end*delta_z_acc) %seeing if the particle's trajectory between the previous and the current time step is at least partially on the fine grid
                        box_acc = min(j_acc(i-1),j_acc_end):-1:j_acc(i); %the list of grid points on the fine scale that the particle traverses in this time step
                    end %if, seeing if the particle's trajectory between the previous and the current time step is at least partially on the fine grid
                end %if, seeing if the particle moved up or down     
                if (h_prev ~= 0) %in this case the particle just underwent a collision!  This accounts for the rebound portion
                    box = 1:1:j(i); %the list of grid points on the coarse scale that the particle traverses in this time step
                    z_down = 0; %the lowest point of the trajectory between the previous and the current time step
                    z_up = z(k,i); %the highest point of the trajectory between the previous and the current time step
                    h_step = t(k,i)-t(k,i-1); %the model time step
                    h_prev = 0;   %h_prev = 0 signals no collision with the surface occurred during this time step
                    collision = true; %keeping track of whether a collision occurred in this time step - in this case it did!
                    box_acc = 1:1:min(j_acc(i),j_acc_end); %the list of grid points on the fine scale that the particle traverses in this time step
                end  %if, in this case the particle just underwent a collision!  This accounts for the rebound portion
                if ((v_z(k,i-1) < 0 && v_z(k,i) > 0) && z(k,i) == d(k)/2) %in this case the particle just underwent a collision!  This accounts for the impact portion
                    box = j(i-1):-1:1; %the list of grid points on the coarse scale that the particle traverses in this time step
                    z_down = 0; %the lowest point of the trajectory between the previous and the current time step
                    z_up = z(k,i-1); %the highest point of the trajectory between the previous and the current time step
                    h_step = t(k,i)-t(k,i-1); %the model time step
                    h_prev = h_step;   %h_prev ~= 0 signals a collision with the surface, and communicates the remainder step size to the next loop
                    collision = true; %keeping track of whether a collision occurred in this time step - in this case it did!
                    speed_x = v_x(k,i-1); %setting the particle speed in the z-direction to that immediately after the collision
                    speed_z = v_z(k,i-1); %setting the particle speed in the x-direction to that immediately after the collision
                    box_acc = min(j_acc(i-1),j_acc_end):-1:1; %the list of grid points on the fine scale that the particle traverses in this time step
                end  %if, in this case the particle just underwent a collision!  This accounts for the impact portion                
                if (i==2) %accounting for the fact that the particle just left the surface at i==2, which counts as a collision
                    collision = true;
                end %if accounting for the fact that the particle just left the surface at i==2, which counts as a collision
                if (z_up==0 && z_down == 0) %to account for the rare case of two successive, very low, bounces
                    z_up = 10^(-6);
                end %if, to account for the rare case of two successive, very low, bounces
            delta_time = 0; delta_time_acc = 0; %resetting delta_time, which holds the list of times that the particle spent in each grid box
            clear speed_x_add; clear speed_z_add; %erasing speed_x_add and speed_z_add, which hold the list of interpolated particle speeds in each grid box the particle traverses between the past and the current model time step
            
            %the below indexes the speed and time spent in grid boxes on the coarse grid
            for box_no=1:1:size(box,2)   %cycling over the boxes that the particle was in between the last time step and the current one
                if (box_no==1)  %this was the first box the particle was in, and it thus didn't traverse the full length of it
                    if (speed_z>0)  %in this case the particle moved upwards                        
                        delta_time(box_no) = h_step*(min(box(box_no)*delta_z,z_up)-z_down)/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                    end %if, in this case the particle moved upwards
                    if (speed_z<0) %in this case the particle moved downwards
                        delta_time(box_no) = h_step*(z_up-(box(box_no)-1)*delta_z)/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                    end %if, in this case the particle moved downwards
                elseif (box_no==size(box,2))  %this is the last box the particle was in, and it thus didn't traverse the full length of it
                    if (speed_z>0)  %in this case the particle moved upwards
                        delta_time(box_no) = h_step*(z_up-(box(box_no)-1)*delta_z)/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                    end %if, in this case the particle moved upwards
                    if (speed_z<0) %in this case the particle moved downwards
                        delta_time(box_no) = h_step*(box(box_no)*delta_z-z_down)/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                    end %if, in this case the particle moved downwards
                elseif (z_up-z_down > delta_z) %in this case, the particle traverses the whole box
                    delta_time(box_no) = h_step*delta_z/(z_up-z_down); %the time the particle spends in the box is the length that it traverses, divided by its speed
                else %if none of the previous conditions are met, then the particle is not properly indexed. This should never happen
                    warning('problem with box assigning loop');
                end %if-elseif loop    
                if (size(box,2)==1) %in this case the particle doesn't leave the box
                    delta_time(box_no) = h_step; %the time the particle spends in the box is in this case the entire time step
                end
                if (speed_z < 0) %setting density and average velocities for a downward-moving particle
                    if (box_no > 1 || speed_change == true) %to prevent double-counting of particles with regards to the average particle speed
                        if (z(k,i) == d(k)/2) %assigning speed_x_add and speed_z_add; in this case the particle is at the surface and thus just underwent a collision
                            if (i>3 && z(k,i-2) ~= d(k)/2 && z(k,i-1) ~= d(k)/2) %extrapolating from v_x(k,i-2) and v_x(k,i-1), because v_x(k,i) describes the velocity upon rebound
                                time_fraction = max(0,min(1,((delta_z*(box(box_no))-z_up)/speed_z)/h_step)); %the fraction of the time step (h) it takes from the previous time (i-1) to reach the grid box in question. Used to interpolate the particle speed in the particular grid boxes on the next 2 lines
                                speed_x_add(k,box(box_no)) = v_x(k,i-1)+time_fraction*(v_x(k,i-1)-v_x(k,i-2)); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                speed_z_add(k,box(box_no)) = v_z(k,i-1)+time_fraction*(v_z(k,i-1)-v_z(k,i-1)); %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                            else %the particle speed can not be extrapolated from v_x(k,i-2_ because a collision occurred at that time
                                speed_x_add(k,box(box_no)) = v_x(k,i-1); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                speed_z_add(k,box(box_no)) = v_z(k,i-1); %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                            end
                        elseif (box_no == 1) %in this case box_no==1, but the particle speed is still indexed because the particle just underwent a change in the direction of its velocity, meaning that its speed in the same box in the previous time step was indexed as an upward-moving speed
                            speed_x_add(k,box(box_no)) = v_x(k,i-1); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                            speed_z_add(k,box(box_no)) = v_z(k,i-1); %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                        elseif (speed_change == true) %in this case the particle just underwent a change in the direction of its velocity, probably because it reached the top of its trajectory, which affects the interpolation of the particle speed in the specific grid boxes
                            speed_x_add(k,box(box_no)) = (v_x(k,i-1)+v_x(k,i))/2; %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                            speed_z_add(k,box(box_no)) = (v_z(k,i-1)+v_z(k,i))/2; %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                        else %this is the 'normal' case when nothing unusual happened (i.e., no collisions or speed changes), and the speed of the particle in each grid box is interpolated using its position in the previous and current time steps
                            time_fraction = ((delta_z*(box(box_no))-z_up)/speed_z)/h_step; %the fraction of the time step (h) it takes from the previous time (i-1) to reach the grid box in question. Used to interpolate the particle speed in the particular grid boxes on the next 2 lines
                            speed_x_add(k,box(box_no)) = v_x(k,i-1)+time_fraction*(v_x(k,i)-v_x(k,i-1)); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                            speed_z_add(k,box(box_no)) = v_z(k,i-1)+time_fraction*(v_z(k,i)-v_z(k,i-1)); %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                        end %assigning speed_x_add and speed_z_add
                        flux_counter_down(k,box(box_no)) = flux_counter_down(k,box(box_no)) + 1; %adding the current particle to the count of the flux of downward moving particle passing bin j on the coarse grid
                        v_x_avg_down(k,box(box_no)) = v_x_avg_down(k,box(box_no)) + speed_x_add(k,box(box_no)); %adding the horizontal speed of the current downward-moving particle to the sum of the horizontal downward-moving particle speeds in this box on the coarse grid
                        v_z_avg_down(k,box(box_no)) = v_z_avg_down(k,box(box_no)) + speed_z_add(k,box(box_no)); %adding the vertical speed of the current downward-moving particle to the sum of the vertical downward-moving particle speeds in this box on the coarse grid
                    end %if, to prevent double-counting of particles with regards to the average particle speed
                    density_down(k,box(box_no)) = density_down(k,box(box_no)) + delta_time(box_no); %calculating the particle density based on the time the particle spends in each grid box
                else %in this case the particle is moving upwards                
                    if (box_no > 1 || collision == true) %to prevent double-counting of particles with regards to the average particle speed
                        if (z(k,i) == d(k)/2) %assigning speed_x_add and speed_z_add; in this case the particle is at the surface and thus just underwent a collision
                            speed_x_add(k,box(box_no)) = v_x(k,i); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                            speed_z_add(k,box(box_no)) = v_z(k,i); %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                        elseif (box_no == 1) %in this case box_no==1, but the particle speed is still indexed because the particle just underwent a change in the direction of its velocity, meaning that its speed in the same box in the previous time step was indexed as an upward-moving speed
                            speed_x_add(k,box(box_no)) = v_x(k,i-1); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                            speed_z_add(k,box(box_no)) = v_z(k,i-1); %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                        elseif (speed_change == true) %in this case the particle just underwent a change in the direction of its velocity, probably because it reached the top of its trajectory, which affects the interpolation of the particle speed in the specific grid boxes
                            speed_x_add(k,box(box_no)) = (v_x(k,i-1)+v_x(k,i))/2; %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                            speed_z_add(k,box(box_no)) = (v_z(k,i-1)+v_z(k,i))/2;    %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed                         
                        else %this is the 'normal' case when nothing unusual happened (i.e., no collisions or speed changes), and the speed of the particle in each grid box is interpolated using its position in the previous and current time steps
                            time_fraction = min(0,max(1,((delta_z*(box(box_no)-1)-z_down)/speed_z)/h_step)); %the fraction of the time step (h) it takes from the previous time (i-1) to reach the grid box in question. Used to interpolate the particle speed in the particular grid boxes on the next 2 lines
                            speed_x_add(k,box(box_no)) = v_x(k,i-1)+time_fraction*(v_x(k,i)-v_x(k,i-1)); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                            speed_z_add(k,box(box_no)) = v_z(k,i-1)+time_fraction*(v_z(k,i)-v_z(k,i-1)); %speed_z_add holds the list of interpolated vertical particle speeds in each grid box the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                        end %if, assigning speed_x_add and speed_z_add
                        flux_counter_up(k,box(box_no)) = flux_counter_up(k,box(box_no)) + 1; %adding the current particle to the count of the flux of upward moving particle passing bin j on the coarse grid
                        v_x_avg_up(k,box(box_no)) = v_x_avg_up(k,box(box_no)) + speed_x_add(k,box(box_no)); %adding the horizontal speed of the current upward-moving particle to the sum of the horizontal upward-moving particle speeds in this box on the coarse grid
                        v_z_avg_up(k,box(box_no)) = v_z_avg_up(k,box(box_no)) + speed_z_add(k,box(box_no)); %adding the vertical speed of the current upward-moving particle to the sum of the vertical upward-moving particle speeds in this box on the coarse grid
                    end %to prevent double-counting                 
                    density_up(k,box(box_no)) = density_up(k,box(box_no)) + delta_time(box_no); %calculating the particle density based on the time the particle spends in each grid box
                end %if, setting density and average velocities
            end %for, cycling over the boxes that the particle was in between the last time step and the current one
            if (abs(1-sum(delta_time)/h_step)>1e-6) %in this case the individual delta_times do not sum up to h_step, which is indicative of a problem.
                warning('problem with box assigning loop');
            end %if, in this case the individual delta_times do not sum up to h_step, which is indicative of a problem.

            %the below indexes the speed and time spent in grid boxes on the fine grid
            if (box_acc~=0) %if box_acc ~= 0, the particle spend time in the bottom two grid boxes
                delta_time_acc_tot = h_step*(min(j_acc_end*delta_z_acc,z_up)-z_down)/(z_up-z_down); %the total time that the particle spends on the fine grid in this time step
                clear speed_x_add_acc; clear speed_z_add_acc; %erasing speed_x_add_acc and speed_z_add_acc, which hold the list of interpolated particle speeds in each grid box on the fine grid that the particle traverses between the past and the current model time step
                for box_no_acc=1:1:size(box_acc,2)   %cycling over the boxes that the particle was in between the last time step and the current one
                    if (box_no_acc==1)  %this was the first box the particle was in, and it thus didn't traverse the full length of it
                        if (speed_z>0)  %in this case the particle moved upwards
                            delta_time_acc(box_no_acc) = h_step*(min(box_acc(box_no_acc)*delta_z_acc,z_up)-z_down)/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                        end %if, in this case the particle moved upwards
                        if (speed_z<0) %in this case the particle moved downwards
                            delta_time_acc(box_no_acc) = h_step*(min(z_up,no_fine_grid*delta_z)-max(z_down,(box_acc(box_no_acc)-1)*delta_z_acc))/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                        end %if, in this case the particle moved downwards
                    elseif (box_no_acc==size(box_acc,2))  %this is the last box the particle was in, and it thus didn't traverse the full length of it
                        if (speed_z>0) %in this case the particle moved upwards
                            delta_time_acc(box_no_acc) = h_step*(min(z_up,no_fine_grid*delta_z)-(box_acc(box_no_acc)-1)*delta_z_acc)/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                        end %if, in this case the particle moved upwards
                        if (speed_z<0) %in this case the particle moved downwards
                            delta_time_acc(box_no_acc) = h_step*(box_acc(box_no_acc)*delta_z_acc-z_down)/(z_up-z_down);  %the time the particle spends in the box is the length that it traverses, divided by its speed
                        end %if, in this case the particle moved downwards
                    elseif (max(z_up,no_fine_grid*delta_z)-z_down > delta_z_acc) %in this case, the particle traverses the whole box
                        delta_time_acc(box_no_acc) = h_step*delta_z_acc/(z_up-z_down); %the time the particle spends in the box is the length that it traverses, divided by its speed
                    else %if none of the previous conditions are met, then the particle is not properly indexed. This should never happen
                        warning('problem with box assigning loop');
                    end %if-elseif loop            
                    if (size(box_acc,2)==1 && z_up<=no_fine_grid*delta_z) %in this case the particle doesn't leave the box
                        delta_time_acc(box_no_acc) = h_step; %the time the particle spends in the box is in this case the entire time step
                    end %in this case the particle doesn't leave the box
                    if (speed_z < 0) %setting density and average velocities for a downward-moving particle on the fine grid
                        if (box_no_acc > 1 || speed_change == true || z_up > j_acc_end*delta_z_acc) %to prevent double-counting of particles with regards to the average particle speed
                            if (z(k,i) == d(k)/2) %assigning speed_x_add_acc and speed_z_add_acc; in this case the particle is at the surface and thus just underwent a collision
                                if (i>3 && z(k,i-2) ~= d(k)/2 && z(k,i-1) ~= d(k)/2) %extrapolating from v_x(k,i-2) and v_x(k,i-1), because v_x(k,i) describes the velocity upon rebound
                                    time_fraction = max(0,min(1,((delta_z_acc*(box_acc(box_no_acc))-z_up)/speed_z)/h_step)); %the fraction of the time step (h) it takes from the previous time (i-1) to reach the grid box in question. Used to interpolate the particle speed in the particular grid boxes on the next 2 lines
                                    speed_x_add_acc(k,box_acc(box_no_acc)) = v_x(k,i-1)+time_fraction*(v_x(k,i-1)-v_x(k,i-2)); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box on the coarse grid that the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                    speed_z_add_acc(k,box_acc(box_no_acc)) = v_z(k,i-1)+time_fraction*(v_z(k,i-1)-v_z(k,i-1)); %speed_z_add_acc holds the list of interpolated vertical particle speeds in each grid box on the fine grid that the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                                else %the particle speed can not be extrapolated from v_x(k,i-2_ because a collision occurred at that time
                                    speed_x_add_acc(k,box_acc(box_no_acc)) = v_x(k,i-1); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box on the coarse grid that the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                    speed_z_add_acc(k,box_acc(box_no_acc)) = v_z(k,i-1); %speed_z_add_acc holds the list of interpolated vertical particle speeds in each grid box on the fine grid that the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                                end %if, extrapolating from v_x(k,i-2) and v_x(k,i-1), because v_x(k,i) describes the velocity upon rebound
                            elseif (box_no_acc == 1) %in this case box_no==1, but the particle speed must still be indexed because the particle just underwent a change in the direction of its velocity, meaning that its speed in the same box in the previous time step was indexed as an upward-moving speed
                                speed_x_add_acc(k,box_acc(box_no_acc)) = v_x(k,i-1); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box on the coarse grid that the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                speed_z_add_acc(k,box_acc(box_no_acc)) = v_z(k,i-1); %speed_z_add_acc holds the list of interpolated vertical particle speeds in each grid box on the fine grid that the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                            else  %this is the 'normal' case when nothing unusual happened (i.e., no collisions or speed changes), and the speed of the particle in each grid box is interpolated using its position in the previous and current time steps
                                time_fraction = max(0,min(1,((delta_z_acc*(box_acc(box_no_acc))-z_up)/speed_z)/h_step)); %the fraction of the time step (h) it takes from the previous time (i-1) to reach the grid box in question. Used to interpolate the particle speed in the particular grid boxes on the next 2 lines
                                speed_x_add_acc(k,box_acc(box_no_acc)) = v_x(k,i-1)+time_fraction*(v_x(k,i)-v_x(k,i-1)); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box on the coarse grid that the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                speed_z_add_acc(k,box_acc(box_no_acc)) = v_z(k,i-1)+time_fraction*(v_z(k,i)-v_z(k,i-1)); %speed_z_add_acc holds the list of interpolated vertical particle speeds in each grid box on the fine grid that the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                            end %if, assigning speed_x_add_acc and speed_z_add_acc
                            flux_counter_down_acc(k,box_acc(box_no_acc)) = flux_counter_down_acc(k,box_acc(box_no_acc)) + 1; %adding the current particle to the count of the flux of downward moving particle passing bin j on the fine grid
                            v_x_avg_down_acc(k,box_acc(box_no_acc)) = v_x_avg_down_acc(k,box_acc(box_no_acc)) + speed_x_add_acc(k,box_acc(box_no_acc)); %adding the horizontal speed of the current downward-moving particle to the sum of the horizontal downward-moving particle speeds in this box on the fine grid
                            v_z_avg_down_acc(k,box_acc(box_no_acc)) = v_z_avg_down_acc(k,box_acc(box_no_acc)) + speed_z_add_acc(k,box_acc(box_no_acc)); %adding the vertical speed of the current downward-moving particle to the sum of the vertical downward-moving particle speeds in this box on the fine grid
                        end %if, to prevent double-counting of particles with regards to the average particle speed
                        density_down_acc(k,box_acc(box_no_acc)) = density_down_acc(k,box_acc(box_no_acc)) + delta_time_acc(box_no_acc); %calculating the particle density based on the time the particle spends in each grid box
                    else %in this case the particle is moving upwards                
                        if (box_no_acc > 1 || collision == true) %to prevent double-counting of particles with regards to the average particle speed
                            if (box_no_acc == 1) %assigning speed_x_add_acc and speed_z_add_acc; in this case box_no==1, but the particle speed must still be indexed because the particle just underwent a change in the direction of its velocity, meaning that its speed in the same box in the previous time step was indexed as an upward-moving speed
                                speed_x_add_acc(k,box_acc(box_no_acc)) = v_x(k,i-1); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box on the coarse grid that the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                speed_z_add_acc(k,box_acc(box_no_acc)) = v_z(k,i-1); %speed_z_add_acc holds the list of interpolated vertical particle speeds in each grid box on the fine grid that the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                            else %this is the 'normal' case when nothing unusual happened (i.e., no collisions or speed changes), and the speed of the particle in each grid box is interpolated using its position in the previous and current time steps
                                time_fraction = max(0,min(1,((delta_z_acc*(box_acc(box_no_acc)-1)-z_down)/speed_z)/h_step)); %the fraction of the time step (h) it takes from the previous time (i-1) to reach the grid box in question. Used to interpolate the particle speed in the particular grid boxes on the next 2 lines
                                speed_x_add_acc(k,box_acc(box_no_acc)) = v_x(k,i-1)+time_fraction*(v_x(k,i)-v_x(k,i-1)); %speed_x_add holds the list of interpolated horizontal particle speeds in each grid box on the coarse grid that the particle traverses between the past and the current model time step and is used to calculate the particle shear stress and the average horizontal particle speed
                                speed_z_add_acc(k,box_acc(box_no_acc)) = v_z(k,i-1)+time_fraction*(v_z(k,i)-v_z(k,i-1)); %speed_z_add_acc holds the list of interpolated vertical particle speeds in each grid box on the fine grid that the particle traverses between the past and the current model time step and is used to calculate average vertical particle speed
                            end %if, assigning speed_x_add_acc and speed_z_add_acc;
                            flux_counter_up_acc(k,box_acc(box_no_acc)) = flux_counter_up_acc(k,box_acc(box_no_acc)) + 1; %adding the current particle to the count of the flux of upward moving particle passing bin j on the fine grid
                            v_x_avg_up_acc(k,box_acc(box_no_acc)) = v_x_avg_up_acc(k,box_acc(box_no_acc)) + speed_x_add_acc(k,box_acc(box_no_acc)); %adding the horizontal speed of the current upward-moving particle to the sum of the horizontal upward-moving particle speeds in this box on the fine grid
                            v_z_avg_up_acc(k,box_acc(box_no_acc)) = v_z_avg_up_acc(k,box_acc(box_no_acc)) + speed_z_add_acc(k,box_acc(box_no_acc)); %adding the vertical speed of the current upward-moving particle to the sum of the vertical upward-moving particle speeds in this box on the fine grid
                        end %to prevent double-counting                                     
                        density_up_acc(k,box_acc(box_no_acc)) = density_up_acc(k,box_acc(box_no_acc)) + delta_time_acc(box_no_acc); %calculating the particle density based on the time the particle spends in each grid box
                    end %if, setting density and average velocities
                end %for, cycling over the boxes that the particle was in between the last time step and the current one        
                if (abs(1-sum(delta_time_acc)/delta_time_acc_tot)>1e-6) %in this case the individual delta_times do not sum up to h_step, which is indicative of a problem.
                    warning('problem with box assigning loop');
                end %if, in this case the individual delta_times do not sum up to h_step, which is indicative of a problem.                
            end %if box_acc ~= 0, the particle spend time in the bottom two grid boxes
        end %if, data points are skipped in the rare instance when t(k,i-1) == t(k,i)
    end %for, cycling over each particle's orbit
    
    t_down = density_down; %saving the time spent in each grid box on the coarse grid for downward-moving particles
    t_up = density_up; %saving the time spent in each grid box on the coarse grid for upward-moving particles
    t_down_acc = density_down_acc; %saving the time spent in each grid box on the fine grid for downward-moving particles
    t_up_acc = density_up_acc; %saving the time spent in each grid box on the fine grid for upward-moving particles    
    total_t_down(k) = sum(t_down(k,:)); %the total time the particle spends moving downward; total_t_down(k) + total_t_up(k) = time_end(k)
    total_t_up(k) = sum(t_up(k,:)); %the total time the particle spends moving downward; total_t_down(k) + total_t_up(k) = time_end(k)
    density_down(k,1:1:j_end)=density_down(k,1:1:j_end)/(total_t_down(k)+total_t_up(k)); %normalizing the particle density on the coarse grid by the total time indexed such that sum(density_down(k,:))+sum(density_up(k,:)) = 1
    density_up(k,1:1:j_end)=density_up(k,1:1:j_end)/(total_t_down(k)+total_t_up(k)); %normalizing the particle density on the coarse grid by the total time indexed such that sum(density_down(k,:))+sum(density_up(k,:)) = 1
    density_down_acc(k,1:1:j_acc_end)=density_down_acc(k,1:1:j_acc_end)/(total_t_down(k)+total_t_up(k)); %normalizing the particle density on the fine grid by the total time indexed.
    density_up_acc(k,1:1:j_acc_end)=density_up_acc(k,1:1:j_acc_end)/(total_t_down(k)+total_t_up(k)); %normalizing the particle density on the fine grid by the total time indexed.
    if (min(min(density_down))<0||min(min(density_up))<0) %checking that the density is not smaller than zero anywhere
        error('density smaller than zero!');
    end %if, checking that the density is not smaller than zero anywhere
    
    for j=1:round(max(z(k,1:i_end(k)))/delta_z+0.5)  %the average velocity is the sum of all velocities times the respective time they spend in the box at that height, divided by the summed time spend by all particles in the box
        if((t_down(k,j)+t_up(k,j))>0) %setting v_x_avg while preventing division by zero
            v_x_avg(k,j) = (v_x_avg_down(k,j)+v_x_avg_up(k,j))/(flux_counter_down(k,j)+flux_counter_up(k,j)); %obtaining the average horizontal particle speed on the coarse grid
        end %if, %setting v_x_avg while preventing division by zero
        if (flux_counter_down(k,j)>0) %to prevent division by zero
            v_x_avg_down(k,j)=v_x_avg_down(k,j)/flux_counter_down(k,j); %obtaining the average horizontal particle speed of downward-moving particles on the coarse grid by dividing by the total number of particles that were indexed
            v_z_avg_down(k,j)=v_z_avg_down(k,j)/flux_counter_down(k,j); %obtaining the average vertical particle speed of downward-moving particles on the coarse grid by dividing by the total number of particles that were indexed
        else %in this case there was no particle flux at j
            v_x_avg_down(k,j)=0; %the average horizontal particle speed of downward-moving particles on the coarse grid is set to zero if no particles reached to this j
            v_z_avg_down(k,j)=0; %the average vertical particle speed of downward-moving particles on the coarse grid is set to zero if no particles reached to this j
        end %if, to prevent division by zero
        if (flux_counter_up(k,j)>0) %to prevent division by zero
            v_x_avg_up(k,j)=v_x_avg_up(k,j)/flux_counter_up(k,j); %obtaining the average horizontal particle speed of upward-moving particles on the coarse grid by dividing by the total number of particles that were indexed
            v_z_avg_up(k,j)=v_z_avg_up(k,j)/flux_counter_up(k,j); %obtaining the average vertical particle speed of upward-moving particles on the coarse grid by dividing by the total number of particles that were indexed
        else %in this case there was no particle flux at j
            v_x_avg_up(k,j)=0; %the average horizontal particle speed of upward-moving particles on the coarse grid is set to zero if no particles reached to this j
            v_z_avg_up(k,j)=0; %the average vertical particle speed of upward-moving particles on the coarse grid is set to zero if no particles reached to this j
        end %if, to prevent division by zero
        mass_flux(k,j) = (density_up(k,j)*v_x_avg_up(k,j)+density_down(k,j)*v_x_avg_down(k,j))*m(k)*n(k); %calculating the horizontal particle mass flux from the particle density and average horizontal particle speed
    end %for, the average velocity/density is the sum of all velocities times the respective time they spend in the box at that height, divided by the summed time spend by all particles in the box
    
    for j_acc=1:1:j_acc_end  %the average velocity is the sum of all velocities times the respective time they spend in the box at that height, divided by the summed time spend by all particles in the box
        if((t_down_acc(k,j_acc)+t_up_acc(k,j_acc))>0) %setting v_x_avg while preventing division by zero
            v_x_avg_acc(k,j_acc) = (v_x_avg_down_acc(k,j_acc)+v_x_avg_up_acc(k,j_acc))/(flux_counter_down_acc(k,j_acc)+flux_counter_up_acc(k,j_acc)); %obtaining the average horizontal particle speed on the fine grid
        end %if, %setting v_x_avg while preventing division by zero        
        if (flux_counter_down_acc(k,j_acc)>0) %to prevent division by zero
            v_x_avg_down_acc(k,j_acc)=v_x_avg_down_acc(k,j_acc)/flux_counter_down_acc(k,j_acc); %obtaining the average horizontal particle speed of downward-moving particles on the fine grid by dividing by the total number of particles that were indexed
            v_z_avg_down_acc(k,j_acc)=v_z_avg_down_acc(k,j_acc)/flux_counter_down_acc(k,j_acc); %obtaining the average vertical particle speed of downward-moving particles on the fine grid by dividing by the total number of particles that were indexed
        else %in this case there was no particle flux at j_acc
            v_x_avg_down_acc(k,j_acc)=0; %the average horizontal particle speed of downward-moving particles on the fine grid is set to zero if no particles reached to this j_acc
            v_z_avg_down_acc(k,j_acc)=0; %the average vertical particle speed of downward-moving particles on the fine grid is set to zero if no particles reached to this j_acc
        end %if, to prevent division by zero
        if (flux_counter_up_acc(k,j_acc)>0) %to prevent division by zero
            v_x_avg_up_acc(k,j_acc)=v_x_avg_up_acc(k,j_acc)/flux_counter_up_acc(k,j_acc); %obtaining the average horizontal particle speed of upward-moving particles on the fine grid by dividing by the total number of particles that were indexed
            v_z_avg_up_acc(k,j_acc)=v_z_avg_up_acc(k,j_acc)/flux_counter_up_acc(k,j_acc); %obtaining the average vertical particle speed of upward-moving particles on the fine grid by dividing by the total number of particles that were indexed
        else %in this case there was no particle flux at j_acc
            v_x_avg_up_acc(k,j_acc)=0; %the average horizontal particle speed of upward-moving particles on the fine grid is set to zero if no particles reached to this j_acc
            v_z_avg_up_acc(k,j_acc)=0; %the average vertical particle speed of upward-moving particles on the fine grid is set to zero if no particles reached to this j_acc
        end %if, to prevent division by zero        
    end %for, the average velocity/density is the sum of all velocities times the respective time they spend in the box at that height, divided by the summed time spend by all particles in the box
    
    if (time_end(k) > 0) %normalizing the particle fluxes by the simulation time
        flux_up(k,1:1:size(t_up,2)) = flux_counter_up(k,1:1:size(flux_counter_up,2))/time_end(k); %normalizing the particle flux by the simulation time to obtain the particle flux per second. This is used in calculating the particle shear stress
        flux_down(k,1:1:size(t_down,2)) = flux_counter_down(k,1:1:size(flux_counter_down,2))/time_end(k); %normalizing the particle flux by the simulation time to obtain the particle flux per second. This is used in calculating the particle shear stress
        flux_up_acc(k,1:1:size(t_up_acc,2)) = flux_counter_up_acc(k,1:1:size(flux_counter_up_acc,2))/time_end(k); %normalizing the particle flux by the simulation time to obtain the particle flux per second. This is used in calculating the particle shear stress
        flux_down_acc(k,1:1:size(t_down_acc,2)) = flux_counter_down_acc(k,1:1:size(flux_counter_down_acc,2))/time_end(k); %normalizing the particle flux by the simulation time to obtain the particle flux per second. This is used in calculating the particle shear stress
    else %this is to prevent division by zero when time_end(k)==0, which should never occur
        warning('time_end is zero!'); %warning that time_end(k)==0, which should never occur
        flux_up(k,1:1:size(t_up,2)) = zeros(1,size(flux_counter_up,2)); %setting the fluxes to zero if the simulation was not performed
        flux_down(k,1:1:size(t_down,2)) = zeros(1,size(flux_counter_down,2));
        flux_up_acc(k,1:1:size(t_up_acc,2)) = zeros(1,size(flux_counter_up_acc,2));
        flux_down_acc(k,1:1:size(t_down_acc,2)) = zeros(1,size(flux_counter_down_acc,2));        
    end %if, normalizing the particle fluxes by the simulation time
end %for, cycling over all particle bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following calculates the average particle flux, speed, and mass flux %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
del = floor(r/delta_z+1); %del is the number of coarse grid boxes spanned by the particle's radius
del_acc = floor(r/delta_z_acc+1); %del_acc is the number of fine grid boxes spanned by the particle's radius
for k=1:part_no %this is to account for the fact that particles don't rebound at z = 0, but at z = r.  On scales smaller than r, the shear stress is thus not realistically represented by a 'point-source' particle.  So the shear stress at z < r is set to the shear stress at z = r
    flux_down(k,1:del(k)) = flux_down(k,1+del(k)); flux_down_acc(k,1:del_acc(k)) = flux_down_acc(k,1+del_acc(k)); %correcting flux_down_acc and flux_down
    flux_up(k,1:del(k)) = flux_up(k,1+del(k)); flux_up_acc(k,1:del_acc(k)) = flux_up_acc(k,1+del_acc(k)); %correcting flux_up_acc and flux_up
    v_x_avg_up_acc(k,1:del_acc(k)) = v_x_avg_up_acc(k,1+del_acc(k)); v_x_avg_up(k,1:del(k)) = v_x_avg_up(k,1+del(k)); %correcting v_x_avg_up_acc and v_x_avg_up
    v_x_avg_down_acc(k,1:del_acc(k)) = v_x_avg_down_acc(k,1+del_acc(k)); v_x_avg_down(k,1:del(k)) = v_x_avg_down(k,1+del(k)); %correcting v_x_avg_down_acc and v_x_avg_down
    v_z_avg_up_acc(k,1:del_acc(k)) = v_z_avg_up_acc(k,1+del_acc(k)); v_z_avg_up(k,1:del(k)) = v_z_avg_up(k,1+del(k)); %correcting v_z_avg_up_acc and v_z_avg_up
    v_z_avg_down_acc(k,1:del_acc(k)) = v_z_avg_down_acc(k,1+del_acc(k)); v_z_avg_down(k,1:del(k)) = v_z_avg_down(k,1+del(k)); %correcting v_z_avg_down_acc and v_z_avg_down
    t_up_acc(k,1:del_acc(k)) = t_up_acc(k,1+del_acc(k)); t_up(k,1:del(k)) = t_up(k,1+del(k)); %correcting t_up_acc and t_up
    t_down_acc(k,1:del_acc(k)) = t_down_acc(k,1+del_acc(k)); t_down(k,1:del(k)) = t_down(k,1+del(k)); %correcting t_down_acc and t_down
end %for, this is to account for the fact that particles don't rebound at z = 0, but at z = r.  On scales smaller than r, the shear stress is thus not realistically represented by a 'point-source' particle.  So the shear stress at z < r is set to the shear stress at z = r
total_mass_flux = sum(sum(mass_flux)); %this calculates the total height-integrated horizontal mass flux of all particle bins
for k=1:1:part_no %cycling over the particle bins
    total_v_x_avg(k) = v_x_avg(k,1:size(v_x_avg,2))*(density_down(k,1:size(v_x_avg,2))+density_up(k,1:size(v_x_avg,2)))'; %finding the final average horizontal particule speed on the coarse grid by normalizing by the particle density
    total_v_x_avg_acc(k) = v_x_avg_acc(k,1:size(v_x_avg_acc,2))*(density_down_acc(k,1:size(v_x_avg_acc,2))+density_up_acc(k,1:size(v_x_avg_acc,2)))'; %finding the final average horizontal particule speed on the fine grid by normalizing by the particle density
end %for, cycling over the particle bins