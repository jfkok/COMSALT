function [Uavg, dU] = calc_windspeed (u, u_acc, d, z, j, j_end_prev, acc_grid)

%these are global parameters
global j_acc_end delta_z delta_z_acc u_fr pi;

z_points = 3; %number of points used on the (spherical) particle to use in integrate the wind speed experienced by the particle
Uavg = 0; %initiating Uavg

if (z<d/2) %checking if the particle is below the surface.  If so, it will recalculated in saltation_trajectory anyways
    Uavg = 0; %setting Uavg to zero below the surface
    dU = 10^(-10); %the small number is there to prevent division by zero
else %normal wind calculation if the particle is not colliding with the surface
    for p=1:1:z_points %looping over the number of points on the spherical particle for which the wind speed is calculated
        if (z_points>1) %determining the points and angles on the particle for which the wind speed is calculated
            z_iter(p) = z - d/2 + (p-1)*d/(z_points-1); %calculating the vertical position of the z_points on the sphere for which the wind speed is calculated
            theta_iter(p) = acos(1-(2*p-1)/(z_points-1)); %calculating the angle of the half-way point between two points for which the wind speed is calculated
        else %in this case, z_iter and theta_iter are trivial
            z_iter(p) = z;
            theta_iter(p) = pi;
        end %if, determining the points and angles on the particle for which the wind speed is calculated
        if (p==1) %calculating the area of the sphere to which each wind speed calculation is applicable
            area_iter(p) = (d/2)^2*(theta_iter(p)-0.5*sin(2*theta_iter(p)));
        elseif (p==z_points)
            area_iter(p) = (d/2)^2*(pi-theta_iter(p-1)+0.5*sin(2*theta_iter(p-1)));
        else
            area_iter(p) = (d/2)^2*(theta_iter(p)-theta_iter(p-1)+0.5*(sin(2*theta_iter(p-1))-sin(2*theta_iter(p))));
        end %if, calculating the area of the sphere to which each wind speed calculation is applicable
        if (acc_grid == true) %checking if the particle is on the fine grid or on the course grid
            j = round(z_iter(p)/delta_z_acc+0.5); %calculating the position on the fine grid
            if (j==1) %in this case the particle is near the surface
                U0 = u_acc(j);  %setting the parameters for the Lagrangian interpolation of the wind
                U1 = u_acc(j+1);
                U2 = u_acc(j+2);
                Z0 = j*delta_z_acc;            
                Z1 = (j+1)*delta_z_acc;
                Z2 = (j+2)*delta_z_acc;                            
            elseif (j>=j_acc_end) %in this case the particle is near the top of the fine grid
                U0 = u_acc(j_acc_end-2); %setting the parameters for the Lagrangian interpolation of the wind
                U1 = u_acc(j_acc_end-1);
                U2 = u_acc(j_acc_end);
                Z0 = (j_acc_end-2)*delta_z_acc;            
                Z1 = (j_acc_end-1)*delta_z_acc;
                Z2 = j_acc_end*delta_z_acc; 
            else %in this case the particle is not near the surface or the top of the fine grid
                U0 = u_acc(j-1); %setting the parameters for the Lagrangian interpolation of the wind
                U1 = u_acc(j);
                U2 = u_acc(j+1);
                Z0 = (j-1)*delta_z_acc;            
                Z1 = j*delta_z_acc;
                Z2 = (j+1)*delta_z_acc;                            
            end %if, setting the parameters for the Lagrangian interpolation of the wind
            U_iter(p) = max(0,(U0*(z_iter(p)-Z1)*(z_iter(p)-Z2)/((Z0-Z1)*(Z0-Z2)) +  U1*(z_iter(p)-Z0)*(z_iter(p)-Z2)/((Z1-Z0)*(Z1-Z2)) + U2*(z_iter(p)-Z0)*(z_iter(p)-Z1)/((Z2-Z0)*(Z2-Z1)))); %this is the Lagrangian interpolation of the wind with the parameters calculated above
        else %in this case the particle is above the fine grid and is treated by the coarse grid
            if (j<j_end_prev-1) %seeing if the particle is below j_end*delta_z
                U0 = u(j-1); %setting the parameters for the Lagrangian interpolation of the wind
                U1 = u(j);
                U2 = u(j+1);
                Z0 = (j-1)*delta_z;
                Z1 = (j)*delta_z;
                Z2 = (j+1)*delta_z;                
                U_iter(p) = max(0,(U0*(z_iter(p)-Z1)*(z_iter(p)-Z2)/((Z0-Z1)*(Z0-Z2)) +  U1*(z_iter(p)-Z0)*(z_iter(p)-Z2)/((Z1-Z0)*(Z1-Z2)) + U2*(z_iter(p)-Z0)*(z_iter(p)-Z1)/((Z2-Z0)*(Z2-Z1)))); %this is the Lagrangian interpolation of the wind with the parameters calculated above
            else %extending the wind profile to above the coarse grid, for the case that the particle is outside of the coarse grid
                U_iter(p) = max(0,u(j_end_prev) + (u_fr/0.4)*(log(z_iter(p)/(j_end_prev*delta_z)))); %this is the Lagrangian interpolation of the wind with the parameters calculated above
            end %if, seeing if the particle is within the coarse grid
        end %if, %checking if the particle is on the fine grid or on the course grid
        Uavg = Uavg + area_iter(p)*U_iter(p)/(pi*d^2/4);  %the average wind speed on the sphere's surface, weighted by cross-sectional area
    end %for, looping over the number of points on the spherical particle for which the wind speed is calculated
    dU = max(10^(-10),(U_iter(z_points)-U_iter(1))/d); %the wind gradient is defined as delta(wind)/delta(z); the small number is there to prevent divition by zero
end %if, checking if the particle is below the surface.  If so, it will recalculated in saltation_trajectory anyways