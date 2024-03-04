function [total_part_no, total_d, d_mean, total_mass_fraction] = PSD_Williams

%the size distribution used in Williams, Sedimentology, 1964

d_mean = 0.000306; % Median diameter ~306 um

total_d(1) = 0.000161; %the size of the particle populations in meters
total_d(2) = 0.000215; 
total_d(3) = 0.000297; 
total_d(4) = 0.000414;
total_d(5) = 0.000550;

total_mass_fraction(1) = 0.05285; %the volume fraction of the soil occupied by a particular particle bin
total_mass_fraction(2) = 0.24073;
total_mass_fraction(3) = 0.36011;
total_mass_fraction(4) = 0.2418;
total_mass_fraction(5) = 0.0507;

total_part_no = size(total_d,2) %the total number of particle bins