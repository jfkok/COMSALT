function [total_part_no, total_d, d_mean, total_mass_fraction] = PSD_Bagnold

%the size distribution described in Bagnold, Proceedings of the Royal
%Society of London Series A - Mathematical and Physical Sciences, 1938

d_mean = 0.000282; % Median dameter ~282 um

total_d(1) = 0.000161; %the size of the particle populations in meters
total_d(2) = 0.0002225; 
total_d(3) = 0.000282; 
total_d(4) = 0.000360;
total_d(5) = 0.000506;

total_mass_fraction(1) = 0.02025; %the volume fraction of the soil occupied by a particular particle bin
total_mass_fraction(2) = 0.20042;
total_mass_fraction(3) = 0.48445;
total_mass_fraction(4) = 0.27356;
total_mass_fraction(5) = 0.02132;

total_part_no = size(total_d,2) %the total number of particle bins