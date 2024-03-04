function [total_part_no, total_d, d_mean, total_mass_fraction] = PSD_Dune

%the size distribution from Lancaster, Journal of Sedimentary Petrology, 1986
%the dune crest with phi = 2.16 +/- 0.49
%1:  phi = 0.91 - 1.41; D = 376 - 532 um; mean D = 448 um; freq. = 6.09 %
%2:  phi = 1.41 - 1.91; D = 266 - 376 um; mean D = 316 um; freq. = 24.94 %
%3:  phi = 1.91 - 2.41; D = 188 - 266 um; mean D = 224 um; freq. = 37.94 %
%4:  phi = 2.41 - 2.91; D = 133 - 188 um; mean D = 158 um; freq. = 24.94 %
%5:  phi = 2.91 - 3.41; D = 94 - 133 um; mean D = 112 um; freq. = 6.09 %

d_mean = 0.000224; % Median diameter ~224 um

total_d(1) = 0.000112; %the size of the particle populations in meters
total_d(2) = 0.000158; 
total_d(3) = 0.000224; 
total_d(4) = 0.000316;
total_d(5) = 0.000448;

total_mass_fraction(1) = 0.0609; %the volume fraction of the soil occupied by a particular particle bin
total_mass_fraction(2) = 0.2494;
total_mass_fraction(3) = 0.3794;
total_mass_fraction(4) = 0.2494;
total_mass_fraction(5) = 0.0609;

total_part_no = size(total_d,2) %the total number of particle bins