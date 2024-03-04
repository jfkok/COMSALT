function [total_part_no, total_d, d_mean, total_mass_fraction] = PSD_Namikas (d_test)

%the mass fraction for each size bin is obtained by integrating the
%relevant Gaussian PSD from the specified phi values (Dp in mm = 2^(-phi)
%). The mean particle size is obtained as D =
%int(G*d^,d_min..d_max)/int(G,d_min...d_max), where G is the function
%(usually a Gaussian) describing the particle size distribution

%the size distribution from Namikas (2003), with phi = 2.00 +/- 0.37 (5 bins)
%1:  phi = -inf - 1.4;   D = 379 - inf um; mean D = 424 um; mass freq. = 5.244 %
%2:  phi = 1.4 - 1.8; D = 287 - 379 um; mean D = 323 um; mass freq. = 24.197 %
%3:  phi = 1.8 - 2.2; D = 218 - 287 um; mean D = 251 um; mass freq. = 41.118 %
%4:  phi = 2.2 - 2.6; D = 165 - 218 um; mean D = 195 um; mass freq. = 24.197 %
%5:  phi = 2.6 - inf;   D = -inf - 165 um; mean D = 149 um; mass freq. = 5.244 %

%the size distribution from Namikas (2003), with phi = 2.00 +/- 0.37 (10 bins)
%1:  phi = -inf - 1.00; D = 500 - inf um; mean D = 542 um; mass freq. = 0.344 %
%2:  phi = 1.00 - 1.25; D = 420 - 500 um; mean D = 449 um; mass freq. = 1.789 %
%3:  phi = 1.25 - 1.50; D = 354 - 420 um; mean D = 380 um; mass freq. = 6.696 %
%4:  phi = 1.50 - 1.75; D = 297 - 354 um; mean D = 321 um; mass freq. = 16.133 %
%5:  phi = 1.75 - 2.00; D = 250 - 297 um; mean D = 272 um; mass freq. = 25.038 %
%6:  phi = 2.00 - 2.25; D = 210 - 250 um; mean D = 230 um; mass freq. = 25.038 %
%7:  phi = 2.25 - 2.50; D = 177 - 210 um; mean D = 195 um; mass freq. = 16.133 %
%8:  phi = 2.50 - 2.75; D = 149 - 177 um; mean D = 165 um; mass freq. = 6.696 %
%9:  phi = 2.75 - 3.00; D = 125 - 149 um; mean D = 139 um; mass freq. = 1.789 %
%10: phi = 3.00 - inf; D = -inf - 125 um; mean D = 116 um; mass freq. = 0.344 %

d_mean = 0.000250; % Median diameter ~250 um
bin5 = false; %if this is set to true, then the soil size distribution will be divided into 5 bins
bin10 = true; %if this is set to true, then the soil size distribution will be divided into 10 bins

if (bin5==true)
    total_d(1) = 0.000149; %the size of the particle populations in meters
    total_d(2) = 0.000195; 
    total_d(3) = 0.000251; 
    total_d(4) = 0.000323;
    total_d(5) = 0.000424;

    total_mass_fraction(1) = 0.0507; %the volume fraction of the soil occupied by a particular particle bin
    total_mass_fraction(2) = 0.2418;
    total_mass_fraction(3) = 0.4150;
    total_mass_fraction(4) = 0.2418;
    total_mass_fraction(5) = 0.0507;

    total_part_no = size(total_d,2) %the total number of particle bins
end %bin 5

if (bin10==true)
    total_d(1) = 0.000116; %the size of the particle populations in meters
    total_d(2) = 0.000139; 
    total_d(3) = 0.000165; 
    total_d(4) = 0.000195;
    total_d(5) = 0.000230;
    total_d(6) = 0.000272; 
    total_d(7) = 0.000321; 
    total_d(8) = 0.000380; 
    total_d(9) = 0.000449;
    total_d(10)= 0.000542;

    total_mass_fraction(1) = 0.00344; %the volume fraction of the soil occupied by a particular particle bin
    total_mass_fraction(2) = 0.01789;
    total_mass_fraction(3) = 0.06696;
    total_mass_fraction(4) = 0.16133;
    total_mass_fraction(5) = 0.25038;
    total_mass_fraction(6) = 0.25038;
    total_mass_fraction(7) = 0.16133;
    total_mass_fraction(8) = 0.06696;
    total_mass_fraction(9) = 0.01789;
    total_mass_fraction(10)= 0.00344;

    total_part_no = size(total_d,2) %the total number of particle bins
end %bin10

if (max(d_test) ~= 0) %this creates a test particle of a particular size. The model calculates the mass flux profile, saltation trajectories, and other relevant parameters for this test particle. Note that the test does not affect the solution, since its concentration is infinitesimal
    total_d(total_part_no+1:total_part_no+size(d_test,2)) = d_test; %appending the test particle size(s) to total_d
    total_mass_fraction = total_mass_fraction*0.99999999;
    total_mass_fraction(total_part_no+1:total_part_no+size(d_test,2)) = 0.00000001/size(d_test,2);
    total_part_no = size(total_d,2);
end %if, this creates a test particle