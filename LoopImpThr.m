function LoopImpThr (t_max, psd_flag, filetext, spare)

%Loop_ufrthrimp uses findImpThr to determine the impact threshold of any size
%distribution under Earth ambient condition, as specified in
%load_parameters. The array of flags for particular size distributions
%should be specified in psd_flag. The result is written to a data file.

fidData = fopen(strcat('LoopImpThr',filetext,'.txt'),'wt'); %opening the data file
fprintf(fidData, 'd  ImpThr       z50\n'); %writing the header for the data file

for i = 1:1:size(psd_flag,2) %cycling over the array of size distributions
    try %preventing potentials errors in findImpThr and saltation_model to crash the model
        if (psd_flag(i)==1) %setting an initial impact threshold; for the size distribution of Namikas (2003), the impact threshold should lie around 0.21 m/s
            u_fr_thr(i) = 0.21;
        elseif (psd_flag(i)==2) %for the size distribution of Bagnold (1938), the impact threshold should lie around 0.22 m/s
            u_fr_thr(i) = 0.22;
        elseif (psd_flag<1) %in this case its a monodisperse size distribution, and the approximate value of that is obtained using Bagnold's empirical relation
            u_fr_thr(i) = (0.47/5.75)*sqrt(((2600-1.2)/1.2)*9.8.*psd_flag(i));
        else %in this case, the size distribution is not monodisperse and not those of Namikas (2003) and Bagnold (1938)
            u_fr_thr(i) = 0.21;
        end %if, setting an initial impact threshold
        [u_fr_thr_imp(i),z50(i)] = findImpThr (t_max, psd_flag(i), u_fr_thr(i), spare); %calling findImpThr to determine the actual impact threshold
        fprintf(fidData, '%1.6f    %1.5f %1.5f\n', psd_flag(i), u_fr_thr_imp(i), z50(i)); %writing the result of the call to findImpThr to a data file
    catch
        warning('error in findImpThr');
    end %try, preventing potentials errors in findImpThr and saltation_model to crash the model
end %for, cycling over the array of size distributions
fclose(fidData); %closing the data file