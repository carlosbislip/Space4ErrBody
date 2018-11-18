function [evolutions,output] = Load_Data(option,p,Folder_Path_List,pop_path,npop_files)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



%% Load lists of analyzed simulations

switch option
    
    case 1
        
        % Figure out output name based on what is in folder path
        if contains(convertCharsToStrings(Folder_Path_List{p}),'OUTPUT') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'T_'),'.mat'));
        elseif contains(convertCharsToStrings(Folder_Path_List{p}),'VALIDATION') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'N_'),'.mat'));
        end
        
        % Open file if it exists. Create variable otherwise.
        if exist(output,'file') == 2
            load(output)
        else
            disp('Run again with option 2 or 3. No .mat file found.');
            
            %             k = npop_files + 1;
            %             evolutions(k).evolution(npop_files)    = nan;
            %             evolutions(k).v_i(npop_files)          = nan;
            %             evolutions(k).gamma_i(npop_files)      = nan;
            %             evolutions(k).chi_i(npop_files)        = nan;
            %             evolutions(k).trajectories(npop_files) = nan;
            %             evolutions(k).fitness(npop_files)      = nan;
            %             evolutions(k).lat(npop_files)          = nan;
            %             evolutions(k).lon(npop_files)          = nan;
            %             evolutions(k).tof(npop_files)          = nan;
        end
        
        idx = [];
        
    case {2, 3}
        
        % Figure out output name based on what is in folder path
        if contains(convertCharsToStrings(Folder_Path_List{p}),'OUTPUT') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'T_'),'.mat'));
        elseif contains(convertCharsToStrings(Folder_Path_List{p}),'VALIDATION') == 1
            output = char(strcat('evolutions_',extractAfter(Folder_Path_List{p},'N_'),'.mat'));
        end
        
        % Open file if it exists. Create variable otherwise.
        if exist(output,'file') == 2
            load(output)
        else
            k = npop_files;
            evolutions(k).evolution       = nan;
            evolutions(k).trajectories.individual.t   = nan;
            evolutions(k).trajectories.individual.time_vector = nan;
            evolutions(k).trajectories.individual.x_R = nan;
            evolutions(k).trajectories.individual.y_R = nan;
            evolutions(k).trajectories.individual.z_R = nan;
            evolutions(k).trajectories.individual.R_R_norm   = nan;
            evolutions(k).trajectories.individual.lat = nan;
            evolutions(k).trajectories.individual.lon = nan;
            evolutions(k).trajectories.individual.h   = nan;
            evolutions(k).trajectories.individual.latitude_angle    = nan;
            evolutions(k).trajectories.individual.longitude_angle   = nan;
            evolutions(k).trajectories.individual.heading_angle     = nan;
            evolutions(k).trajectories.individual.heading_required  = nan;
            evolutions(k).trajectories.individual.heading_error     = nan;
            evolutions(k).trajectories.individual.flight_path_angle = nan;
            evolutions(k).trajectories.individual.angle_of_attack   = nan;
            evolutions(k).trajectories.individual.angle_of_sideslip = nan;
            evolutions(k).trajectories.individual.bank_angle        = nan;
            evolutions(k).trajectories.individual.d_deg             = nan;
            evolutions(k).trajectories.individual.airspeed          = nan;
            evolutions(k).trajectories.individual.E                 = nan;
            evolutions(k).trajectories.individual.interp_E                 = nan;
            evolutions(k).trajectories.individual.interp_angle_of_attack   = nan;
            evolutions(k).trajectories.individual.interp_thrust_angle      = nan;
            evolutions(k).trajectories.individual.interp_throttle_setting  = nan;
%            evolutions(k).trajectories.individual.angle_of_attack   = nan;
            evolutions(k).individuals.v_i     = nan;
            evolutions(k).individuals.gamma_i = nan;
            evolutions(k).individuals.chi_i   = nan;
            evolutions(k).fitness.dif_norm  = nan;
            evolutions(k).fitness.dif_lat   = nan;
            evolutions(k).fitness.dif_lon   = nan;
            evolutions(k).fitness.dif_d_deg = nan;
            evolutions(k).fitness.dif_h     = nan;
            evolutions(k).fitness.tof       = nan;
            evolutions(k).best(3).criteria  = nan;
            evolutions(k).best(3).v_i       = nan;
            evolutions(k).best(3).gamma_i   = nan;
            evolutions(k).best(3).chi_i     = nan;
            evolutions(k).best(3).dif_d_deg = nan;
            evolutions(k).best(3).dif_h     = nan;
            evolutions(k).best(3).tof       = nan;
            
        end
        
        %         idx = nan(1,npop_files);
        %         % Anaysis loop
        %         parfor k = 1:npop_files
        %
        %             % Compare each string to the lists of analyzed simulations. There must
        %             % be a way of comparing entire arrays in one go without looping. Then
        %             % just analyzing the remaining indices. This woudl be done before
        %             % entering this loop. Maybe use:    ~any(X)
        %             discriminant = find(strcmp(pop_path(k,:),{analyzed_simulations.path}) == 1);
        %
        %             if isempty(discriminant) == 1
        %                 idx(1,k) = k;
        %             end
        %
        %         end
        %
        %         idx = idx(~isnan(idx));
        
        
end

