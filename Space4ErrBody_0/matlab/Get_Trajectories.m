function [ evolutions ] = Get_Trajectories(evolutions,prop_path,depvar_path,interp_path,v_i,gamma_i,pop_i,lon_i_rad,lat_f_deg,lon_f_deg, startEpoch)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

population =  numel(prop_path)/(numel(pop_i) + 1);

a = linspace(population,numel(prop_path),numel(pop_i) + 1);

% population =  100;
% pop_i = [ 0 1];
% a = linspace(population,200,numel(pop_i) + 1);

start = [1;a(1:end-1)'+1];
%C = [b a'];
%start = b;
finish = a';
%evolutions(numel(pop_i)).evolution(population,1).individual = [];

for k = 1:(numel(pop_i) + 1)
    
    evolutions(k).evolution = k;
    
    pp = 1;
    %     prop_source = prop_path(start(k):finish(k),:);
    depvar_source = depvar_path(start(k):finish(k),:);
    interp_source = interp_path(start(k):finish(k),:);
    for p = 1:population
        
        % Open this file and read in the contents
        %         fid = fopen(prop_source{p,:});
        %         prop = dlmread(prop_source{p,:},',');
        %         fclose(fid);
        fid = fopen(depvar_source{p,:});
        depvar = dlmread(depvar_source{p,:},',');
        fclose(fid);
        
        fid = fopen(interp_source{p,:});
        interp = dlmread(interp_source{p,:},',');
        fclose(fid);
        
        % analyzed_simulations(k).tof = data(1,1) - data(end,1);
        %data = [2458419.95833333, -1730261.1974541, 3489942.65996152, 5041452.57112492, 1817.13843876246, 7053.76703722487, 474.632964784125]
        %
        %
        %         t   = prop(:,1);
        %         x_I = prop(:,2);
        %         y_I = prop(:,3);
        %         z_I = prop(:,4);
        %         R_I_vect = [x_I y_I z_I];
        %         %R_I_norm = sqrt(x_I.^2+y_I.^2+z_I.^2);
        %         %     %V_I_norm = sqrt(u_I.^2+v_I.^2+w_I.^2);
        %         tau_I_rad = atan2(y_I,x_I);
        %         %     I = find(tau_I_rad<0);
        %         %     tau_I_rad(I) = 2*pi + tau_I_rad(I);
        %         %     tau_I_rad2 = tau_I_rad;
        %         %     %tau_I_rad = [atan2(data(1,3),data(1,2)) atan2(data(end,3),data(end,2))];
        %         %     tau_I_deg = tau_I_rad*180/pi;
        %         %     tau_I_deg2 = tau_I_deg;
        %         %     delta_I_rad = asin(z_I./R_I_norm);
        %         %     I = find(delta_I_rad<0);
        %         %     delta_I_rad(I) = 2*pi + delta_I_rad(I);
        %         %     delta_I_deg = delta_I_rad*180/pi;
        %         %     delta_I_rad2 = delta_I_rad;
        %         %     delta_I_deg2 = delta_I_deg;
        %
        %
        %         %%%%%%%% This section beloy may be mitigated by transforming with TUDAT
        %         %%%%%%%% first and then exporting/printing.
        %         %%%%%%%%
        %         % Adjust Inertial Frame for plotting. TUDAT does all calculations
        %         % in the Inertial Frame, which is a function of the starting epoch.
        %         % TUDAT uses Earth's rotational ephemeris to identify the
        %         % corresponding state of the initial conditions. This creates a
        %         % basic rotational offset that can be easily reversed by
        %         % calculating the longidutindal difference. This difference is then
        %         % used with a transformation matrix that is then applied to the
        %         % entire trajectory.
        %
        %         % Calculate Transformation Matrix from original Inertial Frame to
        %         % Earth-Fixed Frame
        %
        %         initial_longitude_offset = lon_i_rad - tau_I_rad(1);
        %         st_adjust = [ zeros(numel(t),1) zeros(numel(t),1) -initial_longitude_offset*ones(numel(t),1) ];
        %         C_I_to_EF = DCM(st_adjust); % 3x3xN
        %
        %         % Restructure Inertial Frame position vector to 3x1xN
        %         [r,c] = size(R_I_vect);
        %         pages = numel(t);
        %         R_I_vect = permute(reshape(permute(R_I_vect,[2,1]),[r/pages,c,pages]),[2,1,3]);
        %
        %         % Calculate position vector and related vectors/values in Earth-Fixed Frame
        %         R_EF_vect = mmat(C_I_to_EF,R_I_vect,[1 2]);
        %         R_EF_vect = ipermute(squeeze(R_EF_vect),[2,1]);
        
        % x_EF(:,1) = squeeze(R_EF_vect(1,1,:));
        % y_EF(:,1) = squeeze(R_EF_vect(2,1,:));
        % z_EF(:,1) = squeeze(R_EF_vect(3,1,:));
        % R_EF_norm = sqrt(x_EF.^2+y_EF.^2+z_EF.^2);
        %
        %     %clear R_I_vect;
        %     R_I_vect = [x_I y_I z_I];
        %
        %     % Calcualte adjusted longitude and latitude
        %     delta_I_rad = asin(z_I./R_I_norm);
        %     delta_I_deg = delta_I_rad*180/pi;
        %     tau_I_rad = atan2(y_I,x_I);
        %     %tau_I_rad = [atan2(data(1,3),data(1,2)) atan2(data(end,3),data(end,2))];
        %     tau_I_deg = tau_I_rad*180/pi;
        %
        %
        %
        %     % Extract coordinates of final state from simulation data
        %     lat_f_calc_deg = delta_I_deg(end);
        %     lon_f_calc_deg = tau_I_deg(end);
        %
        %
        %
        %     [lat_f_calc_deg,lon_f_calc_deg] = latlon_orientation(lat_f_calc_deg,lon_f_calc_deg)
        %
        %
        %     % Calculate final state offset from IAD coordinates
        %     dif_lat_deg = lat_IAD_deg - lat_f_calc_deg;
        %     dif_lon_deg = lon_IAD_deg - lon_f_calc_deg_pos;
        %
        %
        % Following the rationale from before, throughout the propagation Earth
        % is also rotating. However, plotting this in a fixed image is not
        % possible. Hence, to plot the trajectory it must then be converted to
        % the Rotating Planetocentric Frame. This transformation would consider
        % the translation of the point mass with respect to the rotating
        % surface of the Earth. Unlike the simple parabolic trajectory that is
        % observed in the Inertial frame, this frame would show some odd
        % trajectories. Anyhow, it is a function of the epoch of the state
        % relative to the starting epoch.
        
        %         % Calculate Transformation Matrix from Inertial Frame to Rotational Frame
        %         omega_E = 7.2921150e-5;
        %         del_t = t-t(1);
        %         st_RI = [ zeros(numel(t),1) zeros(numel(t),1) (omega_E*del_t) ];
        %         C_EF_to_R = DCM(st_RI); % 3x3xN
        %
        %         % Restructure Inertial Frame position vector to 3x1xN
        %         [r,c] = size(R_EF_vect);
        %         pages = numel(t);
        %         R_EF_vect = permute(reshape(permute(R_EF_vect,[2,1]),[r/pages,c,pages]),[2,1,3]);
        %
        %         % Calculate position vector and related vectors/values in Rotational Frame
        %         R_R_vect = mmat(C_EF_to_R,R_EF_vect,[1 2]);
        %         R_R_vect = ipermute(squeeze(R_R_vect),[2,1]);
        %         R_R_norm = sqrt(R_R_vect(:,1).^2 + R_R_vect(:,2).^2 + R_R_vect(:,3).^2 );
        
        %%%%%%%%
        
        %         evolutions(k).trajectories(pp).individual.t = t(1:end);
        %         evolutions(k).trajectories(pp).individual.x_R = R_R_vect(1:end,1);
        %         evolutions(k).trajectories(pp).individual.y_R = R_R_vect(1:end,2);
        %         evolutions(k).trajectories(pp).individual.z_R = R_R_vect(1:end,3);
        %         evolutions(k).trajectories(pp).individual.R_R_norm = R_R_norm(1:end);
        
        lat_c_rad = depvar(1:end,6); lon_c_rad = depvar(1:end,7);
        lat_f_rad = deg2rad(lat_f_deg); lon_f_rad =  deg2rad(lon_f_deg);
        lat_dif = lat_f_rad - lat_c_rad; lon_dif = lon_f_rad - lon_c_rad;
        
        d_rad = acos(sin(lat_c_rad).*sin(lat_f_rad) + cos(lat_c_rad).*cos(lat_f_rad).*cos(lon_dif));
        d_deg = rad2deg(d_rad);
        
        
        chi_req_rad_Y = sin(lon_dif).*cos(lat_f_rad);
        chi_req_rad_X = cos(lat_c_rad)*sin(lat_f_rad) - sin(lat_c_rad).*cos(lat_f_rad).*cos(lon_dif);
        chi_req_rad = atan2(chi_req_rad_Y,chi_req_rad_X);
        chi_req_deg = rad2deg( chi_req_rad );
        chi_err_deg = rad2deg(depvar(1:end,11)) - chi_req_deg;
        
        
        
        evolutions(k).trajectories(p).individual.t    = depvar(1:end,1);
        evolutions(k).trajectories(p).individual.time_vector = depvar(1:end,1) - startEpoch;
        evolutions(k).trajectories(p).individual.x_R = depvar(1:end,2);
        evolutions(k).trajectories(p).individual.y_R = depvar(1:end,3);
        evolutions(k).trajectories(p).individual.z_R = depvar(1:end,4);
        evolutions(k).trajectories(p).individual.R_R_norm = depvar(1:end,5);
        evolutions(k).trajectories(p).individual.latitude_angle  = rad2deg(depvar(1:end,6));
        evolutions(k).trajectories(p).individual.longitude_angle  = rad2deg(depvar(1:end,7));
        evolutions(k).trajectories(p).individual.h    = depvar(1:end,8);
        %    evolutions(k).trajectories(pp).individual.latitude_angle    = rad2deg(depvar(1:end,9));
        %    evolutions(k).trajectories(pp).individual.longitude_angle   = rad2deg(depvar(1:end,10));
        evolutions(k).trajectories(p).individual.heading_angle     = rad2deg(depvar(1:end,9));
        evolutions(k).trajectories(p).individual.flight_path_angle = rad2deg(depvar(1:end,10));
        evolutions(k).trajectories(p).individual.angle_of_attack   = rad2deg(depvar(1:end,11));
        evolutions(k).trajectories(p).individual.angle_of_sideslip = rad2deg(depvar(1:end,12));
        evolutions(k).trajectories(p).individual.bank_angle        = rad2deg(depvar(1:end,13));
        evolutions(k).trajectories(p).individual.d_deg             = d_deg;
        evolutions(k).trajectories(p).individual.heading_error     = chi_err_deg;
        evolutions(k).trajectories(p).individual.heading_required  = chi_req_deg;
                evolutions(k).trajectories(p).individual.airspeed  = depvar(1:end,16);
                evolutions(k).trajectories(p).individual.E  = 9.81*depvar(1:end,8) + 0.5*(depvar(1:end,16)).^2;
            evolutions(k).trajectories(p).individual.interp_E  = interp(1:end,1);
            evolutions(k).trajectories(p).individual.interp_angle_of_attack   = interp(1:end,2);
            evolutions(k).trajectories(p).individual.interp_thrust_angle      = interp(1:end,3);
            evolutions(k).trajectories(p).individual.interp_throttle_setting  = interp(1:end,4);
        
        
        evolutions(k).individuals.gamma_i(p) = evolutions(k).trajectories(p).individual.flight_path_angle(1);
        evolutions(k).individuals.chi_i(p)   = evolutions(k).trajectories(p).individual.heading_angle(1);
        evolutions(k).individuals.lat_f_deg(p) = evolutions(k).trajectories(p).individual.latitude_angle(end);
        evolutions(k).individuals.lon_f_deg(p) = evolutions(k).trajectories(p).individual.longitude_angle(end);
        evolutions(k).individuals.tof(p)       = evolutions(k).trajectories(p).individual.time_vector(end);
        
        
        evolutions(k).fitness.dif_lat(p)   = evolutions(k).trajectories(p).individual.latitude_angle(end) - lat_f_deg;
        evolutions(k).fitness.dif_lon(p)   = evolutions(k).trajectories(p).individual.longitude_angle(end) - lon_f_deg;
        evolutions(k).fitness.dif_d_deg(p) = evolutions(k).trajectories(p).individual.d_deg(end) - 0.75;
        evolutions(k).fitness.dif_h(p)     = evolutions(k).trajectories(p).individual.h(end) - 25000;
        evolutions(k).fitness.tof(p)       = evolutions(k).trajectories(p).individual.time_vector(end);
        
        
        pp = pp + 1;
    end
    [ aaa , idx1 ] = min(evolutions(k).fitness.dif_d_deg);
    [ bbb , idx2 ] = min(evolutions(k).fitness.dif_h);
    [ ccc , idx3 ] = min(evolutions(k).fitness.tof);
    I = [ idx1 idx2 idx3 ];
    criteria = [{'d_deg'} {'h'} {'tof'} ];
    %  I_12 = intersect(idx1,idx2);
    %  I_123 = intersect(I_12,idx3);
    %  I_123 = idx1;
            evolutions(k).max_tof       = max(evolutions(k).individuals.tof);

    for i = 1:3
        evolutions(k).best(i).criteria  = criteria(i);
        evolutions(k).best(i).index     = I(i);
        % evolutions(k).best(i).v_i       = evolutions(k).individuals.v_i(I(i));
        evolutions(k).best(i).gamma_i   = evolutions(k).individuals.gamma_i(I(i));
        evolutions(k).best(i).chi_i     = evolutions(k).individuals.chi_i(I(i));
        evolutions(k).best(i).dif_d_deg = evolutions(k).fitness.dif_d_deg(I(i));
        evolutions(k).best(i).dif_h     = evolutions(k).fitness.dif_h(I(i));
        evolutions(k).best(i).tof       = evolutions(k).fitness.tof(I(i));
    end
    
end
end

