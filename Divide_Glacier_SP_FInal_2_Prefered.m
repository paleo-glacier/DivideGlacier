%2D Glacier model extracted form mark kessler's code gc2d on 12/15/2015 
%Based on Kessler et al., 2006

%MOdified Spring 2016 by S. Pendleton and R. Anderson for use on Divide Ice
%Cap Transect Manuscript

%this iteration optimizes a 2-step climate scenario, inital cooling between
%0-995CE and then another cooling between 995-1910CE immediately followed
%by recent warmiing.  This is the simplest and preferred scenario.  


clear all
figure(1)
clf

beep on

%% Constants %%
        %Mass Balance Constants         
        Tsealevel = -12.7;                %Spin up Sea level temperature
        T_initial_advance = -12.89;        %Initial Advance temperature (0-1000CE)
        T_transect_advance = -13.15;       %Transect advance temperature (1000-1910CE)
        warming_rate = 0.028;              %Rate of modern warming (1910-2015CE)
        albedo = 0.5;                   %From Benn and Evans 2010 (Glaciers and Glaciation, pg. 25)
        hrs_melt = 12;                 % hours of melt per day(?)
        Rmelt_rate =  0.0000035;        %m h-1; melt rate to scale radiative forcing, obtained via model calibration
        Rmelt_factor = Rmelt_rate*hrs_melt; %annual melt per Watt
%         t_day = 80;       % spring equinox; average annual solar rad
%         t_day = 171;      %summer solstice; max solar rad
        t_day = 133;      % Average solar rad during melt season (June 21-Aug. 15)


        maxBz= 0.3 ;                    %maximum accumulation (with no melt)
        
%Other Constants
latitude = 67.1; %for calculating solar flux

    %%%% NUMERICAL and PHYSICAL CONSTANTS
        
        %%% Constants
        g = 9.81;                       % grivational acceleration
        rhoI = 917;                     % density of ice
        rhoW = 1000;                    % density of water
        day = 1/365.25 ;                % length of a day in years

        %%% Time
        t = 0 ;                         % set time to zero
        tMax =50000 ;                   % maximum simulation time in years
        dtMax = 180*day ;       % maximum timestep in years
        dtDefault = 30*day ;    % timestep if VARIABLE_DT_TOGGLE==0
               
        %%%% time for plotting
        nplots = 200;
        tplot = tMax/nplots;
        t_lastplot = 0;
        
        %%% Glacier Properties
        MinGlacThick = 0.25 ;              %meters
         
        conserveIce = 0.0;
        iceVolumeLast =0.0;
        
        
        %%% Ice Deformation
%           glensA = 2.16e-16 ;    %Pa-3 yr-1, Patterson, 1994; MacGregor, 2000 (outdated value)
          glensA = 1.1045e-17 ;   % Cuffey and Patterson 2010 for -10  ice
%           glensA = 2.9332e-17 ;   %Cuffey and Patterson 2010 for -5 degree ice
%           glensA = 5.3618e-17 ;   %Cuffey and Patterson 2010 for -2 degree ice
%           glensA = 7.5696e-17 ;   %Cuffey and Patterson 2010 for 0 degree ice
                
%         %%% Attractor Sliding -- only if ICESLIDE_TOGGLE==1 (generally used)
%         UsChar = 200 ;
%         taubChar = 100000 ;        
        
%         %%% Standard Sliding -- used if ICESLIDE_TOGGLE==2 (generally not used)
%         B = 0.0012 ;                    % m/(Pa*yr) -- MacGregor, 2000
%         DepthToWaterTable = 20 ;        % distance from ice surface to water table
%         MaxFloatFraction = 80 ;         % limits water level in ice
%         Hpeff = 20 ;                    % effective pressure (meters of water)
                 
        %%% Avalanching
        angleOfRepose = 30 ;
        avalanchFreq = 3 ;              % average number per year
                
%         %%% Calving
%         seaLevel = -120 ;               % meters
%         calvingCoef = 2 ;               % year^-1
                        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CODE BEHAVIOR TOGGLES

        %% toggles turn on/off segments of the code or select 
        %% between multiple possibilities for a given process
        %% values can be reset in INIT_COND segment
        
        SAVE_TOGGLE         = 1 ;   % saving                            (off|on)
        PLOT_TOGGLE         = 1 ;   % plotting                          (off|on)
        REPORT_TOGGLE       = 1 ;   % plotting                          (off|on)
        
        COMPRESS_TOGGLE     = 0 ;   % only simulate area with ice       (off|on)
        VARIABLE_DT_TOGGLE  = 1 ;   % state dependent time step         (off|on)

        INIT_COND_TOGGLE    = 1 ;   % load DEM and climate              (synth|valley|sheet)
        GENERIC_ICE_TOGGLE  = 0 ;   % start with generic ice surface    (off|on)       
        
        ICEFLOW_TOGGLE      = 1 ;   % ice motion by deformation         (off|on)
        ICESLIDE_TOGGLE     = 0 ;   % ice motion by sliding             (off|on|select)        
        
        THERMAL_TOGGLE      = 0 ;   % temp dependance of flow           (off|on)
        FREEZEON_TOGGLE     = 0 ;   % basal ice freeze to bed           (off|on)

        AVALANCH_TOGGLE     = 0 ;   % avalanch off steep surfaces       (off|on)
        ERODE_TOGGLE        = 0 ;   % erode the bed                     (off|on|select)
        CALVING_TOGGLE      = 0 ;   % calving front                     (off|on)
        
        
        %%% Available Mass Balance
        ZERO_BALANCE        = 1 ;   % Constant Ice Flux
        CONSTANT_ELA        = 2 ;   % Ice Free Boundary
        ELA_LOWERING        = 3 ;   % Zero Ice Flux
        ELA_TIME_SERIES     = 4 ;   % Continuous Ice Surface Slope
        EXTERNAL_FUNC       = 5 ;   % Constant Surface Elevation
        ELA_LOWERING2       = 6 ;   % Zero Ice Flux
        
        MASS_BALANCE_TOGGLE = CONSTANT_ELA ;    % select climate scenerio   (off|on|select)
        
        
        %%% Available Boundary Conditions
        ICE_FREE_BOUND      = 1 ;   % Ice Free Boundary
        ZERO_FLUX_BOUND     = 2 ;   % Zero Ice Flux
        CONST_FLUX_BOUND    = 3 ;   % Constant Ice Flux
        SURF_ELEV_BOUND     = 4 ;   % Constant Surface Elevation
        SURF_SLOPE_BOUND    = 5 ;   % Continuous Ice Surface Slope
        
        WEST_BC_TOGGLE = ICE_FREE_BOUND    ;   % boundary condition    (no ice|reflect|no flow)
        EAST_BC_TOGGLE = ICE_FREE_BOUND    ;   % boundary condition    (no ice|reflect|no flow)
        SOUTH_BC_TOGGLE = ICE_FREE_BOUND    ;  % boundary condition    (no ice|reflect|no flow)
        NORTH_BC_TOGGLE = ICE_FREE_BOUND    ;  % boundary condition    (no ice|reflect|no flow)
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OUTPUT BEHAVIOR
    
        plotInterval =  5;            % seconds
        saveInterval = 10 ;             % whole years
        reportInterval = 20 ;           % seconds
        
        nextPlot = 0 ;                  % initialize to plot on first timestep
        nextSave = 0 ;                  % initialize to save on first timestep
        nextReport = 0 ;                % initialize to report on first timestep
        
        figureNumber = 0 ;
        zScale = 1.5 ;
        plotContent = 'Extent' ;
        
        outputFile = 'savetmp' ;

                               
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% INITIALIZE COUNTERS
    
         numTimeSteps = 0 ;
         timeSteps = zeros(10000000,1) ;
        
    Thickness_save = [];
    Volume_save = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% INITIALIZE BED and ICE TOPOGRAPHY, and CLIMATE VARIABLES
        
        %%% must define topo, cellsize, dx, and dy
        
            
%             if ( INIT_COND_TOGGLE == 1 )    %% valley glaciers

                filenameDEM = 'Divide_IceStrip_2_Smooth_Clip' ;
                
                load( filenameDEM ) ;
                
                dx = 60 ;  % set a new dx, meters
                dy = dx ;
                
                %%% Mass Balance                 
                initELA = 1200;
                ELA = initELA;
                ELA_last = ELA;

                gradBz = 0.0023 ;

                
%             elseif ( INIT_COND_TOGGLE == 2 )    %% ice sheets
%             
%                 filenameDEM = 'Baffin200d' ;
%                 
%                 load( filenameDEM ) ;
%                             
%                 dx = 2000 ;  % set a new dx
%                 dy = dx ;
%                 
%                 UsChar = 100 ;
%                 taubChar = 50000 ;
%                 
%                 %%% Mass Balance
%                 initELA = 3200 ;
%                 maxBz= 0 ;
%                 gradBz = 1/100 ;
%                 
%                 Hbound = 2000 ;
%                 
%                 Elev0 = 0 ;             % reference elevation
%                 To = -30 ;              % temperature at Elev0
%                 lapseRate = -0.0065 ;   % degrees per meter
%                    
%                 COMPRESS_TOGGLE         = 1 ;
%                 GENERIC_ICE_TOGGLE      = 1 ;
%                 MASS_BALANCE_TOGGLE     = ZERO_BALANCE ;
%                 CALVING_TOGGLE          = 1 ;
%                 ERODE_TOGGLE            = 0 ;
%         
%                 THERMAL_TOGGLE          = 1 ;
%                 FREEZEON_TOGGLE         = 1 ;
%                 HORZTL_ADVECT_TOGGLE    = 1 ;
%                 GEOTHERMAL_HEAT_TOGGLE  = 1 ;
%                 STRAIN_HEAT_TOGGLE      = 1 ;
%                 SLIDING_HEAT_TOGGLE     = 1 ;
%                 SURFACE_HEAT_FLUX_TOGGLE= 1 ;
%                 THERMAL_3D_TOGGLE       = 0 ;
%         
%                 WEST_BC_TOGGLE      = SURF_ELEV_BOUND ;
%                 EAST_BC_TOGGLE      = ZERO_FLUX_BOUND ;
%                 SOUTH_BC_TOGGLE     = ZERO_FLUX_BOUND ;
%                 NORTH_BC_TOGGLE     = ZERO_FLUX_BOUND ;
%                                 
%             end 
     
             
            [rws,cls] = size(topo) ;
            if ( ~exist('easting') )
                easting = 1:cls ;
            end
            if ( ~exist('northing') )
                northing = 1:rws ;
            end
            
            
            % if we need to flip n-s or e-w and dont resize dx do the
            % following
 %               x = (0:clsNew-1)*dx ;
 %                x = fliplr(x);
%                 y = (0:rwsNew-1)*dy ;
%                 y = fliplr(y); % to alter the clearcreek dem
 %               [X, Y] = meshgrid(x,y);
 
  %%%%%Outputs %Simon Pendleton 2/09/16
  
  
    Ela_save = zeros(10000000,1) ;
    Thickness_save = zeros(10000000,1) ;
    Volume_save = zeros(10000000,1) ;
    Coverage_save = [];
    dice_rate = zeros(10000000,1) ;
    Percent_change = zeros(10000000,1) ;
    threshold_save = [];
    Thickness_2015margin = 0;
    Thickness_trimline = 0;
    Thickness_2ka = 0;
    Volume_change = zeros(10000000,1) ;
    areaIce = zeros(10000000,1) ;
    Temp_save = zeros(10000000,1) ;
    T_steadystate = 0;
    Phase_1 = 1;
    Phase_2 = 0;
    Phase_3 = 0;
    Phase_4 = 0;
    Phase_5 = 0;
    Phase_6 = 0;
    Phase_7 = 0;
    T_trimline =0;
    Advance_time = 0;
            %% resample DEM at new node spacing
            if ( cellsize ~= dx )
                
                [rws,cls] = size(topo) ;
                xOld = (0:cls-1)*cellsize ;
                yOld = (0:rws-1)*cellsize ;
                [XOld,YOld] = meshgrid(xOld,yOld);
                
                if ( rem(max(xOld),dx) == 0 & rem(max(yOld),dy) == 0 )
                    clsNew = max(xOld)/dx + 1 ;
                    rwsNew = max(yOld)/dy + 1 ;
                else
                    clsNew = ceil( xOld(end) / dx ) ;
                    rwsNew = ceil( yOld(end) / dy ) ;
                end
                    
                x = (0:clsNew-1)*dx ;
%                 x = fliplr(x);
                y = (0:rwsNew-1)*dy ;
%                 y = fliplr(y); % to alter the clearcreek dem
                [X, Y] = meshgrid(x,y);
            
                topo = interp2( XOld, YOld, topo, X, Y ) ;
%                 easting = interp1( xOld, easting, x ) ;
%                 northing = interp1( yOld, northing, y ) ;
                cellsize = dx ;
                
            end
            
            % Set the bed elevation to 'topo'
            Zb = topo ;
            Zb = Zb-18;
            Zb = double(Zb);
            Zb = transpose(Zb);
            Zb = rot90(Zb,1); %rotates array A counterclockwise by k*90 degrees, where k is an integer.
            initZb = Zb;
            H = zeros(size(Zb)) ;
            Zi = H + Zb ; %turn me back on if you want to run from zero state
            
            %% Calculate PDD
            
            PDD_meltfactor = (6.3)/1000; %melt factor in m day-1
            lapse_rate = -0.0049;           %C m-1 (Gardner et al., 2009)
            Tbar = Tsealevel + lapse_rate*Zi;
%             Tbar = reshape(Tbar, 30798,1);

            DT=20; %annual Temperture delta
            PDD_t = 1:1:365;
            P = 365; %days; annual temperature period
            Temp_curve = DT*sin(2*pi*PDD_t/P);
            
            %analytical solution
            term1 = Tbar.*P.*(1-(1/pi).*acos(Tbar./DT));
            term2 = (DT*P/pi).*sin(acos(Tbar/DT));
            PDD = term1+term2;
            PDD(Tbar>DT) = Tbar(Tbar >DT)*365;
            PDD(-Tbar>DT) = 0;
            Melt_season_length = zeros(size(PDD));
            
            Melt_season_length(PDD > 0) = (log(PDD(PDD>0))/0.03012);
            Melt_season_length(Melt_season_length < 0) = 0;
%             Melt_season_length(Melt_season_length<0) = 0;
            
%             jmax = length(Tbar);
%             for j = 1:jmax
%                 Telev = Tbar(j)+Temp_curve;
%                 PDD(j) = sum(Telev(Telev>0));>DT
%             end
%                 PDD = reshape(PDD, 118,261);
 
            
            %% Calculate aspect, slope, and Solar Radiation
            refvec = [2000 67.1718770638110954 -64.9112130328888810];
            [aspect,slope,gradN,gradE] = gradientm(Zi,refvec);
            ASPECT = 180-aspect; %degrees away from south
            SLOPE = slope; %ice surface slope
            lat = ones(size(Zb)); 
            lat = lat*latitude; %array of latitude values
            elev = Zi; %Ice surface elevations
            
            [Solar_Surf] = Solar_Rad(t_day,lat, elev, ASPECT, SLOPE);
            Solar_Surf = max(Solar_Surf,0);
            
            %% Calculate initial Mass Balance Modulated by Solar Radiation

            Bs_1 = PDD_meltfactor*PDD;
            %= (melt_factor*max((Tsealevel+(lapse_rate*elev)),0));
            Bs_2 = zeros(size(Solar_Surf));
            
            Bs_2(Bs_1 > 0) = (Rmelt_factor.*(Melt_season_length(Bs_1 > 0))*(1-albedo).*Solar_Surf(Bs_1 > 0)); %summer mass balance
            Bs = Bs_1 + Bs_2;
            
            %Winter Mass Balance
            Bw = ones(size(Zb));
            Bw = Bw*maxBz;
            
            %Annual Mass Balance
            Bxy = Bw-Bs;
            
            t_surf_calc = 0;


            clear topo
            
            [rws,cls] = size(Zb) ;
            x = 0:dx:(cls-1)*dx ;
%             x = fliplr(x);
            y = 0:dy:(rws-1)*dy ;
%             y = fliplr(y);
            [X,Y] = meshgrid(x,y);
            
            max(max(Zb));
            
            figure(1);
            clf
           
            surf(X,Y,Zb);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALCULATE MASS BALANCE TIME SERIES
    
        %% load in Delta O18 record and create ELA time series
%         if ( MASS_BALANCE_TOGGLE == ELA_TIME_SERIES )
%             
%             maxELA = 4000 ;
%             minELA = 3350 ;
%             load delO18record
%             dO18 = delO18record(:,2) ;
%             yrBP = delO18record(:,3) + 2006 - 1949 ;
% 
%             n = 5 ;
%             dO18f = filter2( ones(n,1)/n, dO18, 'valid' ) ;
%             ELArecord = flipud(minELA + (maxELA-minELA) * ...
%                     ((dO18f-min(dO18f))/max(dO18f-min(dO18f)))) ;
%             trecord = flipud(yrBP(end) - yrBP) ;
%             trecord = trecord(3:end-2) ;
%             trecord = trecord - min(trecord) ;
% 
%             % spin-up model to steady state
%             tmin = 500 ;
%             trecord = [0 trecord + tmin] ;
%             ELArecord = [ELArecord(1) ELArecord] ;
%                 
%         end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% INITIALIZE RANDOM NUMBER GENERATOR
    
        rand('state',1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SAVE INITIAL SIMULATION STATE
    
        rndstate = rand('state') ;
%         initFilename = 'yearInit0' ;
%         save( initFilename ) ;
        
        
    %% Changes to toggles, parameters or variables
        
        nextPlot = 0 ;      % initialize to plot on next timestep
        nextReport = 0 ;    % initialize to report on next timestep
        
%         load inputArgs      % reload the input args in case they were
%                             % overwritten when inputFile was loaded
        
    
    %% Set random number generator to previous position
    
        rand( 'state', rndstate );
        
       
    %% Save simulation state   
    
%         initFilename = sprintf( 'yearInit%d', round(t) ) ;
%         save( initFilename ) ;
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              START THE TIME LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_last = -50;
nframe = 0;
tic
% 
% % % while Thickness_trimline <0.01 || Thickness_trimline >0.5 && Thickness_trimline_N <0.01 || Thickness_trimline_N >0.5
% while Thickness_trimline <0.01 || Thickness_trimline >1
% while Thickness_2015margin <0.01 || Thickness_2015margin >1
% while Thickness_2ka <0.01 || Thickness_2ka >1
% while Advance_time > 610 || Advance_time < 590
    
    while t<tMax
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS
        
            %%% DEFAULT BOUNDARY CONDITION IS ZERO FLUX
            H_ext = [H(1,1) H(1,:) H(1,end)
                    H(:,1) H H(:,end)
                    H(end,1) H(end,:) H(end,end)] ;
                 
            Zb_ext = [Zb(1,1) Zb(1,:) Zb(1,end)
                    Zb(:,1) Zb Zb(:,end)
                    Zb(end,1) Zb(end,:) Zb(end,end)] ;
        
            Zi_ext = [Zi(1,1) Zi(1,:) Zi(1,end)
                    Zi(:,1) Zi Zi(:,end)
                    Zi(end,1) Zi(end,:) Zi(end,end)] ;
        
            %%% WESTERN BOUNDARY CONDTION
            if ( WEST_BC_TOGGLE == SURF_ELEV_BOUND )            % Constant Ice Surface Height
                ZiBound = mean(Zb(:,1)) + Hbound ;
                H_ext(:,1) = ZiBound - Zb_ext(:,1)  ;
            elseif ( WEST_BC_TOGGLE == CONST_FLUX_BOUND )         % Constant Ice Flux B.C.
            elseif ( WEST_BC_TOGGLE == SURF_SLOPE_BOUND )       % Constant Ice Surface Slope
                Zi_ext(:,1) = 2*Zi_ext(:,2) - Zi_ext(:,3) ;
                H_ext(:,1) = Zi_ext(:,1) - Zb_ext(:,1) ;
                H_ext(:,1) = max( 0, H_ext(:,1) ) ;
            elseif( WEST_BC_TOGGLE == ICE_FREE_BOUND )          % Ice Free Boundary
                H_ext(:,1) = 0 ;
            end

            %%% EASTERN BOUNDARY CONDTION
            if ( EAST_BC_TOGGLE == SURF_ELEV_BOUND )            % Constant Ice Surface Height
                ZiBound = mean(Zb(:,end)) + Hbound ;
                H_ext(:,end) = ZiBound - Zb_ext(:,end)  ;
            elseif ( EAST_BC_TOGGLE == CONST_FLUX_BOUND )         % Constant Ice Flux B.C.
            elseif ( EAST_BC_TOGGLE == SURF_SLOPE_BOUND )       % Constant Ice Surface Slope
                Zi_ext(:,end) = 2*Zi_ext(:,end-1) - Zi_ext(:,end-2) ;
                H_ext(:,end) = Zi_ext(:,end) - Zb_ext(:,end) ;
                H_ext(:,end) = max( 0, H_ext(:,end) ) ;
            elseif( EAST_BC_TOGGLE == ICE_FREE_BOUND )          % Ice Free Boundary
                H_ext(:,end) = 0 ;
            end
            
            %%% SOUTHERN BOUNDARY CONDTION
            if ( SOUTH_BC_TOGGLE == SURF_ELEV_BOUND )           % Constant Ice Surface Height
                ZiBound = mean(Zb(1,:)) + Hbound ;
                H_ext(1,:) = ZiBound - Zb_ext(1,:)  ;
            elseif ( SOUTH_BC_TOGGLE == CONST_FLUX_BOUND )        % Constant Ice Flux B.C.
            elseif ( SOUTH_BC_TOGGLE == SURF_SLOPE_BOUND )      % Constant Ice Surface Slope
                Zi_ext(1,:) = 2*Zi_ext(2,:) - Zi_ext(3,:) ;
                H_ext(1,:) = Zi_ext(1,:) - Zb_ext(1,:) ;
                H_ext(1,:) = max( 0, H_ext(1,:) ) ;
            elseif( SOUTH_BC_TOGGLE == ICE_FREE_BOUND )          % Ice Free Boundary
                H_ext(1,:) = 0 ;
            end
            
            %%% NORTHERN BOUNDARY CONDTION
            if ( NORTH_BC_TOGGLE == SURF_ELEV_BOUND )           % Constant Ice Surface Height
                ZiBound = mean(Zb(end,:)) + Hbound ;
                H_ext(end,:) = ZiBound - Zb_ext(end,:)  ;
            elseif ( NORTH_BC_TOGGLE == CONST_FLUX_BOUND )        % Constant Ice Flux B.C.
            elseif ( NORTH_BC_TOGGLE == SURF_SLOPE_BOUND )      % Constant Ice Surface Slope
                Zi_ext(end,:) = 2*Zi_ext(end-1,:) - Zi_ext(end-2,:) ;
                H_ext(end,:) = Zi_ext(end,:) - Zb_ext(end,:) ;
                H_ext(end,:) = max( 0, H_ext(end,:) ) ;
            elseif( NORTH_BC_TOGGLE == ICE_FREE_BOUND )         % Ice Free Boundary
                H_ext(end,:) = 0 ;
            end
        
            Zi_ext = Zb_ext + H_ext ;
                        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE THE BASAL SHEAR STRESS
        

            
            % forward differences
            dZidxX_ext = ( Zi_ext(:,2:end) - Zi_ext(:,1:end-1) ) / dx ;
            dZidyY_ext = ( Zi_ext(2:end,:) - Zi_ext(1:end-1,:) ) / dy ;
            dZidxX = dZidxX_ext(2:end-1,:) ;
            dZidyY = dZidyY_ext(:,2:end-1) ;

            HX_ext = ( H_ext(:,2:end) + H_ext(:,1:end-1) ) / 2 ;
            HY_ext = ( H_ext(2:end,:) + H_ext(1:end-1,:) ) / 2 ;
            HX = HX_ext(2:end-1,:) ;
            HY = HY_ext(:,2:end-1) ;
            
            taubxX_ext = -rhoI * g * HX_ext .* dZidxX_ext ;
            taubyY_ext = -rhoI * g * HY_ext .* dZidyY_ext ;
            
            taubxX = taubxX_ext(2:end-1,:) ;
            taubyY = taubyY_ext(:,2:end-1) ;
            
            taubxY = ( taubxX_ext(1:end-1,1:end-1) + taubxX_ext(1:end-1,2:end) + ...
                       taubxX_ext(2:end,1:end-1) + taubxX_ext(2:end,2:end) ) / 4;
            
            taubyX = ( taubyY_ext(1:end-1,1:end-1) + taubyY_ext(1:end-1,2:end) + ...
                       taubyY_ext(2:end,1:end-1) + taubyY_ext(2:end,2:end) ) / 4;
            
            taubX = sqrt( taubxX.^2 + taubyX.^2 ) ;
            taubY = sqrt( taubxY.^2 + taubyY.^2 ) ;
       
            taubX = taubX .* ( HX > 0 ) ;
            taubY = taubY .* ( HY > 0 ) ;
            
            % fill in zero values with 1 for use in division
            taubOnesX = taubX ;
%             taubX(imag(taubX) ~= 0) = 0; %remove complex numbers
            taubNot_indX = find(~taubX) ;
            taubOnesX(taubNot_indX) = 1 ;
            taubOnesY = taubY ;
%             taubY(imag(taubY) ~= 0) = 0; %remove complex numbers
            taubNot_indY = find(~taubY) ;
            taubOnesY(taubNot_indY) = 1 ;
            
            xcmpnt = taubxX./taubOnesX ;
            xcmpnt(taubNot_indX) = 0 ;
            
            ycmpnt = (taubyY./taubOnesY) ;
            ycmpnt(taubNot_indY) = 0 ;
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE ICE VELOCITY DUE TO DEFORMATION
        
        if ( ICEFLOW_TOGGLE )
            
                AX = glensA ;
                AY = glensA ;

            %% here's the guts of calculating the depth averaged velocity
            UdxX = abs( (2/5) * AX .* (taubX.*taubX.*taubX) .* HX ) .* xcmpnt ;
            UdyY = abs( (2/5) * AY .* (taubY.*taubY.*taubY) .* HY ) .* ycmpnt ;
                    
        else  % need variables filled with zero values
        
            [rws,cls] = size(Zb) ;
            UdxX = zeros(rws,cls+1) ;
            UdyY = zeros(rws+1,cls) ;
            
        end %- ICEFLOW_TOGGLE -%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE SLIDING VELOCITY
        
%         if ( ICESLIDE_TOGGLE == 1  )        %%% ATTRACTOR FORMULATION
%        
%             if ( THERMAL_TOGGLE && FREEZEON_TOGGLE )
%             
%                 notFrozen = ( Tb_ext > -0.5 | Zb_ext < seaLevel );
%                 notFrozenX = ( notFrozen(2:end-1,1:end-1) + notFrozen(2:end-1,2:end) )/2 ;
%                 notFrozenY = ( notFrozen(1:end-1,2:end-1) + notFrozen(2:end,2:end-1) )/2 ;
%                
%                 %% here's the guts of calculating the sliding velocity 
%                 UsxX = notFrozenX .* UsChar .* exp(1 - taubChar./taubOnesX) .* xcmpnt ;
%                 UsyY = notFrozenY .* UsChar .* exp(1 - taubChar./taubOnesY) .* ycmpnt ;
%                 
%             else
%             
%                 %% here's the guts of calculating the sliding velocity
%                 UsxX = UsChar .* exp(1 - taubChar./taubOnesX) .* xcmpnt ; 
%                 UsyY = UsChar .* exp(1 - taubChar./taubOnesY) .* ycmpnt ; 
%                                                   
%             end             
%                         
%         else  % need variables filled with zero values
        
            [rws,cls] = size(Zb) ;
            UsxX = zeros(rws,cls+1) ;
            UsyY = zeros(rws+1,cls) ;

%         end %- ICESLIDE_TOGGLE -%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% MASS CONSERVATION -- CONTINUITY

            % ensure that no ice is drawn from the rock
            CLASS = ( H_ext >= MinGlacThick ) ;
            
            DCLASSx = ( CLASS(2:end-1,2:end) - CLASS(2:end-1,1:end-1) ) .* ...
                        sign( dZidxX ) ;
                     
            DCLASSy = ( CLASS(2:end,2:end-1) - CLASS(1:end-1,2:end-1) ) .* ...
                        sign( dZidyY ) ;
            
            % sum all contributions to ice motion
            UxX = UdxX + UsxX ;
            UyY = UdyY + UsyY ;
                        
            ind = find( DCLASSx == -1 ) ;
            UxX(ind) = 0;
                        
            ind = find( DCLASSy == -1 ) ;
            UyY(ind) = 0;
            
            % calculate both components of the ice flux
            qxX = UxX .* HX ;
            qyY = UyY .* HY ;
            
%             if ( WEST_BC_TOGGLE == CONST_FLUX_BOUND )
% 				qxX(:,1) = BoundaryFlux ;
%             end
%             
%             if ( EAST_BC_TOGGLE == CONST_FLUX_BOUND )
% 				qxX(:,end) = BoundaryFlux ;
%             end
%             
%             if ( SOUTH_BC_TOGGLE == CONST_FLUX_BOUND )
% 				qyY(1,:) = BoundaryFlux ;
%             end
%             
%             if ( NORTH_BC_TOGGLE == CONST_FLUX_BOUND )
% 				qyY(end,:) = BoundaryFlux ;
%             end
            
            % here's the guts of the continuity equation
            dqdxX = ( qxX(:,2:end) - qxX(:,1:end-1) ) / dx ;
            dqdyY = ( qyY(2:end,:) - qyY(1:end-1,:) ) / dy ;
            dHdt = -dqdxX -dqdyY ;
            
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE MASS BALANCE
            
            %% the imposed mass balance is the imposed climate
            %% there are many possibilities, here are only a few
            %% all must populate the 2D matrix Bxy of size = size(Zb)
            %% with values of net precip/melt rate in m/yr
            %% define the scalar, ELA (m), for plotting
            
            if t>t_surf_calc +1;
                t_surf_calc = t;
                
            [Solar_Surf] = Solar_Rad(t_day, lat, elev, ASPECT, SLOPE);
            Solar_Surf = max(Solar_Surf,0);
            end
            
            %% Calculate PDD
            Tbar = Tsealevel + lapse_rate*Zi;
%             Tbar = reshape(Tbar,30798,1);
%             jmax = length(Tbar);

            %analytical solution
            term1 = Tbar.*P.*(1-(1/pi).*acos(Tbar./DT));
            term2 = (DT*P/pi).*sin(acos(Tbar/DT));
            PDD = term1+term2;
            PDD(Tbar>DT) = Tbar(Tbar>DT)*365;
            PDD(-Tbar>DT) = 0;
            Melt_season_length = zeros(size(PDD));
            Melt_season_length(PDD > 0) = (log(PDD(PDD>0))/0.03012);
            Melt_season_length(Melt_season_length < 0) = 0;
            
%             for j = 1:jmax
%                 Telev = Tbar(j)+Temp_curve;
%                 PDD(j) = sum(Telev(Telev>0));
%             end
%                 PDD = reshape(PDD, 118,261);
            %%Mass Balance Modulated by Solar Radiation

            Bs_1 = PDD_meltfactor*PDD;
%             Bs_1 = real(Bs_1);
            Bs_2 = zeros(size(Solar_Surf));
            
            Bs_2(Bs_1 > 0) = (Rmelt_factor.*Melt_season_length(Bs_1 > 0)*(1-albedo).*Solar_Surf(Bs_1 > 0)); %summer mass balance
            Bs = Bs_1 + Bs_2;
            
            %Winter Mass Balance
            Bw = ones(size(Zb));
            Bw = Bw*maxBz;
            
            Bxy = Bw-Bs;
%             end
%             else
%                 Bxy = Bxy;
%             end
            
%             ELA = (elev(max(find(abs(Bxy)< 0.001)))+elev(min(find(abs(Bxy)< 0.001))))/2; %mean ELA
            
%             Bs_south = Bs(find(abs(Solar_Surf)<5))   Zi(find(abs(Solar_Surf)<5))
            
%             if ( MASS_BALANCE_TOGGLE == CONSTANT_ELA )
            % Simple ELA, maxBz, gradBz
%             ELA = initELA ;
%             Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
                    
%                     %Aspect Scaling
%                if t > 10 & ELA == ELA_last;
%                    Bxy = Bxy;
%                else
%                 ELA = initELA ;
%                 Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
%                 Bxy = Bxy - (abs(Bxy).* aspect);
%                end
            
%             elseif ( MASS_BALANCE_TOGGLE == ELA_LOWERING )
%             % ELA changing with time experiment
%             
%                 % ELAStepSize = -10 ;       % positive/negative values raise/lower ELA
%                 % ELAStepInterval = 500 ;
%                 
%                 ELA = initELA + ELAStepSize * max( 0, floor( t/ELAStepInterval ) ) ;
%                 Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
%             
%             elseif ( MASS_BALANCE_TOGGLE == ELA_LOWERING2 )
%             % ELA changing with time experiment
%             
%                 tau = 30 ;              % intrinsic timescale of ice dynamics 
%                 tmin = 0 ;              % time to begin ELA modification
%                 initELA = 1380 ;        % initial ELA
%                 stepSize = -20 ;        % positive/negative values raise/lower ELA
%                 dELAdt = -0.1 ;
%                 
%                 ELA = initELA + stepSize * max( 0, floor( (t-tmin) / (8*tau) ) ) ;
%                 Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
%             
%             elseif ( MASS_BALANCE_TOGGLE == EXTERNAL_FUNC )
%             % external mass balance function
%                 
%                 if ( ~exist('Bxy') || t >= nextGetBxy )
%             
%                     % Mass Balance 2D Must Return Bxy (2d Matrix)
%                     Bxy = mass_balance_gc2d( t, cellsize, Zi ) ;
%                     nextGetBxy = t + getBxyInterval ;
%                     
%                 end
%             
%             elseif ( MASS_BALANCE_TOGGLE == ELA_TIME_SERIES )
%             % ELA time series
% 
%                 ELA = interp1( trecord, ELArecord ) ;
%                 Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
%         
%             elseif ( MASS_BALANCE_TOGGLE == ZERO_BALANCE )
%             
%                 ELA = 0 ;
%                 Bxy = zeros( size(Zb) ) ;
%                 
%             else
%             
%                 error('Unrecognized Mass Balance')
%                 
%             end    %- MASS_BALANCE_TOGGLE -%
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE TIMESTEP
        
        if ( VARIABLE_DT_TOGGLE == 1 )
            
            %% now that we know the rate of change in ice surface heights due to  
            %% ice motion and due to precipitation or melt we need to know over 
            %% what period of time we can project forward with these rates and 
            %% maintain stability of the ice surface.  The basic idea here is that
            %% we don't want to take a timestep any longer then it would take to 
            %% reverse the ice surface slope between two cells, such that ice 
            %% should be flowing in the other direction.  In fact, let's make our 
            %% timestep much less then that.
            
            %% this calculation sets the timestep such that the change
            %% in ice surface elevation nowhere exceeds a set fraction
            %% of the local standard deviation in ice surface elevations
            
            % include ice changes by precip and melt
            dHdtTot = dHdt + Bxy ;
            adHdt = abs(dHdtTot) ;
            
            % something like standard deviation of 3x3 cell areas around each cell
            filt = [1 1 1;1 1 1;1 1 1] / 9 ;
            ZiMean = filter2( filt, Zi_ext, 'valid' ) ;
            dHmax = sqrt( filter2( filt, (ZiMean - Zi).^2 ) ) ;
            
            % only consider cells with ice thickness > 10 m
            isGlac = H > 10 ;
            
            % find limiting timestep for each considered cell
            ind = find( adHdt & dHmax & isGlac ) ;    
            dtLimits = dHmax(ind)./adHdt(ind) ;
            [dt, idt] = min( dtLimits ) ;
            
            % locate the x and y position of limiting cell for plotting
            [rwDT,clDT] = ind2sub( size(adHdt), ind(idt) ) ; 
            
            % limit timestep to dtMax or some fraction of the calculated timestep
            dt = min( dtMax, dt ) ;
            
            % catch an error, (e.g. if H<10 in all cells )
            if isempty(dt)  
                dt = dtDefault ;
            end
            
        else

            dt = dtDefault ;
                        
        end %- VARIABLE_DT_TOGGLE -%
        
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% UPDATE the TIME and ICE THICKNESS
        
            %% time update
            t = t + dt ;
            numTimeSteps = numTimeSteps + 1 ;
            timeSteps(numTimeSteps) = dt ;
            t_surf_calc = t;

            
            %% increase in ice thicknesses due to precip
            Bxy_pos = Bxy .* ( Bxy > 0 ) ;
            H = H + (Bxy_pos * dt) ;
            
            %% change ice thicknesses due to ice motion
            H = H + (dHdt*dt) ;
            
            %% decrease in ice thicknesses due to melt
            Bxy_neg = Bxy .* ( Bxy < 0 ) ;
            Bxy_neg = -min( H, -Bxy_neg ) ;
%             Bxy_neg = Bxy_neg .* aspect;
            H = H + (Bxy_neg * dt) ;
            
            %% record ice addition or removal by climate
            snowFall = ( Bxy_neg + Bxy_pos ) * dt ;
            conserveIce = conserveIce + sum(sum(snowFall));
            
            %% record ice flux through boundaries
            qbound = sum(qyY(1,:)) - sum(qyY(end,:)) + sum(qxX(:,1)) - sum(qxX(:,end)) ;
            conserveIce = conserveIce + dt * qbound / dx ;
            
            Zi = Zb + max( H, 0 );

            
            %Calculate aspect with new ice surface
            if t>t_surf_calc +1;
                [aspect,slope,gradN,gradE] = gradientm(Zi,refvec);
                ASPECT = 180-aspect; %degrees away from south
                SLOPE = slope; %ice surface slope
                elev = Zi; %Ice surface elevations
            end
                
            if ( isnan(Zi) )
                save workspacedump
                error('NaN in ice thickness') ;
            end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% THIS IS THE END OF THE CONTINUUM CALCULATIONS 
    %%% NOW SIMULATE PROCESSES FOR WHICH WE HAVE NO EQUATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% AVALANCHE SNOW OFF OF STEEP SURFACES
        
        if ( AVALANCH_TOGGLE && ( rand < dt*avalanchFreq ) )
            
            %% move ice downslope until the ice surface is everywhere
            %% less then or near the angle of repose
            
            [rws,cls] = size(Zb) ;
            dHRepose = dx*tan(angleOfRepose*pi/180) ;
            Ho = max( H, 0 ) ;
            
            while ( 1 )

                dZidx_down = zeros(rws,cls);
                dZidx_down(:,2:end) = max( 0, ( Zi(:,2:end) - Zi(:,1:end-1) ) ) ;
                dZidx_up = zeros(rws,cls);
                dZidx_up(:,1:end-1) = max( 0, ( Zi(:,1:end-1) - Zi(:,2:end) ) ) ;
                dZidx = max( dZidx_up, dZidx_down ) ;

                dZidy_left = zeros(rws,cls);
                dZidy_left(2:end,:) = max( 0, ( Zi(2:end,:) - Zi(1:end-1,:) ) ) ;
                dZidy_right = zeros(rws,cls);
                dZidy_right(1:end-1,:) = max( 0, ( Zi(1:end-1,:) - Zi(2:end,:) ) ) ;
                dZidy = max( dZidy_left, dZidy_right ) ;

                grad = sqrt( dZidx.^2 + dZidy.^2 );
                gradT =  dZidy_left + dZidy_right + dZidx_down + dZidx_up ;
                gradT(find(gradT==0)) = 1;
                grad( find(Ho < 0.1) ) = 0 ;

                mxGrad = max(max( grad ) ) ;
        
                if ( mxGrad <= 1.1*dHRepose )
                    break ;
                end

                delH = max( 0, ( grad - dHRepose ) / 3 ) ;
        
                Htmp = Ho ;        
                Ho = max( 0, Htmp - delH );
                delH = Htmp - Ho ;
        
                delHdn = zeros(rws,cls) ; delHup = zeros(rws,cls) ;
                delHlt = zeros(rws,cls) ; delHrt = zeros(rws,cls) ;
        
                delHup(:,2:end) = delH(:,1:end-1) .* dZidx_up(:,1:end-1)./gradT(:,1:end-1) ;
                delHdn(:,1:end-1) = delH(:,2:end) .* dZidx_down(:,2:end)./gradT(:,2:end) ;
                delHrt(2:end,:) = delH(1:end-1,:) .* dZidy_right(1:end-1,:)./gradT(1:end-1,:) ;
                delHlt(1:end-1,:) = delH(2:end,:) .* dZidy_left(2:end,:)./gradT(2:end,:) ;
        
                Ho = max( 0, Ho + delHdn + delHup + delHlt + delHrt ) ;
        
                Zi = Zb + Ho ;
        
            end
            
            H = Ho + (H<0).*H ;
            
        end %- AVALANCH_TOGGLE -%
        
               iceVolume = sum(sum(H.*(H>0)))*dx*dy ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REPORT SOME STUFF
        
        if ( REPORT_TOGGLE &  toc >= nextReport )
                
            disp( sprintf('\nelapsed time: %1.2f s', toc ) ) ;
            disp( sprintf('simulation time: %1.2f yr', t) ) ;
            disp( sprintf('timestep: %1.2e yr', dt) ) ;
            disp( sprintf('ELA: %1.2e m', ELA) ) ;
            
            %% fractional ice conservation
            
            disp( sprintf('total ice: %1.2e km^3', iceVolume/1e9 ) ) ;
            disp( sprintf('excess ice: %1.2f m^3', (( iceVolume - conserveIce*dx*dy )) ) ) ;
            disp( sprintf('ice change: %f m^3', iceVolume - iceVolumeLast ) ) ;
            disp( sprintf('ice conservation (%%): %1.15f', 100 - 100*( iceVolume - conserveIce*dx*dy ) / iceVolume ) ) ;
            disp( sprintf('max ice thickness: %1.2e m', max(max(H)) ) ) ;
            iceVolumeLast = iceVolume ;
            
            
            nextReport = toc + reportInterval ;
            
        end %- REPORT -%
        
        
         %%%%append variables to arrays

         areaIce(numTimeSteps) = (dx*dy)*(length(H(H>MinGlacThick)));       
        Ela_save(numTimeSteps) = ELA;
        Thickness_save(numTimeSteps) = H(88,158);
%         if Volume_save(end) ~= icevolume
        Volume_save(numTimeSteps) = iceVolume;
        Timestep_save(numTimeSteps) = t;
        Temp_save(numTimeSteps) = Tsealevel;

                Elv_range = (1150:1:1600);
    nmax = length(Elv_range);
            for n = 1:nmax; 
                Mass_balance(n) = maxBz*((Elv_range(n)-ELA)*gradBz); 
            end
            
            Mass_balance(Mass_balance>maxBz)=maxBz;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PLOT SOME STUFF
        tsince = t-t_lastplot;
        if( toc > nextPlot | t>tMax)
            
        ix = X(X==(60*157));%210
        iy = Y(X==(60*157)); 
        iz = Zi(X==(60*157));
        izb = Zb(X==(60*157));
        
%         if (t == dt | tsince>tplot | t>tMax)
            nframe = nframe+1;
        figure(2) %3d Surface plot
            clf
%             subplot(2,2,1)
            titleText = 'Map View' ;
            colorscheme = othercolor('Blues9',200);
            colormap(colorscheme);

            hold on
            plot3(X(84,158)/1000,Y(84,158)/1000,Zi(84,158)+200,'s','MarkerSize',10,'Color','k','MarkerFaceColor','w');%2015 Trimline
            plot3(X(88,158)/1000,Y(88,158)/1000,Zi(88,158)+200,'d','MarkerSize',10,'Color','k','MarkerFaceColor','w');%2015 Trimline
            plot3(X(60,60)/1000,Y(60,60)/1000,Zi(60,60)+200,'o','MarkerSize',10,'Color','k','MarkerFaceColor','w');%2 ka position
%             plot3(ix/1000,iy/1000,iz+100,'Color', 'k','LineWidth', 2); %Transect line
            mesh(X/1000,Y/1000, Zi) %Ice surface
            mesh(X/1000,Y/1000, Zb,'Edgecolor',[0.7843,0.7843,0.7843]) %Bedrock Surface
            [C,h] = contour3(X/1000,Y/1000,Zb,30,'LineColor', [0 0 0]); %Bedrock Contours
            box
            hold off

            view([0,90]); % for view of divide ice cap
%             camroll(45);
            set(gca,'Visible','on')
            set(gca,'FontSize',16)
            set(gca,'FontSize',16)
            xlabel('km','fontsize',21)
            ylabel('km','fontsize',21)


            h = title(titleText,'fontname','arial','fontsize',24) ;
            set( h, 'visible', 'on' ) ;
%             ht=text(.01*max(max(X))/1000,.95*max(max(Y))/1000,6*max(max(Zb)),['Simulation Time = ',num2str(floor(t)), ' years'],'fontsize',21,'BackgroundColor', 'w'); %Year counter
%             text(.01*max(max(X))/1000,.85*max(max(Y))/1000,6*max(max(Zb)),['M_{R} = ',num2str(Rmelt_rate*12), ' mm day^{-1}'],'fontsize',16,'BackgroundColor', 'w'); %melt rate text
            text(1.01*max(max(ix))/1000,0.03*max(max(iy))/1000,3*max(max(Zb)),['S'],'fontsize',21)% North end of Transect
            text(1.01*max(max(ix))/1000,0.95*max(max(iy))/1000,3*max(max(Zb)),['N'],'fontsize',21)% North end of Transect
            
            
%             figure(1) %Cross-section
%                 subplot(2,2,3)
%                 set(gca,'FontSize',16)
%                 set(gca,'FontSize',16)
%                 titleText = 'Transect Cross-section' ;
%                 title(titleText,'fontname','arial','fontsize',24) ;
%                 xlabel('Meters', 'fontsize', 21);
%                 ylabel('Elevation (m)', 'fontsize', 21);
%                 hold on
%                 plot(iy,iz, 'color', 'b')
%                 plot(iy,izb, 'color', 'k')
%                 area(iy, iz,'FaceColor', 'b')
%                 area(iy,izb,'FaceColor', [0.8 0.8 0.8]);
%                 plot(4980, Zb(84,158),'s','MarkerSize',10,'Color','k','MarkerFaceColor','w'); %2015 Ice margin
%                 plot(5220, Zb(88,158),'d','MarkerSize',10,'Color','k','MarkerFaceColor','w');%Local trim line
%                 xlim([0 7000])
%                 ylim([1150 1600])
%                 text(max(iy)-800,(max(max(izb))+115),['North'],'fontsize',21)% North end of Transect
%                 text(min(iy)+200,(max(max(izb))+115),['South'],'fontsize',21)% North end of Transect
%                 text(0.70*(max(max(iy))),(max(max(izb))+75),['M_{R} = ',num2str(Rmelt_rate*1000*12), ' mm day^{-1}'],'fontsize',16,'BackgroundColor', 'w');
%                 text(0.70*(max(max(iy))),(max(max(izb))+30),['M_{T} = ',num2str(PDD_meltfactor*1000), ' mm day^{-1} ^{\circ}C^{-1}'],'fontsize',16,'BackgroundColor', 'w');
%                 text(0.70*(max(max(iy))),(max(max(izb))-10),['Maxbz = ',num2str(maxBz), ' mwe year^{-1}'],'fontsize',16,'BackgroundColor', 'w');
%                 box
%                 
%                 hold off
%                 
%                 %%%build cumulative time series
%                 % now calculated above
%                 Timestep_save = cumsum(timeSteps(1:sum(timeSteps(:,:)~=0)));
% 
% %                 figure(1) %ELA and Ice Volume over time
% %                     subplot(2,3,3)
% %                     hold on
% %                     x1 = Timestep_save;
% %                     y1 = Ela_save(1:numTimeSteps);
% %                     y2 = Thickness_save(1:numTimeSteps);
% %                     y3 = Volume_save(1:numTimeSteps)/1e9;
% %                     [ax h1 h2] =plotyy(x1,y1,x1,y3);
% %                     titleText = 'ELA and Ice Volume' ;
% %                     title(titleText,'fontname','arial','fontsize',24) ;
% %                     ylabel(ax(1),'ELA (m)', 'fontsize', 18);
% %                     ylabel(ax(2),'Ice Volume (km^3)', 'fontsize', 18);
% %                     xlabel('Time (yr)','fontsize',18,'fontname','arial')
% %                     set(gca,'fontsize',14,'fontname','arial')
% %                     hold off
%                     
% 
%                     
                figure(1) %Mass Balance
%                     subplot(2,2,4)
                    plot(Bw(find(abs(ASPECT)>175)), Zi(find(abs(ASPECT)>175)),'Color', 'k'); %North slopes Winter balance (same everywhere)
                    set(gca,'FontSize',16)
                    hold on
                    %Annual Mass Balance Components
                    titleText = 'Mass Balance' ;
                    title(titleText,'fontname','arial','fontsize',24) ;
                    
                    scatter(Bxy(find(abs(ASPECT)>179)), Zi(find(abs(ASPECT)>179)),'filled'); %Annual Balance North Slopes
                    scatter(Bxy(find(abs(ASPECT)<1)), Zi(find(abs(ASPECT)<1)),'filled','d'); %Annual Balance South Slopes
                    line(zeros(length(Elv_range)), Elv_range,'linewidth',2); %Zero mass balance
                    legend({'Winter Accumulation','North Slopes','South Slopes'},'Location','NorthWest','FontSize',18);
%                     plot(Mass_balance,ELA*ones(size(Mass_balance)),'g--','linewidth',1);
                    xlabel('Mass Balance (m/yr)', 'fontsize', 21);
                    ylabel('Elevation (m)', 'fontsize', 18);
                    title('Annual Mass Balance','fontname','arial','fontsize',24);
%                     text((min(min(Mass_balance))*0.9),(ELA+10),'ELA', 'fontsize', 14);
                ylim([1150 1600])
                xlim([-1 0.5])
%                 box
                    hold off
%         if t>2            
%       figure(1) %ice volume and area change
%             subplot(2,2,2)
%             set(gca,'FontSize',16)
%             set(gca,'FontSize',16)
%             hold on
%             titleText = 'Volume Change and Temperature' ;
%             title(titleText,'fontname','arial','fontsize',24) ;
%             x1 = Timestep_save(1:(numTimeSteps-1));
%             y1 = dice_rate(1:(numTimeSteps-1))/1000000;
% %             y2 = areaIce(1:(numTimeSteps-1))/1000000; %ice area
%             y2 = Temp_save(1:(numTimeSteps-1)); %temp history
%             [axx, hh1, hh2] = plotyy(x1,y1,x1,y2);
%         %       plot(Timestep_save(1:(numTimeSteps-1)),dice_rate(1:(numTimeSteps-1)));
%         %       plot(Timestep_save(1:(numTimeSteps-1)),areaIce(1:(numTimeSteps-1)));
% 
%             ylabel(axx(1),'Volume Change (1000 km^{3} yr^{-1})','fontname','arial','fontsize',18)
% %             ylabel(axx(2),'Ice Area (km^{2})','fontname','arial','fontsize',18) %ice area
%             ylabel(axx(2),'Sea Level MAT (^{\circ}C)','fontname','arial','fontsize',18)
%             xlabel('Time (yr)', 'fontsize', 18);
% 
% 
%             hold off
%             box
%         %       set(gca,'fontsize',14,'fontname','arial')
%         end

      %drawnow
            
%            hold on
            
%             if(t>100)
%             colorscheme2 = othercolor('PuBu3',200);
%             %set(gcf,'DefaultAxesColorOrder',colorscheme); 
%              %surf(X/1000,Y/1000,Zi,'EdgeColor','none','FaceColor','texturemap');
%             surf(X(H>0.5)/1000,Y(H>0.5)/1000,Zi(H>0.5))
%             %view([65,75]);
%             view([3,76]); % for view up arkansas river
%             colormap(colorscheme2);
%             set(gca,'Visible','off')
%             end
            
%             CONTOUR_INTERVAL = 100;
%             V = min(min(Zb)):CONTOUR_INTERVAL*zScale:max(max(Zi)) ;
%             c = contourc( X(1,:)/1000, Y(:,1)/1000, Zi, V ) ;
%             mask = H > MinGlacThick ;
%             mask = double(mask);
%             i = 1 ;
%             while i < length(c)
%             
%                 cz = c(1,i) ;
%                 num = c(2,i) ;
%                 cnext = i+num+1;
%                 
%                 glac = interp2( X/1000, Y/1000, mask, c(1,i+1:i+num), c(2,i+1:i+num),'spline' ) ;
%                 if( sum(glac) >=2 )
%                     for j = 1:num-1
%                         i = i+1 ;
%                         if ( glac(j)+glac(j+1) == 2 )
% 	                        plot3( c(1,i:i+1), c(2,i:i+1), [cz cz], 'y', 'linewidth', 1 ) ;
%                         end
%                     end
%                 end
%                 
%                 i = cnext ;
%             end
            drawnow
        %        
            M(:,nframe) = getframe(gcf); 
            
            nextPlot = toc + plotInterval ;
            t_lastplot = t;
            
        end
        

        
        
%         if ( PLOT_TOGGLE & toc >= nextPlot)
%             
%             
%             figure(1)
%             clf
%             surf(X,Y,H);
%             
%             rndstate = rand('state') ;
%             save( outputFile )
%             [az,el] = view ;
%             plot_gc2d( 'savetmp', 'FIGURE_NUMBER', figureNumber, 'zScale', zScale, ...
%                         'SUBFIGURE', 1, 'PLOT_CONTENT', plotContent ) ;
%             
%             nextPlot = toc + plotInterval ;
%             
%         end %- PLOT -%
        
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SAVE SOME STUFF 
        %%continue or break loop if SS is achieved (by change in ice volume)
        
    Thickness_2015margin =(Zi(84,158)-Zb(84,158));
    Thickness_trimline = (Zi(88,158)-Zb(88,158));
    Thickness_trimline_N = (Zi(95,158)-Zb(95,158));
    Thickness_2ka = (Zi(60,60)-Zb(60,60));
        
  if t > 1
    Volume_change(numTimeSteps) = ((Volume_save(numTimeSteps)-Volume_save(numTimeSteps-1)));
    dice_rate(numTimeSteps) = Volume_change(numTimeSteps)/dt;
    Percent_change(numTimeSteps) = abs(Volume_change(numTimeSteps))/Volume_save(numTimeSteps); %scales against an assumed max volume of ice
  end
  
  %%Finding SS%%
%     if t >200 && abs(dice_rate(numTimeSteps)) < 0.01*max(max(dice_rate(1:numTimeSteps))) && Phase_1 == 1;
% %         T_steadystate = t;
% %         Phase_1 = 0;
% %         Phase_2 =1;
% %         Tsealevel = Tsealevel+2;
%           break
%     end
%     
%     if Thickness_2015margin <0.5 &&Phase_2 == 1
%         T_2015margin = t;
%         Retreat_time = T_2015margin - T_steadystate
%         break
%     end

  %% Advance from 2015 margin SS to Trimline %%
%     if t >200 && abs(dice_rate(numTimeSteps)) < 0.01*max(max(dice_rate(1:numTimeSteps))) && Phase_1 == 1;
%         T_steadystate = t;
%         Phase_1 = 0;
%         Phase_2 = 1;
%         Tsealevel = Tsealevel-0.22;
% 
%     end
%    
%    if Phase_2 ==1 && Thickness_trimline >0.01
%         T_trimline = t;
%         Phase_2 = 0;
%         Phase_3 = 1;
%         SS_Advance_time = T_trimline-T_steadystate
%         break
%    end
   
   %% Simulated actual advance %%
%       if t >1
%           disp(['Ice Thickness at 2ka location:',num2str(Thickness_2ka)]);
%        break
%       end
   
%        if t >200 && abs(dice_rate(numTimeSteps)) < 0.01*max(max(dice_rate(1:numTimeSteps))) && Phase_1 == 1;
if t> 1 && Phase_1 ==1
    T_begin = t;
        Phase_1 = 0;
        Phase_2 = 1;
        preceeding_temp = Tsealevel;
        Tsealevel = T_initial_advance; %cooling for first ice advance
end

if Thickness_2ka > 0.05 && Phase_2 ==1
      T_2ka = t;
      Phase_2 = 0;
      Phase_3 = 1;
      filename_2ka = ['Fig_2ka_config_', num2str(maxBz),'_',num2str(Rmelt_rate),'.png'];
      saveas(figure(2), filename_2ka)
end
  
  if Thickness_2015margin > 0.05 && Phase_3 ==1
      T_2015margin = t;
      Tsealevel = T_transect_advance; %cooling for second ice advance
      Phase_3 = 0;
      Phase_4 = 1;
      filename_1000 = ['Fig_1ka_Config_', num2str(maxBz),'_',num2str(Rmelt_rate),'.png'];
      saveas(figure(2), filename_1000)
  end
%   
  if Thickness_trimline >0.01 && Phase_4 ==1 
        T_trimline = t;
        Phase_4 = 0;
        Phase_5 = 1;
        T_trimlineend = t;
        filename_lia = ['Fig_LIA_Config_', num2str(maxBz),'_',num2str(Rmelt_rate),'.png'];
        saveas(figure(2), filename_lia)
  end
  
  if Phase_5 == 1 && Thickness_2015margin > 0.1
% %       Tsealevel = Tsealevel +(0.0867*dt); %Qik 1995-2009
% %       Tsealevel = Tsealevel +(0.0141*dt); %Dewar Lakes 1959-2015
      Tsealevel = Tsealevel + (warming_rate * dt); %model calibrated
      Phase_6 = 1;
  end
      
  if Phase_6 ==1 && Thickness_2015margin < 0.1
          T_modern = t;
          Modern_2ka_thickness = Thickness_2ka;
          Temp_modern = Tsealevel;
          Phase_6 = 0;
          Phase_7 = 1;
          Modern_2ka_thickness = Thickness_2ka;
          filename_modern = ['Fig_modern_Config_', num2str(maxBz),'_',num2str(Rmelt_rate),'.png'];
          saveas(figure(2), filename_modern)
  end
  
  if Phase_7 == 1 && iceVolume >100
            Tsealevel = Tsealevel + (warming_rate * dt); %model calibrated
  end
      
        
  if Phase_7 == 1 && iceVolume <100
           T_gone = t;
           T_2ka;
           Initial_advance_time = T_2015margin-T_2ka;
           Transect_Advance_Time = T_trimline - T_2015margin;
           Occupation_time = T_trimlineend-T_trimline;
           Retreat_time = T_modern-T_trimlineend;
           Total_time = T_modern-T_2ka;
           Ice_Dissappearance = T_gone - T_modern;
       break
  end
      
    
    end %- TIME LOOP -%
    

% if Advance_time> 610
%     Tsealevel = Tsealevel - 0.05
% elseif Advance_time <590
%     Tsealevel = Tsealevel + 0.05
% else
%     break
% end
% end
%     if Thickness_trimline > 1;
%         Tsealevel = Tsealevel+0.01;
%         t = 0 ;
%     elseif Thickness_trimline <=0.01;
%         Tsealevel = Tsealevel-0.01
%         t = 0 ;
%     else
%         disp('You Got It');
%     end
%     end
%     
    diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ONE LAST REPORT AND PLOT BEFORE EXIT

if ( REPORT_TOGGLE )
    disp( sprintf('\nelapsed time: %1.2f s', toc ) ) ;
    disp( sprintf('simulation time: %1.2f yr', t) ) ;
    disp( sprintf('timestep: %1.2e yr', dt) ) ;
    %disp( sprintf('ELA: %1.2e m', ELA) ) ;

    %% fractional ice conservation
    iceVolume = sum(sum(H.*(H>0)))*dx*dy ;
    disp( sprintf('total ice: %1.2e km^3', iceVolume/1e9 ) ) ;
    %disp( sprintf('excess ice: %1.2f m^3', (( iceVolume - conserveIce*dx*dy )) ) ) ;
    %disp( sprintf('ice conservation (%%): %1.15f', 100 - 100*( iceVolume - conserveIce*dx*dy ) / iceVolume ) ) ;
    disp( sprintf('max ice thickness: %1.2e km', max(max(H))/1000 ) ) ;
    
    disp( sprintf('Preceeding Temp: %1.2f C', preceeding_temp ) ) ;
    disp( sprintf('Initial Advance Temp: %1.2f C', T_initial_advance ) ) ;
    disp( sprintf('Transect Advance Temp: %1.2f C', T_transect_advance ) ) ;
    disp( sprintf('Warming Rate: %1.2f C/yr', warming_rate ) ) 
    disp( sprintf('Total Warming: %1.2f C', warming_rate*Retreat_time ) )     
    disp( sprintf('Modern Temp: %1.3f C', Temp_modern ) ) ;
    disp( sprintf('Initial advance time: %1.2f yr', Initial_advance_time ) ) ;
    disp( sprintf('Transect_Advance_Time: %1.2f yr', Transect_Advance_Time ) ) ;
    disp( sprintf('Occupation_time: %1.2f yr', Occupation_time ) ) ;
    disp( sprintf('Retreat_time: %1.2f yr', Retreat_time ) ) ;
    disp( sprintf('Modeled Modern 2ka Ice Thickness: %1.2f m', Modern_2ka_thickness ) ) ;
    disp( sprintf('Total_time: %1.2f yr', Total_time ) ) ;
    disp( sprintf('Ice_Dissappearance: %1.2f yr', Ice_Dissappearance ) ) ;
    
end

beep
