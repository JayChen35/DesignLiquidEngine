%% Simulate Liquid Engine Performance
% Jason Chen, Project Caelus 501(c)(3) (https://projectcaelus.org)
% Heavily based on Nick Gloria's PropSim (github.com/ngloria1/PropSim)
% Driver for running PerformaceCode.m

function [] = PropSimIntegrated(json_path)
    %% Ingest PropSimInput (JSON)
    data = jsondecode(fileread(json_path));
    
    %% Unit Conversion
    psi_to_Pa = 6894.75729; % 1 psi in Pa
    in_to_m = 0.0254; % 1 in in m
    mm_to_m = 1e-3; % 1 mm in m
    lbf_to_N = 4.44822162; % 1 lbf in N
    lbm_to_kg = 0.453592; % 1 lbm in kg
    atm_to_Pa = 101325; % 1 atm in Pa
    L_to_m3 = 1e-3; % 1 L in m^3


    %% Options
    % Import tests data and plot against simulation data
    test_data.test_plots_on = data.test_data.test_plots_on; 
    % File from which to import test data
    test_data.test_data_file = data.test_data.test_data_file; 
    % Time offset of test data wrt simulation data [s]
    test_data.t_offset = data.test_data.t_offset; 

    % 1: Simulate combustion (hot fire), 0: no combustion (cold flow)
    mode.combustion_on = data.mode.combustion_on;
    % 1: Simulate flight conditions (i.e. acceleration head), 0: ground test
    mode.flight_on = data.mode.flight_on;

    %% Input Parameters
    %-------------Gases------------------

    helium = Gas();
    helium.c_v = 3.12*1e03; % J/kg*K
    helium.molecular_mass = 4.0026e-3; % kg/mol

    nitrogen = Gas();
    nitrogen.c_v = 0.743e3; % J/kg*K
    nitrogen.molecular_mass = 2*14.0067e-3; % kg/mol

    %-------Injector Properties----------

    %Injector Exit Area
    inputs.ox.injector_area = data.ox.injector_area; % m^2
    inputs.fuel.injector_area = data.fuel.injector_area; % m^2

    % Fuel line valve Cv values (flow coefficients), assumed in series
    % Leave array as empty or zero if there are no propellant valves
    inputs.fuel.valve_cvs = data.fuel.valve_cvs; 
    inputs.ox.valve_cvs = data.ox.valve_cvs;

    % Ball Valve Time to Injector Area (s)
    inputs.dt_valve_open = data.dt_valve_open;

    % Discharge Coefficient
    inputs.ox.Cd_injector = data.ox.Cd_injector;
    inputs.fuel.Cd_injector = data.fuel.Cd_injector;

    %-------Rocket Properties--------
    % Rocket Dry Mass
    inputs.mass_dry_rocket = data.mass_dry_rocket*lbm_to_kg;

    %-------Oxidizer Properties--------
    % Tank Volume
    inputs.ox.V_tank = data.ox.V_tank*L_to_m3; 

    % Nitrous Volume
    inputs.ox.V_l = data.ox.V_l*L_to_m3; % 3.59

    % Tank Inner Diameter
    inputs.ox.tank_id = data.ox.tank_id*in_to_m;

    % Distance from Bottom of Tank to Injector
    inputs.ox.h_offset_tank = data.ox.h_offset_tank*in_to_m;

    % Main Flow Line Diameter
    inputs.ox.d_flowline = data.ox.d_flowline*in_to_m;

    % Tank Temperature (K)
    inputs.ox.T_tank = data.ox.T_tank;

    %-------Oxidizer Pressurant Properties--------

    inputs.ox_pressurant = Pressurant('oxidizer');
    if data.ox_pressurant.gas_properties == 0
        inputs.ox_pressurant.gas_properties = nitrogen;
    else
        inputs.ox_pressurant.gas_properties = helium;
    end
    inputs.ox_pressurant.set_pressure = ...
        data.ox_pressurant.set_pressure*psi_to_Pa;
    inputs.ox_pressurant.storage_initial_pressure = ...
        data.ox_pressurant.storage_initial_pressure*psi_to_Pa;
    inputs.ox_pressurant.tank_volume = ...
        data.ox_pressurant.tank_volume*L_to_m3;
    inputs.ox_pressurant.flow_CdA = data.ox_pressurant.flow_CdA; 
    % inputs.ox_pressurant.flow_CdA = 8*mm_to_m^2;

    %Are You Supercharging? (0 for 'NO' 1 for 'YES')
    inputs.ox_pressurant.active = data.ox_pressurant.active;

    %-------Fuel Properties--------
    % Tank Volume
    inputs.fuel.V_tank = data.fuel.V_tank*L_to_m3; 

    % Fuel Volume
    inputs.fuel.V_l = data.fuel.V_l*L_to_m3;

    % Tank Inner Diameter
    inputs.fuel.tank_id = data.fuel.tank_id*in_to_m;

    % Distance from Bottom of Tank to Injector
    inputs.fuel.h_offset_tank = data.fuel.h_offset_tank*in_to_m;

    % Main Flow Line Diameter(in)
    inputs.fuel.d_flowline = data.fuel.d_flowline*in_to_m;

    inputs.fuel.rho = data.rho_f; %Kg/m^3

    %-------Fuel Pressurant Properties--------

    inputs.fuel_pressurant = Pressurant('fuel');
    if data.fuel_pressurant.gas_properties == 0
        inputs.fuel_pressurant.gas_properties = nitrogen;
    else
        inputs.fuel_pressurant.gas_properties = helium;
    end
    inputs.fuel_pressurant.set_pressure = ...
        data.fuel_pressurant.set_pressure*psi_to_Pa;
    inputs.fuel_pressurant.storage_initial_pressure = ...
        data.fuel_pressurant.storage_initial_pressure*psi_to_Pa;
    inputs.fuel_pressurant.tank_volume = ...
        data.fuel_pressurant.tank_volume*L_to_m3;
    inputs.fuel_pressurant.flow_CdA = data.fuel_pressurant.flow_CdA; 

    % Are You Supercharging? (0 for 'NO' 1 for 'YES')
    inputs.fuel_pressurant.active = data.fuel_pressurant.active;

    %-------Other Properties--------
    % Combustion chamber dimensions
    inputs.length_cc = data.length_cc; % m
    inputs.d_cc = data.d_cc; % m

    % Estimated nozzle efficiency
    inputs.nozzle_efficiency = data.nozzle_efficiency;
    inputs.nozzle_correction_factor = data.nozzle_correction_factor;

    % Estimated combustion efficiency
    inputs.c_star_efficiency = data.c_star_efficiency;

    % Nozzle Throat diameter
    inputs.d_throat = data.d_throat;

    % Expansion Ratio
    inputs.exp_ratio = data.exp_ratio;

    % Ambient Temperature (K)
    inputs.T_amb = data.T_amb;

    % Ambient Pressure (Pa)
    inputs.p_amb = data.P3;

    % Load Combustion Data
    inputs.comb_data = load('CombustionData_T1_N2O.mat'); 
    inputs.comb_data = inputs.comb_data.CombData;

    %-------Other Options--------
    options.t_final   =  60;   % Integration time limit
    options.dt        = 0.01;  % Timestep [s]
    options.output_on = true;  % Whether or not to print/plot

    %% Run Performance Code
    PerformanceCode(inputs, mode, test_data, options);
end
    