%% Simulate Liquid Engine Performance
% Jason Chen, Project Caelus 501(c)(3) (https://projectcaelus.org)
% Heavily based on Nick Gloria's PropSim (github.com/ngloria1/PropSim)
% Driver for running PerformaceCode.m

%% Preprocessing
function [] = PropSimIntegrated()
    %% Unit Conversion
    psi_to_Pa = 6894.75729; % 1 psi in Pa
    in_to_m = 0.0254; % 1 in in m
    mm_to_m = 1e-3; % 1 mm in m
    lbf_to_N = 4.44822162; % 1 lbf in N
    lbm_to_kg = 0.453592; % 1 lbm in kg
    atm_to_Pa = 101325; % 1 atm in Pa
    L_to_m3 = 1e-3; % 1 L in m^3

    %% Options
    % 1: simulate combustion (hot fire), 0: no combustion (cold flow)
    mode.combustion_on = 1;
    % 1: simulate flight conditions (i.e. acceleration head), 0: ground test
    mode.flight_on = 0;

    %-------Gases-----------------------
    helium = Gas();
    helium.c_v = 3.12; % J/kg*K
    helium.molecular_mass = 4.0026e-3; % kg/mol

    nitrogen = Gas();
    nitrogen.c_v = 0.743e3; % J/kg*K
    nitrogen.molecular_mass = 2*14.0067e-3; % kg/mol

    %-------Injector Properties----------
    % Injector Exit Area
    inputs.ox.injector_area = 2.571e-05; % m^2
    inputs.fuel.injector_area = 6.545e-06; % 4.571e-06 m^2

    % Fuel line valve Cv values (flow coefficients), assumed in series
    inputs.fuel.valve_cvs = [0.9 0.7]; % Leave array empty/zero if no propellant valves
    inputs.ox.valve_cvs = [0.9 0.7];

    % Ball Valve Time to Injector Area (s)
    inputs.dt_valve_open = 0.01;

    % Discharge Coefficient
    inputs.ox.Cd_injector = 0.9;
    inputs.fuel.Cd_injector = 0.88;

    %-------Rocket Properties--------
    % Rocket Dry Mass
    inputs.mass_dry_rocket = 30*lbm_to_kg;

    %-------Oxidizer Properties--------
    % Tank Volume
    inputs.ox.V_tank = 11.564*L_to_m3; 

    % Nitrous Volume
    inputs.ox.V_l = 4.208*L_to_m3; % 3.59

    % Tank Inner Diameter
    inputs.ox.tank_id = 5*in_to_m;

    % Distance from Bottom of Tank to Injector
    inputs.ox.h_offset_tank = 2; % 0 m

    % Main Flow Line Diameter
    inputs.ox.d_flowline = .5*in_to_m;

    % Tank Temperature (K)
    inputs.ox.T_tank = 300;

    %-------Oxidizer Pressurant Properties--------
    inputs.ox_pressurant = Pressurant('oxidizer');
    inputs.ox_pressurant.gas_properties = helium;
    inputs.ox_pressurant.set_pressure = 750*psi_to_Pa;
    inputs.ox_pressurant.storage_initial_pressure = 0*psi_to_Pa;
    inputs.ox_pressurant.tank_volume = 0*L_to_m3;
    inputs.ox_pressurant.flow_CdA = 8*mm_to_m^2;

    % Are You Supercharging? (0 for 'NO' 1 for 'YES')
    inputs.ox_pressurant.active = 0;

    %-------Fuel Properties--------
    % Tank Volume
    inputs.fuel.V_tank = 11.564*L_to_m3; 

    % Fuel Volume
    inputs.fuel.V_l = 4*L_to_m3; % 0.881

    % Tank Inner Diameter
    inputs.fuel.tank_id = 5*in_to_m;

    % Distance from Bottom of Tank to Injector
    inputs.fuel.h_offset_tank = 24*in_to_m; % 24*in_to_m

    % Main Flow Line Diameter(in)
    inputs.fuel.d_flowline = .5*in_to_m;

    inputs.fuel.rho = 789; %Kg/m^3

    %-------Fuel Pressurant Properties--------
    inputs.fuel_pressurant = Pressurant('fuel');
    inputs.fuel_pressurant.gas_properties = nitrogen;
    inputs.fuel_pressurant.set_pressure = 750*psi_to_Pa; % 326 psi
    inputs.fuel_pressurant.storage_initial_pressure = 0*psi_to_Pa;
    inputs.fuel_pressurant.tank_volume = 0.0*L_to_m3;
    inputs.fuel_pressurant.flow_CdA = 8*mm_to_m^2;

    % Are You Supercharging? (0 for 'NO' 1 for 'YES')
    inputs.fuel_pressurant.active = 1;

    %-------Other Properties--------
    % Combustion chamber dimensions
    inputs.length_cc = 0.258; % m
    inputs.d_cc = 0.08; % m

    % Estimated nozzle efficiency
    inputs.nozzle_efficiency = 0.95;
    inputs.nozzle_correction_factor = 0.9830;

    % Estimated combustion efficiency
    inputs.c_star_efficiency = 0.85;

    % Nozzle Throat diameter
    inputs.d_throat = 30.46e-3;

    % Expansion Ratio
    inputs.exp_ratio = 3.26;

    % Ambient Temperature
    inputs.T_amb = 292; % K

    % Ambient Pressure
    inputs.p_amb = 9.554e04; % Pa

    % Load Combustion Data
    inputs.comb_data = load('CombustionData_T1_N2O.mat'); 
    inputs.comb_data = inputs.comb_data.CombData;

    %-------Other Options--------
    options.t_final   =  60;   % Integration time limit
    options.dt        = 0.01;  % Timestep [s]
    options.output_on = true;  % Whether or not to print/plot

    %% Run Performance Code
    PerformanceCode(inputs, mode, options);
end