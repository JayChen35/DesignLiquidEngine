function record = PerformanceCode(inputs, mode, options)
% Take inputs defining motor characteristics, take test data, and calculate
% and plot performance. 
%   INPUTS:
%       - inputs: see Integration.m
%       - mode: see Integration.m
%       - test_data: structure of data related to test data for plotting
%       purposes
%           - test_plots_on: 1 for plot test data, 0 for nothing
%           - test_data_file: filename of test data within "Test Data"
%           directory
%               - All data should be in physical values (not voltages) in
%               SI units of m, kg, K, s, Pa, N, etc.
%           - t_offset: time at which ignition occurs in test data]
%       - options: set of options
%           - t_final: simulation time limit
%           - dt: output time resolution
%           - output_on: true for plots, false for no plots
%   OUTPUTS:
%       - Isp: specific impulse [s]
%       - time: time vector over burn for F_thrust
%       - F_thrust: vector of thrust over burn
%       - F_thrust_RASAERO: vector of thrust points used for RASAero thrust
%       curve output
%       - impulse: total impulse
%       - pressure_drop: average pressure drop over burn duration
%       - F_thrust_RASAERO.txt file in "Outputs" with thrust curve that can 
%       be adap_oxtanked for RASAero or open rocket. 

% To Do:
%   - Change to higher-order integration method (e.g. Runge-Kutta 4th order
%   or matlab ODE45 function
%   - Add feedline fricitional losses to two-phase injector flow calculator
%   - Investigate differences between simulation and test data for
%   supercharged oxidizer tank case

%% Constants
R_u = 8.3144621; % universal gas constant [J/mol*K]
g_0 = 9.80665; % standard gravitational constant [m/s^2]

% Unit Conversion
psi_to_Pa = 6894.75729; % 1 psi in Pa
in_to_m = 0.0254; % 1 in in m
lbf_to_N = 4.44822162; % 1 lbf in N
atm_to_Pa = 101325; % 1 atm in Pa

%% Integration Parameters
default_options.t_final  =  60;    % Integration time limit
default_options.dt      = 0.01;  % Timestep [s]
default_options.output_on = true;
if nargin < 4
    options = default_options;
else
    if ~isfield(options, 't_final')
        options.t_final = default_options.t_final;
    end
    if ~isfield(options, 'dt')
        options.dt = default_options.dt;
    end
    if ~isfield(options, 'output_on')
        options.output_on = default_options.output_on;
    end
end
tspan = 0:options.dt:options.t_final;

%% Calculate initial properties of the nitrous in the tank
N2O = N2O_Properties(inputs.ox.T_tank);

%Our integration variables are oxidizer mass and liquid oxidizer volume
Mox = N2O.rho_l*(inputs.ox.V_l) + N2O.rho_g*(inputs.ox.V_tank - inputs.ox.V_l);
if options.output_on
    fprintf('Initial oxidizer mass: %.2f kg\n', Mox);
end

%% Create Injector Area Function
function [F] = Throttle(t)
    F = zeros(length(t));
    for i = 1:length(t)
        if t(i) > inputs.dt_valve_open
            F(i) = 1;
        else
            F(i) = (t(i)/inputs.dt_valve_open);
        end
    end
end
inputs.Throttle = @Throttle;

%% Run Integration Function
tic;
[time, record] = Integration(inputs,mode,tspan);
if options.output_on
    toc;
end

F_thrust = record.F_thrust;
p_cc = record.p_cc;
p_oxtank = record.p_oxtank;
p_oxpresstank = record.p_oxpresstank;
p_fueltank = record.p_fueltank;
p_fuelpresstank = record.p_fuelpresstank;
p_oxmanifold = record.p_oxmanifold;
T_oxtank = record.T_oxtank;
T_cc = record.T_cc;
area_core = record.area_core;
OF = record.OF_i;
gamma_ex = record.gamma_ex;
m_dot_ox = record.m_dot_ox;
m_dot_fuel = record.m_dot_fuel;
p_crit = record.p_crit;
m_dot_ox_crit = record.m_dot_ox_crit;
M_e = record.M_e;
p_exit = record.p_exit;
p_shock = record.p_shock;

% Use trapezoidal integration
impulse = trapz(time, F_thrust);
Mox_initial = record.m_ox(1);
Mox = record.m_ox(1) - record.m_ox(end);

Mfuel_initial = record.m_fuel(1);
Mfuel = record.m_fuel(1) - record.m_fuel(end);

fprintf('Pressurant Mass: %.3f kg\n', record.m_press)
fprintf(['Impulse: %.2f kN*s\t\t\nOxidizer Mass Spent: '...
    '%.2f kg\t\tOxidizer Mass Remaining: %.2f kg\nFuel Mass Spent: %.2f kg\t\t' ...
    'Fuel Mass Remaining: %.2f kg\n OF ratio: %.2f \n']...
    , impulse/1000, Mox, Mox_initial-Mox, Mfuel, Mfuel_initial - Mfuel, ...
    Mox/Mfuel);
fprintf('Isp: %.1f s\t\tC*: %.0f m/s\t\tC_f: %.2f\n', ...
    record.Isp/g_0, record.c_star, record.c_f)
save('PropSimOutput.mat', 'record')
end
