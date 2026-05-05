%% Lunar Mission Architecture Trade Study: Delta-V & Mass Analysis

clear all; clc; close all;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Constants

g0 = 9.80665 / 1000; % Standard gravity (km/s^2)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Load Data Files from the Patched Conic Runs
%...-- Ensure these match the exact 'Run_Name' strings --

Scenario1 = load('Export_HEO_to_Equatorial.mat');
Scenario2 = load('Export_HEO_to_NRHO.mat');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Define the Mass Sweep Space (IMLEO from 1,000 kg to 30,000 kg)
%...Represents testing different spacecraft sizes (e.g., small probe vs. Orion capsule)

mass_sweep = linspace(1000, 30000, 500); 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Required Propellant Mass for Scenario 1 (Equatorial LLO)

prop_TLI_1 = mass_sweep .* (1 - exp(-Scenario1.DV_TLI / (Scenario1.Isp * g0)));
mass_after_TLI_1 = mass_sweep - prop_TLI_1;
prop_LOI_1 = mass_after_TLI_1 .* (1 - exp(-Scenario1.DV_LOI / (Scenario1.Isp * g0)));
Total_Propellant_1 = prop_TLI_1 + prop_LOI_1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Calculate Required Propellant Mass for Scenario 2 (NRHO)

prop_TLI_2 = mass_sweep .* (1 - exp(-Scenario2.DV_TLI / (Scenario2.Isp * g0)));
mass_after_TLI_2 = mass_sweep - prop_TLI_2;
prop_LOI_2 = mass_after_TLI_2 .* (1 - exp(-Scenario2.DV_LOI / (Scenario2.Isp * g0)));
Total_Propellant_2 = prop_TLI_2 + prop_LOI_2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Generate Trade Study Visualizations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure('Position', [100, 100, 1200, 500], 'Name', 'Mission Architecture Comparison' ...
    , 'NumberTitle', 'off');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Plot 1: Total Propellant vs Initial Mass

subplot(1,2,1); hold on; grid on;
plot(mass_sweep, Total_Propellant_1, 'b-', 'LineWidth', 2, 'DisplayName' ...
    , strrep(Scenario1.Name, '_', ' '));
plot(mass_sweep, Total_Propellant_2, 'r--', 'LineWidth', 2, 'DisplayName' ...
    , strrep(Scenario2.Name, '_', ' '));
title('Total Propellant Required vs. Initial Mass in LEO');
xlabel('Initial Spacecraft Mass in LEO (kg)');
ylabel('Total Propellant Required (kg)');
legend('Location', 'northwest');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%...Plot 2: Delta-V Budget Bar Chart

subplot(1,2,2); hold on; grid on;
categories = categorical({'TLI Burn', 'LOI Burn', 'Total Mission'});
categories = reordercats(categories, {'TLI Burn', 'LOI Burn', 'Total Mission'});
dv_data = [Scenario1.DV_TLI, Scenario2.DV_TLI; 
           Scenario1.DV_LOI, Scenario2.DV_LOI; 
           Scenario1.DV_Total, Scenario2.DV_Total];

bar_handle = bar(categories, dv_data);
bar_handle(1).FaceColor = 'b';
bar_handle(2).FaceColor = 'r';
title('Mission \DeltaV Cost Comparison');
ylabel('\DeltaV (km/s)');
legend(strrep(Scenario1.Name, '_', ' '), ...
    strrep(Scenario2.Name, '_', ' '), 'Location', 'northwest');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Created by Cameron Edwards, Florida Institute of Technology (April - May 2026)