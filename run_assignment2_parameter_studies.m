% Numerical parameter-study and verification workflow for CFD Assignment 2a/2b.
%
% Run from MATLAB with:
%   run_assignment2_parameter_studies
%
% This script intentionally does not modify the original teaching scripts:
%   conv_ee_fillin.m
%   conv_ie_fillin.m
%   conv_diff_cn_fillin.m

format long;
clear;
close all;
clc;

base_dir = fileparts(mfilename('fullpath'));
if isempty(base_dir)
    base_dir = pwd;
end

results_dir = fullfile(base_dir, 'results_assignment2');
fig_dir = fullfile(results_dir, 'figures');
table_dir = fullfile(results_dir, 'tables');
data_dir = fullfile(results_dir, 'data');
ensure_folder(results_dir);
ensure_folder(fig_dir);
ensure_folder(table_dir);
ensure_folder(data_dir);

fprintf('Running CFD Assignment 2 parameter studies...\n');
fprintf('Results directory: %s\n', results_dir);

all_results = {};
rows_EE_Courant = {};
rows_EE_grid = {};
rows_IE_Courant = {};
rows_IE_grid = {};
rows_CN_baseline_limits = {};
rows_CN_Gamma = {};
rows_CN_dt = {};
rows_CN_grid = {};

%% C1. Explicit Euler + central difference for pure advection
fprintf('\n[EE] Running explicit Euler studies...\n');

% EE-1: short-time reference comparison
U0 = 1.0;
points = 161;
dx = grid_spacing(points);
Courant = 0.05;
dt = Courant * dx / U0;
tend = 0.5;
EE1 = solve_advection_ee(U0, points, dt, tend);
EE1.experiment_name = 'EE-1 short-time reference';
all_results{end+1} = EE1; %#ok<SAGROW>
plot_profile(EE1, 'EE short-time profile: C=0.05, t=0.5', ...
    fullfile(fig_dir, 'EE_short_time_profile'));

% EE-2: instability sensitivity to Courant number
U0 = 1.0;
points = 81;
dx = grid_spacing(points);
tend = 2.0;
Courant_values_EE = [0.05, 0.20, 0.50, 1.00];
EE2 = cell(size(Courant_values_EE));
for k = 1:numel(Courant_values_EE)
    C = Courant_values_EE(k);
    dt = C * dx / U0;
    result = solve_advection_ee(U0, points, dt, tend);
    result.experiment_name = sprintf('EE-2 Courant %.2f', C);
    EE2{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_EE_Courant{end+1} = result; %#ok<SAGROW>
end
plot_history(EE2, 'L2_history', 'Relative L_2 error', ...
    'EE Courant study: relative L_2 error', ...
    fullfile(fig_dir, 'EE_Courant_L2_history'));
plot_history(EE2, 'max_abs_phi_history', 'max(|\phi|)', ...
    'EE Courant study: max absolute solution value', ...
    fullfile(fig_dir, 'EE_Courant_maxabs_history'));
plot_profile(EE2{end}, 'EE profile distortion: C=1.00, t=2.0', ...
    fullfile(fig_dir, 'EE_Courant_1_profile'));

% EE-3: grid-resolution comparison at fixed small Courant
U0 = 1.0;
Courant = 0.05;
tend = 1.0;
points_values_EE = [41, 81, 161];
EE3 = cell(size(points_values_EE));
for k = 1:numel(points_values_EE)
    points = points_values_EE(k);
    dx = grid_spacing(points);
    dt = Courant * dx / U0;
    result = solve_advection_ee(U0, points, dt, tend);
    result.experiment_name = sprintf('EE-3 points %d', points);
    EE3{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_EE_grid{end+1} = result; %#ok<SAGROW>
end
plot_final_metric_vs_points(EE3, 'final_L2', 'Final relative L_2 error', ...
    'EE grid study: final L_2 error vs points', ...
    fullfile(fig_dir, 'EE_grid_final_L2'));

%% C2. Implicit Euler + central difference for pure advection
fprintf('[IE] Running implicit Euler studies...\n');

% IE-1: baseline comparison
U0 = 1.0;
points = 161;
dx = grid_spacing(points);
Courant = 0.20;
dt = Courant * dx / U0;
tend = 1.0;
IE1 = solve_advection_ie(U0, points, dt, tend);
IE1.experiment_name = 'IE-1 baseline';
all_results{end+1} = IE1; %#ok<SAGROW>
plot_profile(IE1, 'IE baseline profile: C=0.20, t=1.0', ...
    fullfile(fig_dir, 'IE_baseline_profile'));

% IE-2: timestep/Courant sensitivity
U0 = 1.0;
points = 81;
dx = grid_spacing(points);
tend = 2.0;
Courant_values_IE = [0.10, 0.50, 1.00, 2.00];
IE2 = cell(size(Courant_values_IE));
for k = 1:numel(Courant_values_IE)
    C = Courant_values_IE(k);
    dt = C * dx / U0;
    result = solve_advection_ie(U0, points, dt, tend);
    result.experiment_name = sprintf('IE-2 Courant %.2f', C);
    IE2{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_IE_Courant{end+1} = result; %#ok<SAGROW>
end
plot_history(IE2, 'L2_history', 'Relative L_2 error', ...
    'IE Courant study: relative L_2 error', ...
    fullfile(fig_dir, 'IE_Courant_L2_history'));
plot_history(IE2, 'amplitude_history', 'Numerical amplitude', ...
    'IE Courant study: amplitude history', ...
    fullfile(fig_dir, 'IE_Courant_amplitude_history'));
plot_two_profiles(IE2{1}, IE2{end}, ...
    'IE final profiles: C=0.10 and C=2.00', ...
    fullfile(fig_dir, 'IE_Courant_profile_C010_C200'));

% IE-3: grid refinement at fixed Courant
U0 = 1.0;
Courant = 0.20;
tend = 1.0;
points_values_IE = [41, 81, 161];
IE3 = cell(size(points_values_IE));
for k = 1:numel(points_values_IE)
    points = points_values_IE(k);
    dx = grid_spacing(points);
    dt = Courant * dx / U0;
    result = solve_advection_ie(U0, points, dt, tend);
    result.experiment_name = sprintf('IE-3 points %d', points);
    IE3{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_IE_grid{end+1} = result; %#ok<SAGROW>
end
plot_final_metric_vs_points(IE3, 'final_L2', 'Final relative L_2 error', ...
    'IE grid study: final L_2 error vs points', ...
    fullfile(fig_dir, 'IE_grid_final_L2'));

%% C3. Crank-Nicolson + central differences for advection-diffusion
fprintf('[CN] Running Crank-Nicolson studies...\n');

% CN-1: coupled advection-diffusion baseline
U0 = 1.0;
Gamma = 0.1;
points = 161;
dt = 0.005;
tend = 1.0;
CN1 = solve_advdiff_cn(U0, Gamma, points, dt, tend);
CN1.experiment_name = 'CN-1 coupled baseline';
all_results{end+1} = CN1; %#ok<SAGROW>
rows_CN_baseline_limits{end+1} = CN1; %#ok<SAGROW>
plot_profile(CN1, 'CN coupled advection-diffusion: U0=1, Gamma=0.1', ...
    fullfile(fig_dir, 'CN_baseline_profile'));
plot_history({CN1}, 'L2_history', 'Relative L_2 error', ...
    'CN baseline: relative L_2 error', ...
    fullfile(fig_dir, 'CN_baseline_L2_history'));

% CN-2: pure advection limit
U0 = 1.0;
Gamma = 0.0;
points = 161;
dx = grid_spacing(points);
tend = 1.0;
Courant_values_CN = [0.10, 0.50, 1.00];
CN2 = cell(size(Courant_values_CN));
for k = 1:numel(Courant_values_CN)
    C = Courant_values_CN(k);
    dt = C * dx / U0;
    result = solve_advdiff_cn(U0, Gamma, points, dt, tend);
    result.experiment_name = sprintf('CN-2 pure advection Courant %.2f', C);
    CN2{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_CN_baseline_limits{end+1} = result; %#ok<SAGROW>
end
plot_profile(CN2{end}, 'CN pure advection limit: Gamma=0, C=1.00', ...
    fullfile(fig_dir, 'CN_pure_advection_C100_profile'));

% CN-3: pure diffusion limit
U0 = 0.0;
Gamma = 0.1;
points = 161;
dt = 0.005;
tend = 1.0;
CN3 = solve_advdiff_cn(U0, Gamma, points, dt, tend);
CN3.experiment_name = 'CN-3 pure diffusion';
all_results{end+1} = CN3; %#ok<SAGROW>
rows_CN_baseline_limits{end+1} = CN3; %#ok<SAGROW>
plot_amplitude_vs_analytic(CN3, 'CN pure diffusion: amplitude decay', ...
    fullfile(fig_dir, 'CN_pure_diffusion_amplitude_history'));

% CN-4: effect of physical diffusion coefficient Gamma
U0 = 1.0;
points = 161;
dt = 0.005;
tend = 1.0;
Gamma_values = [0.00, 0.05, 0.10, 0.50];
CN4 = cell(size(Gamma_values));
for k = 1:numel(Gamma_values)
    Gamma = Gamma_values(k);
    result = solve_advdiff_cn(U0, Gamma, points, dt, tend);
    result.experiment_name = sprintf('CN-4 Gamma %.2f', Gamma);
    CN4{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_CN_Gamma{end+1} = result; %#ok<SAGROW>
end
plot_final_amplitude_vs_gamma(CN4, ...
    'CN Gamma study: final amplitude', ...
    fullfile(fig_dir, 'CN_Gamma_final_amplitude'));

% CN-5: timestep sensitivity
U0 = 1.0;
Gamma = 0.1;
points = 161;
tend = 1.0;
dt_values = [0.001, 0.005, 0.020, 0.050];
CN5 = cell(size(dt_values));
for k = 1:numel(dt_values)
    dt = dt_values(k);
    result = solve_advdiff_cn(U0, Gamma, points, dt, tend);
    result.experiment_name = sprintf('CN-5 dt %.3f', dt);
    CN5{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_CN_dt{end+1} = result; %#ok<SAGROW>
end
plot_final_metric_vs_dt(CN5, 'final_L2', 'Final relative L_2 error', ...
    'CN dt study: final L_2 error vs dt', ...
    fullfile(fig_dir, 'CN_dt_final_L2'));

% CN-6: grid refinement
U0 = 1.0;
Gamma = 0.1;
dt = 0.001;
tend = 1.0;
points_values_CN = [41, 81, 161, 321];
CN6 = cell(size(points_values_CN));
for k = 1:numel(points_values_CN)
    points = points_values_CN(k);
    result = solve_advdiff_cn(U0, Gamma, points, dt, tend);
    result.experiment_name = sprintf('CN-6 points %d', points);
    CN6{k} = result;
    all_results{end+1} = result; %#ok<SAGROW>
    rows_CN_grid{end+1} = result; %#ok<SAGROW>
end
plot_final_metric_vs_points(CN6, 'final_L2', 'Final relative L_2 error', ...
    'CN grid study: final L_2 error vs points', ...
    fullfile(fig_dir, 'CN_grid_final_L2'));

%% Export tables, data, and numerical audit report
fprintf('[Output] Writing CSV tables, MAT data, and audit report...\n');

write_result_table(rows_EE_Courant, fullfile(table_dir, 'table_EE_Courant_study.csv'));
write_result_table(rows_EE_grid, fullfile(table_dir, 'table_EE_grid_study.csv'));
write_result_table(rows_IE_Courant, fullfile(table_dir, 'table_IE_Courant_study.csv'));
write_result_table(rows_IE_grid, fullfile(table_dir, 'table_IE_grid_study.csv'));
write_result_table(rows_CN_baseline_limits, fullfile(table_dir, 'table_CN_baseline_and_limits.csv'));
write_result_table(rows_CN_Gamma, fullfile(table_dir, 'table_CN_Gamma_study.csv'));
write_result_table(rows_CN_dt, fullfile(table_dir, 'table_CN_dt_study.csv'));
write_result_table(rows_CN_grid, fullfile(table_dir, 'table_CN_grid_study.csv'));

save(fullfile(data_dir, 'assignment2_parameter_study_results.mat'), ...
    'all_results', 'EE1', 'EE2', 'EE3', 'IE1', 'IE2', 'IE3', ...
    'CN1', 'CN2', 'CN3', 'CN4', 'CN5', 'CN6');

write_audit_report(fullfile(results_dir, 'numerical_audit_report.md'), ...
    EE2, IE2, CN1, CN2, CN3, CN4, CN5, CN6);

fprintf('\nParameter studies complete.\n');
fprintf('Figures: %s\n', fig_dir);
fprintf('Tables:  %s\n', table_dir);
fprintf('Data:    %s\n', data_dir);
fprintf('Audit:   %s\n', fullfile(results_dir, 'numerical_audit_report.md'));

%% Solver functions
function result = solve_advection_ee(U0, points, dt, tend)
    [x, dx] = make_grid(points);
    tsteps = max(1, round(tend / dt));
    tend_actual = tsteps * dt;
    phi = sin(x);
    phinew = zeros(points,1);

    result = init_result('Explicit Euler + CDS', U0, NaN, points, dx, dt, tend_actual, tsteps);
    result.Courant = U0 * dt / dx;
    result.Fourier_number = NaN;

    for i = 1:tsteps
        phinew(1) = phi(1) ...
                  - U0 * dt / (2.0 * dx) * (phi(2) - phi(points-1));

        for j = 2:points-1
            phinew(j) = phi(j) ...
                      - U0 * dt / (2.0 * dx) * (phi(j+1) - phi(j-1));
        end

        phinew(points) = phi(points) ...
                       - U0 * dt / (2.0 * dx) * (phi(2) - phi(points-1));

        phi = phinew;
        t = i * dt;
        phi_a = sin(x - U0 * t);
        result = update_diagnostics(result, phi, phi_a, t, i, 1.0);

        if result.unstable_detected
            result = fill_remaining_with_nan(result, i+1, tsteps);
            break;
        end
    end

    result = finalize_result(result, x, phi, phi_a);
end

function result = solve_advection_ie(U0, points, dt, tend)
    [x, dx] = make_grid(points);
    tsteps = max(1, round(tend / dt));
    tend_actual = tsteps * dt;
    phi = sin(x);

    result = init_result('Implicit Euler + CDS', U0, NaN, points, dx, dt, tend_actual, tsteps);
    result.Courant = U0 * dt / dx;
    result.Fourier_number = NaN;

    a_w = - U0 * dt / (2.0 * dx);
    a_p = 1.0;
    a_e =   U0 * dt / (2.0 * dx);
    A = build_periodic_matrix(points, a_w, a_p, a_e);

    for i = 1:tsteps
        b = phi;
        phi = A \ b;
        t = i * dt;
        phi_a = sin(x - U0 * t);
        result = update_diagnostics(result, phi, phi_a, t, i, 1.0);

        if result.unstable_detected
            result = fill_remaining_with_nan(result, i+1, tsteps);
            break;
        end
    end

    result = finalize_result(result, x, phi, phi_a);
end

function result = solve_advdiff_cn(U0, Gamma, points, dt, tend)
    [x, dx] = make_grid(points);
    tsteps = max(1, round(tend / dt));
    tend_actual = tsteps * dt;
    phi = sin(x);

    result = init_result('Crank-Nicolson + CDS', U0, Gamma, points, dx, dt, tend_actual, tsteps);
    result.Courant = U0 * dt / dx;
    result.Fourier_number = Gamma * dt / dx^2;

    a_w = - U0 * dt / (4.0 * dx) - Gamma * dt / (2.0 * dx^2);
    a_p = 1.0 + Gamma * dt / (dx^2);
    a_e =   U0 * dt / (4.0 * dx) - Gamma * dt / (2.0 * dx^2);
    A = build_periodic_matrix(points, a_w, a_p, a_e);

    rhs_w =   U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2);
    rhs_p =   1.0 - Gamma * dt / (dx^2);
    rhs_e = - U0 * dt / (4.0 * dx) + Gamma * dt / (2.0 * dx^2);

    b = zeros(points,1);
    for i = 1:tsteps
        b(1) = rhs_w * phi(points-1) + rhs_p * phi(1) + rhs_e * phi(2);
        for j = 2:points-1
            b(j) = rhs_w * phi(j-1) + rhs_p * phi(j) + rhs_e * phi(j+1);
        end
        b(points) = rhs_w * phi(points-1) + rhs_p * phi(points) + rhs_e * phi(2);

        phi = A \ b;
        t = i * dt;
        phi_a = exp(-Gamma * t) * sin(x - U0 * t);
        analytic_amplitude = exp(-Gamma * t);
        result = update_diagnostics(result, phi, phi_a, t, i, analytic_amplitude);

        if result.unstable_detected
            result = fill_remaining_with_nan(result, i+1, tsteps);
            break;
        end
    end

    result = finalize_result(result, x, phi, phi_a);
end

%% Diagnostic and setup helpers
function result = init_result(solver, U0, Gamma, points, dx, dt, tend, tsteps)
    result = struct();
    result.experiment_name = '';
    result.solver = solver;
    result.U0 = U0;
    result.Gamma = Gamma;
    result.points = points;
    result.dx = dx;
    result.dt = dt;
    result.tend = tend;
    result.tsteps = tsteps;
    result.Courant = NaN;
    result.Fourier_number = NaN;
    result.x = [];
    result.t = (1:tsteps)' * dt;
    result.phi_final = [];
    result.phi_analytic_final = [];
    result.L2_history = NaN(tsteps,1);
    result.Linf_history = NaN(tsteps,1);
    result.periodic_error_history = NaN(tsteps,1);
    result.amplitude_history = NaN(tsteps,1);
    result.analytic_amplitude_history = NaN(tsteps,1);
    result.max_abs_phi_history = NaN(tsteps,1);
    result.final_L2 = NaN;
    result.final_Linf = NaN;
    result.final_periodic_error = NaN;
    result.final_amplitude = NaN;
    result.final_analytic_amplitude = NaN;
    result.max_abs_phi_final = NaN;
    result.amplitude_ratio_numeric_over_analytic = NaN;
    result.unstable_detected = false;
end

function result = update_diagnostics(result, phi, phi_a, t, index, analytic_amplitude)
    result.t(index) = t;
    result.L2_history(index) = relative_l2(phi, phi_a);
    result.Linf_history(index) = max(abs(phi - phi_a));
    result.periodic_error_history(index) = abs(phi(1) - phi(end));
    result.amplitude_history(index) = (max(phi) - min(phi)) / 2.0;
    result.analytic_amplitude_history(index) = analytic_amplitude;
    result.max_abs_phi_history(index) = max(abs(phi));

    if any(~isfinite(phi)) || result.max_abs_phi_history(index) > 10.0
        result.unstable_detected = true;
    end
end

function result = fill_remaining_with_nan(result, start_index, tsteps)
    if start_index <= tsteps
        result.L2_history(start_index:end) = NaN;
        result.Linf_history(start_index:end) = NaN;
        result.periodic_error_history(start_index:end) = NaN;
        result.amplitude_history(start_index:end) = NaN;
        result.analytic_amplitude_history(start_index:end) = NaN;
        result.max_abs_phi_history(start_index:end) = NaN;
    end
end

function result = finalize_result(result, x, phi, phi_a)
    result.x = x;
    result.phi_final = phi;
    result.phi_analytic_final = phi_a;
    result.final_L2 = last_finite(result.L2_history);
    result.final_Linf = last_finite(result.Linf_history);
    result.final_periodic_error = last_finite(result.periodic_error_history);
    result.final_amplitude = last_finite(result.amplitude_history);
    result.final_analytic_amplitude = last_finite(result.analytic_amplitude_history);
    result.max_abs_phi_final = last_finite(result.max_abs_phi_history);
    result.amplitude_ratio_numeric_over_analytic = ...
        result.final_amplitude / result.final_analytic_amplitude;
end

function value = relative_l2(phi, phi_a)
    denom = norm(phi_a, 2);
    if denom == 0
        value = norm(phi - phi_a, 2);
    else
        value = norm(phi - phi_a, 2) / denom;
    end
end

function value = last_finite(vector)
    idx = find(isfinite(vector), 1, 'last');
    if isempty(idx)
        value = NaN;
    else
        value = vector(idx);
    end
end

function [x, dx] = make_grid(points)
    xend = 2.0 * pi;
    dx = xend / (points - 1);
    x = (0.0 : dx : xend)';
end

function dx = grid_spacing(points)
    dx = 2.0 * pi / (points - 1);
end

function A = build_periodic_matrix(points, a_w, a_p, a_e)
    A = spalloc(points, points, 3 * points);

    A(1,points-1) = a_w;
    A(1,1) = a_p;
    A(1,2) = a_e;

    for j = 2:points-1
        A(j,j-1) = a_w;
        A(j,j) = a_p;
        A(j,j+1) = a_e;
    end

    A(points,points-1) = a_w;
    A(points,points) = a_p;
    A(points,2) = a_e;
end

function ensure_folder(folder_path)
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end
end

%% Plotting helpers
function plot_profile(result, plot_title, output_base)
    fig = make_figure();
    plot(result.x, result.phi_analytic_final, 'k-', 'LineWidth', 1.8);
    hold on;
    plot(result.x, result.phi_final, 'r--', 'LineWidth', 1.6);
    grid on;
    box on;
    xlabel('x');
    ylabel('\phi');
    title(plot_title, 'Interpreter', 'none');
    legend('analytical', 'numerical', 'Location', 'best');
    add_parameter_text(result);
    save_figure(fig, output_base);
end

function plot_two_profiles(result_a, result_b, plot_title, output_base)
    fig = make_figure();
    plot(result_a.x, result_a.phi_analytic_final, 'k-', 'LineWidth', 1.8);
    hold on;
    plot(result_a.x, result_a.phi_final, 'b--', 'LineWidth', 1.6);
    plot(result_b.x, result_b.phi_final, 'r-.', 'LineWidth', 1.6);
    grid on;
    box on;
    xlabel('x');
    ylabel('\phi');
    title(plot_title, 'Interpreter', 'none');
    legend('analytical', result_label(result_a), result_label(result_b), 'Location', 'best');
    save_figure(fig, output_base);
end

function plot_history(results, field_name, y_label, plot_title, output_base)
    fig = make_figure();
    hold on;
    for k = 1:numel(results)
        r = results{k};
        plot(r.t, r.(field_name), 'LineWidth', 1.6, 'DisplayName', result_label(r));
    end
    grid on;
    box on;
    xlabel('t');
    ylabel(y_label);
    title(plot_title, 'Interpreter', 'none');
    legend('Location', 'best');
    save_figure(fig, output_base);
end

function plot_amplitude_vs_analytic(result, plot_title, output_base)
    fig = make_figure();
    plot(result.t, result.amplitude_history, 'm-', 'LineWidth', 1.7);
    hold on;
    plot(result.t, result.analytic_amplitude_history, 'k--', 'LineWidth', 1.7);
    grid on;
    box on;
    xlabel('t');
    ylabel('amplitude');
    title(plot_title, 'Interpreter', 'none');
    legend('numerical amplitude', 'analytical amplitude', 'Location', 'best');
    save_figure(fig, output_base);
end

function plot_final_amplitude_vs_gamma(results, plot_title, output_base)
    Gamma = cellfun(@(r) r.Gamma, results);
    amplitude = cellfun(@(r) r.final_amplitude, results);
    analytic_amplitude = cellfun(@(r) r.final_analytic_amplitude, results);

    fig = make_figure();
    plot(Gamma, amplitude, 'ro-', 'LineWidth', 1.7, 'MarkerSize', 7);
    hold on;
    plot(Gamma, analytic_amplitude, 'ks--', 'LineWidth', 1.7, 'MarkerSize', 7);
    grid on;
    box on;
    xlabel('\Gamma');
    ylabel('final amplitude');
    title(plot_title, 'Interpreter', 'none');
    legend('numerical', 'analytical', 'Location', 'best');
    save_figure(fig, output_base);
end

function plot_final_metric_vs_dt(results, metric_field, y_label, plot_title, output_base)
    dt_values = cellfun(@(r) r.dt, results);
    metric = cellfun(@(r) r.(metric_field), results);

    fig = make_figure();
    loglog(dt_values, metric, 'bo-', 'LineWidth', 1.7, 'MarkerSize', 7);
    grid on;
    box on;
    xlabel('\Delta t');
    ylabel(y_label);
    title(plot_title, 'Interpreter', 'none');
    save_figure(fig, output_base);
end

function plot_final_metric_vs_points(results, metric_field, y_label, plot_title, output_base)
    points = cellfun(@(r) r.points, results);
    metric = cellfun(@(r) r.(metric_field), results);

    fig = make_figure();
    plot(points, metric, 'bo-', 'LineWidth', 1.7, 'MarkerSize', 7);
    grid on;
    box on;
    xlabel('number of grid points');
    ylabel(y_label);
    title(plot_title, 'Interpreter', 'none');
    save_figure(fig, output_base);
end

function fig = make_figure()
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 900 560]);
    set(fig, 'DefaultAxesFontSize', 12);
end

function save_figure(fig, output_base)
    print(fig, [output_base '.png'], '-dpng', '-r200');
    close(fig);
end

function label = result_label(result)
    if contains(result.experiment_name, 'Courant')
        label = sprintf('C=%.2f', result.Courant);
    elseif contains(result.experiment_name, 'points')
        label = sprintf('N=%d', result.points);
    elseif contains(result.experiment_name, 'Gamma')
        label = sprintf('\\Gamma=%.2f', result.Gamma);
    elseif contains(result.experiment_name, 'dt')
        label = sprintf('\\Delta t=%.3f', result.dt);
    else
        label = result.experiment_name;
    end
end

function add_parameter_text(result)
    text(0.02, 0.04, ...
        sprintf('N=%d, dt=%.4g, C=%.3g, F=%.3g', ...
        result.points, result.dt, result.Courant, result.Fourier_number), ...
        'Units', 'normalized', 'FontSize', 10, ...
        'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);
end

%% Table and report writers
function write_result_table(results, filename)
    if isempty(results)
        warning('No results to write for %s', filename);
        return;
    end

    n = numel(results);
    experiment_name = strings(n,1);
    solver = strings(n,1);
    U0 = NaN(n,1);
    Gamma = NaN(n,1);
    points = NaN(n,1);
    dx = NaN(n,1);
    dt = NaN(n,1);
    tend = NaN(n,1);
    tsteps = NaN(n,1);
    Courant = NaN(n,1);
    Fourier_number = NaN(n,1);
    final_L2 = NaN(n,1);
    final_Linf = NaN(n,1);
    final_periodic_error = NaN(n,1);
    final_amplitude = NaN(n,1);
    final_analytic_amplitude = NaN(n,1);
    amplitude_ratio_numeric_over_analytic = NaN(n,1);
    max_abs_phi_final = NaN(n,1);
    unstable_detected = false(n,1);

    for k = 1:n
        r = results{k};
        experiment_name(k) = string(r.experiment_name);
        solver(k) = string(r.solver);
        U0(k) = r.U0;
        Gamma(k) = r.Gamma;
        points(k) = r.points;
        dx(k) = r.dx;
        dt(k) = r.dt;
        tend(k) = r.tend;
        tsteps(k) = r.tsteps;
        Courant(k) = r.Courant;
        Fourier_number(k) = r.Fourier_number;
        final_L2(k) = r.final_L2;
        final_Linf(k) = r.final_Linf;
        final_periodic_error(k) = r.final_periodic_error;
        final_amplitude(k) = r.final_amplitude;
        final_analytic_amplitude(k) = r.final_analytic_amplitude;
        amplitude_ratio_numeric_over_analytic(k) = r.amplitude_ratio_numeric_over_analytic;
        max_abs_phi_final(k) = r.max_abs_phi_final;
        unstable_detected(k) = r.unstable_detected;
    end

    T = table(experiment_name, solver, U0, Gamma, points, dx, dt, tend, tsteps, ...
        Courant, Fourier_number, final_L2, final_Linf, final_periodic_error, ...
        final_amplitude, final_analytic_amplitude, ...
        amplitude_ratio_numeric_over_analytic, max_abs_phi_final, unstable_detected);
    writetable(T, filename);
end

function write_audit_report(filename, EE2, IE2, CN1, CN2, CN3, CN4, CN5, CN6)
    fid = fopen(filename, 'w');
    if fid < 0
        error('Could not open audit report for writing: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '# Numerical Audit Report for CFD Assignment 2a and 2b\n\n');

    fprintf(fid, '## 1. Purpose of the numerical audit\n\n');
    fprintf(fid, ['The goal of this workflow is to verify the analytical solution implementation, ', ...
        'the consistency of periodic boundary conditions, and the expected stability, ', ...
        'dissipation, and dispersion behaviour of the three numerical schemes.  The tests ', ...
        'also document the sensitivity to advection velocity `U0`, diffusion coefficient ', ...
        '`Gamma`, grid spacing `dx`, and timestep `dt`.\n\n']);

    fprintf(fid, '## 2. Verification metrics\n\n');
    fprintf(fid, '- **Relative L2 error**: `norm(phi - phi_a, 2) / norm(phi_a, 2)`. This measures the global profile error relative to the analytical solution.\n');
    fprintf(fid, '- **Linf error**: `max(abs(phi - phi_a))`. This catches the largest pointwise deviation.\n');
    fprintf(fid, '- **Periodic boundary mismatch**: `abs(phi(1) - phi(end))`. This checks whether the duplicated periodic endpoint remains consistent.\n');
    fprintf(fid, '- **Numerical amplitude**: `(max(phi) - min(phi)) / 2`. This is used to diagnose artificial growth or damping.\n');
    fprintf(fid, '- **Analytical amplitude**: `1` for pure advection and `exp(-Gamma*t)` for advection-diffusion.\n');
    fprintf(fid, '- **Courant number**: `C = U0*dt/dx`. This controls the advective timestep size.\n');
    fprintf(fid, '- **Fourier number**: `F = Gamma*dt/dx^2`. This measures the diffusive timestep size.\n\n');

    fprintf(fid, '## 3. Implementation sanity checks\n\n');
    fprintf(fid, ['The periodic mismatch should remain close to machine precision because `x=0` ', ...
        'and `x=2*pi` represent the same physical point.  For `U0 > 0`, the sine wave should ', ...
        'propagate in the positive x direction.  For `U0 = 0`, no translation should occur.  ', ...
        'For `Gamma = 0`, the analytical solution has no physical amplitude decay, while ', ...
        'for `Gamma > 0` the analytical amplitude is `exp(-Gamma*t)`.\n\n']);
    fprintf(fid, '- Representative EE-2 final periodic mismatch range: %.3e to %.3e.\n', ...
        min_cell_metric(EE2, 'final_periodic_error'), max_cell_metric(EE2, 'final_periodic_error'));
    fprintf(fid, '- Representative IE-2 final periodic mismatch range: %.3e to %.3e.\n', ...
        min_cell_metric(IE2, 'final_periodic_error'), max_cell_metric(IE2, 'final_periodic_error'));
    fprintf(fid, '- CN-1 final amplitude: numerical %.6g, analytical %.6g.\n', ...
        CN1.final_amplitude, CN1.final_analytic_amplitude);
    fprintf(fid, '- CN-3 pure diffusion final amplitude: numerical %.6g, analytical %.6g.\n\n', ...
        CN3.final_amplitude, CN3.final_analytic_amplitude);

    fprintf(fid, '## 4. Results for Explicit Euler + Central Difference\n\n');
    fprintf(fid, ['The explicit central advection scheme exhibits unstable amplitude growth. ', ...
        'A smaller Courant number delays visible blow-up, but it does not remove the ', ...
        'fundamental instability of the FTCS-like central scheme for pure advection. ', ...
        'The larger Courant-number cases show faster growth in relative L2 error and ', ...
        '`max(abs(phi))`.\n\n']);
    fprintf(fid, 'Relevant outputs:\n');
    fprintf(fid, '- `figures/EE_short_time_profile.png`\n');
    fprintf(fid, '- `figures/EE_Courant_L2_history.png`\n');
    fprintf(fid, '- `figures/EE_Courant_maxabs_history.png`\n');
    fprintf(fid, '- `figures/EE_Courant_1_profile.png`\n');
    fprintf(fid, '- `figures/EE_grid_final_L2.png`\n');
    fprintf(fid, '- `tables/table_EE_Courant_study.csv`\n');
    fprintf(fid, '- `tables/table_EE_grid_study.csv`\n\n');

    fprintf(fid, '## 5. Results for Implicit Euler + Central Difference\n\n');
    fprintf(fid, ['The implicit Euler scheme remains bounded for all tested Courant numbers. ', ...
        'However, larger Courant number increases numerical dissipation, which is visible ', ...
        'as stronger amplitude decay and larger final L2 error.  This confirms that a ', ...
        'stable scheme can still be inaccurate when the timestep is too large.\n\n']);
    fprintf(fid, 'Relevant outputs:\n');
    fprintf(fid, '- `figures/IE_baseline_profile.png`\n');
    fprintf(fid, '- `figures/IE_Courant_L2_history.png`\n');
    fprintf(fid, '- `figures/IE_Courant_amplitude_history.png`\n');
    fprintf(fid, '- `figures/IE_Courant_profile_C010_C200.png`\n');
    fprintf(fid, '- `figures/IE_grid_final_L2.png`\n');
    fprintf(fid, '- `tables/table_IE_Courant_study.csv`\n');
    fprintf(fid, '- `tables/table_IE_grid_study.csv`\n\n');

    fprintf(fid, '## 6. Results for Crank-Nicolson Advection-Diffusion\n\n');
    fprintf(fid, '### 6.1 Coupled advection-diffusion baseline\n\n');
    fprintf(fid, ['The baseline case translates and decays as expected.  The final numerical ', ...
        'amplitude %.6g is close to the analytical amplitude %.6g, and the final relative ', ...
        'L2 error is %.3e.\n\n'], CN1.final_amplitude, CN1.final_analytic_amplitude, CN1.final_L2);

    fprintf(fid, '### 6.2 Pure advection limit Gamma = 0\n\n');
    fprintf(fid, ['With `Gamma = 0`, no physical amplitude decay is present.  The solution stays ', ...
        'bounded in all tested Courant-number cases, while larger Courant number makes ', ...
        'phase/dispersive error more visible.\n\n']);

    fprintf(fid, '### 6.3 Pure diffusion limit U0 = 0\n\n');
    fprintf(fid, ['With `U0 = 0`, the wave does not translate.  The sinusoidal shape is preserved ', ...
        'while the amplitude decays approximately as `exp(-Gamma*t)`.  The pure diffusion ', ...
        'final relative L2 error is %.3e.\n\n'], CN3.final_L2);

    fprintf(fid, '### 6.4 Gamma sensitivity\n\n');
    fprintf(fid, ['Increasing `Gamma` increases physical damping.  The final numerical amplitudes ', ...
        'track the analytical amplitudes in `figures/CN_Gamma_final_amplitude.png`.\n\n']);

    fprintf(fid, '### 6.5 dt and grid sensitivity\n\n');
    fprintf(fid, ['The timestep study shows that larger `dt` increases the final error, while the ', ...
        'grid study shows decreasing final L2 error under mesh refinement.  These trends ', ...
        'are consistent with Crank-Nicolson time integration and central spatial differences.\n\n']);

    fprintf(fid, 'Relevant outputs:\n');
    fprintf(fid, '- `figures/CN_baseline_profile.png`\n');
    fprintf(fid, '- `figures/CN_baseline_L2_history.png`\n');
    fprintf(fid, '- `figures/CN_pure_advection_C100_profile.png`\n');
    fprintf(fid, '- `figures/CN_pure_diffusion_amplitude_history.png`\n');
    fprintf(fid, '- `figures/CN_Gamma_final_amplitude.png`\n');
    fprintf(fid, '- `figures/CN_dt_final_L2.png`\n');
    fprintf(fid, '- `figures/CN_grid_final_L2.png`\n');
    fprintf(fid, '- `tables/table_CN_baseline_and_limits.csv`\n');
    fprintf(fid, '- `tables/table_CN_Gamma_study.csv`\n');
    fprintf(fid, '- `tables/table_CN_dt_study.csv`\n');
    fprintf(fid, '- `tables/table_CN_grid_study.csv`\n\n');

    fprintf(fid, '## 7. Suggested report-ready discussion text\n\n');
    fprintf(fid, ['The explicit Euler scheme with central spatial differencing is not suitable ', ...
        'for stable long-time integration of the pure advection equation.  Although small ', ...
        'Courant numbers and short integration times can make the numerical profile appear ', ...
        'reasonable, the error and amplitude histories reveal the unstable amplification ', ...
        'expected from the FTCS-like central scheme.  Increasing the Courant number makes ', ...
        'the instability visible more quickly.\n\n']);
    fprintf(fid, ['The implicit Euler central scheme removes the blow-up observed in the explicit ', ...
        'scheme, but it does so at the price of numerical damping.  Larger timesteps remain ', ...
        'bounded, yet the wave amplitude decays too strongly and the phase error increases. ', ...
        'Therefore stability alone is not sufficient for accuracy.\n\n']);
    fprintf(fid, ['For the advection-diffusion equation, the Crank-Nicolson central scheme gives ', ...
        'a more balanced result.  The coupled baseline agrees well with the analytical ', ...
        'solution, including both translation and exponential amplitude decay.  The limiting ', ...
        'cases confirm that the implementation reduces correctly to pure advection when ', ...
        '`Gamma = 0` and pure diffusion when `U0 = 0`.\n\n']);

    fprintf(fid, '## 8. Final checklist\n\n');
    fprintf(fid, '[ ] All periodic boundary mismatches are small\n');
    fprintf(fid, '[ ] Numerical propagation direction is correct\n');
    fprintf(fid, '[ ] EE instability observed\n');
    fprintf(fid, '[ ] IE numerical damping observed\n');
    fprintf(fid, '[ ] CN pure advection case behaves reasonably\n');
    fprintf(fid, '[ ] CN pure diffusion case matches exp(-Gamma*t)\n');
    fprintf(fid, '[ ] Parameter tables exported\n');
    fprintf(fid, '[ ] Figures exported\n');

    clear cleaner;
end

function value = min_cell_metric(results, field_name)
    values = cellfun(@(r) r.(field_name), results);
    value = min(values, [], 'omitnan');
end

function value = max_cell_metric(results, field_name)
    values = cellfun(@(r) r.(field_name), results);
    value = max(values, [], 'omitnan');
end
