%% Replicates Table 1A and 1B of Smets & Wouters (2007)
clear; clc;

%% Loading Dynare results
baseName = 'sw_model';
load(fullfile(baseName, 'Output', [baseName '_results.mat']));

param_names = bayestopt_.name;
nParams = length(param_names);

%% Extracting posterior statistics
post_mean = zeros(nParams, 1);
post_std  = zeros(nParams, 1);
hpd_low   = zeros(nParams, 1);
hpd_high  = zeros(nParams, 1);

for i = 1:nParams
    pname = param_names{i};
    if isfield(oo_.posterior_mean.parameters, pname)
        post_mean(i) = oo_.posterior_mean.parameters.(pname);
        post_std(i)  = oo_.posterior_std.parameters.(pname);
        hpd_low(i)   = oo_.posterior_hpdinf.parameters.(pname);
        hpd_high(i)  = oo_.posterior_hpdsup.parameters.(pname);
    elseif isfield(oo_.posterior_mean.shocks_std, pname)
        post_mean(i) = oo_.posterior_mean.shocks_std.(pname);
        post_std(i)  = oo_.posterior_std.shocks_std.(pname);
        hpd_low(i)   = oo_.posterior_hpdinf.shocks_std.(pname);
        hpd_high(i)  = oo_.posterior_hpdsup.shocks_std.(pname);
    end
end

%% Table 1A - Structural parameters
table1A_params = {'phi','sigma_c','lambda','xi_w','sigma_l','xi_p','iota_w','iota_p',...
    'psi','phi_p','r_pi','rho','r_y','r_dy','pi_bar','beta_const','l_bar','alpha',...
    'rho_a','rho_b','rho_g','rho_i','rho_r','rho_p','rho_w','Mu_p','Mu_w','rho_ga'};

fprintf('\n===============================================================\n');
fprintf('  TABLE 1A: Structural & Shock Process Parameters (No Growth)\n');
fprintf('===============================================================\n');
fprintf('%-12s %10s %10s %10s %10s\n', 'Parameter', 'Mean', 'Std', 'HPD5%', 'HPD95%');
fprintf('---------------------------------------------------------------\n');
for j = 1:length(table1A_params)
    idx = find(strcmp(param_names, table1A_params{j}));
    if ~isempty(idx)
        fprintf('%-12s %10.4f %10.4f %10.4f %10.4f\n', ...
            table1A_params{j}, post_mean(idx), post_std(idx), hpd_low(idx), hpd_high(idx));
    end
end

%% Table 1B - Shock standard deviations
table1B_params = {'eta_a','eta_b','eta_g','eta_i','eta_r','eta_p','eta_w'};

fprintf('\n===============================================================\n');
fprintf('  TABLE 1B: Standard Deviations of Structural Shocks\n');
fprintf('===============================================================\n');
fprintf('%-12s %10s %10s %10s %10s\n', 'Shock', 'Mean', 'Std', 'HPD5%', 'HPD95%');
fprintf('---------------------------------------------------------------\n');
for j = 1:length(table1B_params)
    idx = find(strcmp(param_names, table1B_params{j}));
    if ~isempty(idx)
        fprintf('%-12s %10.4f %10.4f %10.4f %10.4f\n', ...
            table1B_params{j}, post_mean(idx), post_std(idx), hpd_low(idx), hpd_high(idx));
    end
end
fprintf('---------------------------------------------------------------\n');
