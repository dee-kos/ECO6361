%% Generates IRF figures for 3 shocks x 7 variables x 6 calibrations
clc; clear; close all;

baseName = 'sw_model';
modTxt = fileread([baseName '.mod']);
modelEnd = strfind(modTxt, 'end;');
modelStart = strfind(modTxt, 'model(');
modelEnd = modelEnd(find(modelEnd > modelStart, 1)) + 3;
modelTxt = modTxt(1:modelEnd);

%% Loading posterior draws
mcmcPath = fullfile(baseName, 'metropolis');
S1 = load(fullfile(mcmcPath, [baseName '_mh1_blck1.mat']), 'x2', 'logpo2');
S2 = load(fullfile(mcmcPath, [baseName '_mh1_blck2.mat']), 'x2', 'logpo2');
draws = [S1.x2; S2.x2];
logpo = [S1.logpo2; S2.logpo2];
nParams = size(draws, 2);

%% Parameter ordering
ordered = {'phi','sigma_c','lambda','xi_w','sigma_l','xi_p','iota_w','iota_p',...
    'psi','phi_p','r_pi','rho','r_y','r_dy','pi_bar','beta_const','l_bar','alpha',...
    'eta_a','eta_b','eta_g','eta_i','eta_r','eta_p','eta_w',...
    'rho_a','rho_b','rho_g','rho_i','rho_r','rho_p','rho_w','Mu_p','Mu_w','rho_ga'};
shock_std = {'eta_a','eta_b','eta_g','eta_i','eta_r','eta_p','eta_w'};

%% Prior values
theta_prior = [4,1.5,0.7,0.5,2,0.5,0.5,0.5,0.5,1.25,1.5,0.75,0.125,0.125,0.625,0.25,0,0.3,...
    0.1,0.1,0.1,0.1,0.1,0.1,0.1, 0.5,0.5,0.5,0.5,0.5,0.5,0.5, 0.5,0.5,0.5];

%% Building 6 calibrations
[~, iMode] = max(logpo);
CAL = struct('label', {'Post Mean','Post Median','Post Mode','Prior Mean','Prior Median','Prior Mode'}, ...
    'theta', {mean(draws,1), median(draws,1), draws(iMode,:), theta_prior, theta_prior, theta_prior});

%% IRF
vars = {'dy','dc','dinve','labobs','pinfobs','dw','robs'};
varLabels = {'Output Growth','Consumption Growth','Investment Growth','Hours Worked','Inflation','Wage Growth','Interest Rate'};
shocks = {'eta_a','eta_i','eta_r'};
shockLabels = {'Technology Shock','Investment Shock','Monetary Policy Shock'};
T = 40;

% Initializing storage
IRF = struct();
for s = 1:3, for v = 1:7, IRF.(shocks{s}).(vars{v}) = zeros(T,6); end, end

%% Computing IRFs for each calibration
for c = 1:6
    fprintf('Running calibration %d/6: %s\n', c, CAL(c).label);
    theta = CAL(c).theta;
    
    %temp mod file
    fid = fopen('tmp.mod', 'w');
    fprintf(fid, '%s\n\n', modelTxt);
    for p = 1:nParams
        if ~ismember(ordered{p}, shock_std)
            fprintf(fid, '%s = %.12f;\n', ordered{p}, theta(p));
        end
    end
    fprintf(fid, '\nshocks;\n');
    for s = 1:7
        idx = find(strcmp(ordered, shock_std{s}));
        fprintf(fid, '    var %s; stderr %.12f;\n', shock_std{s}, theta(idx));
    end
    fprintf(fid, 'end;\n\nstoch_simul(order=1, irf=%d, nograph, noprint);\n', T);
    fclose(fid);
    
    % Running dynare
    dynare('tmp', 'noclearall', 'nolog');
    
    % Extracting IRFs
    for s = 1:3
        for v = 1:7
            f = [vars{v} '_' shocks{s}];
            if isfield(oo_.irfs, f), IRF.(shocks{s}).(vars{v})(:,c) = oo_.irfs.(f)(1:T); end
        end
    end
end

% Cleaning
delete('tmp.mod'); rmdir('tmp', 's');

%% Plotting figures
colors = {[0 .45 .74],[.85 .33 .1],[.47 .67 .19],[.5 .5 .5],[.93 .69 .13],[.49 .18 .56]};
styles = {'-','-','-','--','--','--'};
widths = [2 2 2 1.5 1.5 1.5];

for s = 1:3
    figure('Position', [50 50 1400 900], 'Color', 'w');
    for v = 1:7
        subplot(3,3,v); hold on;
        plot(1:T, zeros(T,1), 'k-', 'LineWidth', 0.5);
        for c = 1:6
            plot(1:T, IRF.(shocks{s}).(vars{v})(:,c), 'Color', colors{c}, 'LineStyle', styles{c}, 'LineWidth', widths(c));
        end
        xlabel('Quarters'); ylabel('% deviation'); title(varLabels{v}); xlim([1 T]); grid on;
    end
    subplot(3,3,8); axis off; hold on;
    for c = 1:6, plot(NaN, NaN, 'Color', colors{c}, 'LineStyle', styles{c}, 'LineWidth', widths(c)); end
    legend({CAL.label}, 'Location', 'best', 'FontSize', 10);
    sgtitle(['IRFs to ' shockLabels{s}], 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, ['IRF_' shocks{s} '.png']);
end
fprintf('Done! Figures saved.\n');
