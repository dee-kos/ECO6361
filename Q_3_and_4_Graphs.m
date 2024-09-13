% Question 3 - Using tiledlayout function to see IR to shock e

% List of variable names and their corresponding titles
variables = {'y', 'c', 'inv', 'n', 'r', 'w', 'a'};
titles = {'Output', 'Consumption', 'Investment', 'Labor Supply', ...
          'Interest Rate', 'Wage', 'Technology'};

shock_name = '_e';

% Creating a tiled layout for plotting IRFs
figure;
tiledlayout(3,3);

for i = 1:length(variables)
    irf_field = [variables{i}, shock_name];

        irf_data = oo_.irfs.(irf_field);

        % Plotting the IRF
        nexttile;
        plot(irf_data, 'LineWidth', 1.5);
        title(titles{i});
        xlabel('Periods');
        ylabel('IRF');
end

sgtitle('Impulse Response Functions for Endogenous Variables');


% Question 4 - Using tiledlayout function to graph autocorrelations

% Creating a tiled layout for plotting autocorrelations
figure;
tiledlayout(3,3);

for i = 1:length(variables)
    % Finding the index of the variable in the Dynare endogenous names list
    var_index = strmatch(variables{i}, M_.endo_names, 'exact');
    
    % Extracting the simulated data for the variable
    data = oo_.endo_simul(var_index, :)';
    
    % Calculating autocorrelation for the first 6 lags
    acf_data = autocorr(data,6);

    % Plotting the autocorrelation function
    nexttile;
    stem(0:6, acf_data, 'filled'); % Plot with lags 0 to 6
    title(titles{i});
    xlabel('Lags');
    ylabel('Autocorrelation');
end

sgtitle('Autocorrelation Functions for Endogenous Variables');
