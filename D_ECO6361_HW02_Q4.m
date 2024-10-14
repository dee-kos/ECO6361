% Running the NK simulation
dynare D_ECO6361_HW02_Q2.mod;

% Storing NK model IRFs in a separate structure
irfs_nk = oo_.irfs;

% Running the RBC simulation
dynare D_ECO6361_HW02_Q4.mod;

% Storing RBC model IRFs in a separate structure
irfs_rbc = oo_.irfs;

% Creating a 2x2 tiled layout for comparing key variables
figure;
tiledlayout(2,2);

% Comparing Output (y) responses to technology shock
nexttile;
plot(irfs_nk.y_epsilon_a); hold on;
plot(irfs_rbc.y_epsilon_a);
title('Output (y) response to Technology Shock');
legend('NK', 'RBC');

% Comparing Consumption (c) responses to technology shock
nexttile;
plot(irfs_nk.c_epsilon_a); hold on;
plot(irfs_rbc.c_epsilon_a);
title('Consumption (c) response to Technology Shock');
legend('NK', 'RBC');

% Comparing Investment (inv) responses to technology shock
nexttile;
plot(irfs_nk.inv_epsilon_a); hold on;
plot(irfs_rbc.inv_epsilon_a);
title('Investment (inv) response to Technology Shock');
legend('NK', 'RBC');

% Comparing Labor (n) responses to technology shock
nexttile;
plot(irfs_nk.n_epsilon_a); hold on;
plot(irfs_rbc.n_epsilon_a);
title('Labor (n) response to Technology Shock');
legend('NK', 'RBC');
