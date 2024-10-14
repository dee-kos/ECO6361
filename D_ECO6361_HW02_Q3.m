% Creating a 2x3 tiled layout for multiple graphs
figure;
tiledlayout(2,3);

% Plotting Output (y) in response to technology shock
nexttile;
plot(oo_.irfs.y_epsilon_a);
title('Y to Tec Shock');

% Plotting Consumption (c) in response to technology shock
nexttile;
plot(oo_.irfs.c_epsilon_a);
title('C to Tec Shock');

% Plotting Investment (inv) in response to technology shock
nexttile;
plot(oo_.irfs.inv_epsilon_a);
title('Inv to Tec Shock');

% Plotting Interest rate (r) in response to monetary shock
nexttile;
plot(oo_.irfs.r_epsilon_M);
title('IR to Monetary Shock');

% Plotting Inflation (Pi) in response to monetary shock
nexttile;
plot(oo_.irfs.Pii_epsilon_M);
title('Inf to Monetary Shock');

% Plotting Wages (w) in response to technology shock
nexttile;
plot(oo_.irfs.w_epsilon_a);
title('Wages to Tec Shock');
