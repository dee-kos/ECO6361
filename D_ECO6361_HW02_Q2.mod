// Declaring endogenous variables
var y c inv n r Pii k a M w;
varexo epsilon_a epsilon_M;

// Declaring parameters
parameters alphaa betaa deltaa phii rho_a rho_M rho_i zeta_pi zeta_y lambdaa omegaa psii S_double_prime xi_p xi_w eta_w kappaa;

// Assigning parameter values
alphaa = 0.35;
betaa = 0.99;
deltaa = 0.025;
phii = 1.25;
rho_a = 0.92;
rho_M = 0.4;
rho_i = 0.78;
zeta_pi = 1.5;
zeta_y = 0.1;
lambdaa = 0.75;
omegaa = 1.25;
psii = 0.6;
S_double_prime = 5;
xi_p = 0.6;
xi_w = 0.9;
eta_w = 2.1;
kappaa = 0.1;

// Model equations (log-linearized)
model;
    // 1. Euler equation
    c = c(+1) - (r - Pii(+1));

    // 2. Capital accumulation 
    k(+1) = (1 - deltaa) * k + inv;

    // 3. Firm production function
    y = alphaa * k + (1 - alphaa) * n + a;

    // 4. Wage equation (log-linearized marginal product of labor)
    w = (1 - alphaa) * y - n;

    // 5. New Keynesian Phillips Curve (Price setting equation)
    Pii = betaa * Pii(+1) + kappaa * y;

    // 6. Taylor Rule (Monetary Policy Rule)
    r = rho_i * r(-1) + zeta_pi * Pii + zeta_y * y + epsilon_M;

    // 7. AR(1) process for technology shock
    a = rho_a * a(-1) + epsilon_a;

    // 8. AR(1) process for monetary shock
    M = rho_M * M(-1) + epsilon_M;

    // 9. Goods market equilibrium
    y = c + inv;

    // 10. Labor supply equation (relating labor to wages)
    n = phii * w;
end;

// Shocks variances
shocks;
    var epsilon_a = 0.01^2;  // Variance of technology shock
    var epsilon_M = 0.01^2;  // Variance of monetary policy shock
end;

// Steady state and stochastic simulations
steady;
stoch_simul(order=1, irf=50);
