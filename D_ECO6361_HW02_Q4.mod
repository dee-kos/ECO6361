// RBC Model (without nominal rigidities)

// Declaring endogenous variables
var y c inv n r k a w;
varexo epsilon_a;

// Declaring parameters
parameters alphaa betaa deltaa phii rho_a lambdaa omegaa psii;

// Assiging given parameter values
alphaa = 0.35;   
betaa = 0.99;    
deltaa = 0.025;  
phii = 1.25;     
rho_a = 0.92;   
lambdaa = 0.75;  
omegaa = 1.25;   
psii = 0.6;      

// Modeling equations
model;
    // 1. Euler equation (Household consumption smoothing)
    c = c(+1) - (r);

    // 2. Capital accumulation equation
    k(+1) = (1 - deltaa) * k + inv;

    // 3. Firm production function
    y = alphaa * k + (1 - alphaa) * n + a;

    // 4. Wage equation (log-linearized marginal product of labor)
    w = (1 - alphaa) * y - n;

    // 5. Interest rate determined by marginal product of capital
    r = alphaa * y / k;

    // 6. AR(1) process for technology shock
    a = rho_a * a(-1) + epsilon_a;

    // 7. Goods market equilibrium
    y = c + inv;

    // 8. Labor supply equation (relating labor to wages)
    n = phii * w;
end;

// Initial guesses for steady-state values to reach SS
initval;
    y = 2;        // Adjusted steady-state output
    c = 1.6;      // Adjusted steady-state consumption
    inv = 0.4;    // Adjusted steady-state investment
    n = 0.5;      // Adjusted steady-state labor
    k = 8;        // Adjusted steady-state capital
    r = 0.05;     // Adjusted steady-state interest rate
    w = 1.0;      // Adjusted steady-state wage
    a = 0;        // Technology shock initially 0
end;

shocks;
    var epsilon_a = 0.01^2;  // Variance of technology shock
end;

// Steady state and stochastic simulations
steady;
stoch_simul(order=1, irf=50);
