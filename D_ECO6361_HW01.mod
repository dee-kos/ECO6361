// Question 1 - full system of equations
//I'm in the second half of the students, so non-linear form is used

//Endogenous variables
var c n k w r a y inv pi;

//exogenous shock
varexo e;

// Parameters
parameters beta gamma phi delta alpha rho sigma Psi;

// Calibration
beta = 0.99;       // Discount factor
gamma = 1;         
phi = 2;           // Labor supply elasticity
delta = 0.025;     // Depreciation rate
alpha = 0.35;      // Capital share in production
rho = 0.92;        // Persistence of technology shock
sigma = 0.04;      // Std. dev. of technology shock

//Variabe definitions at steady state (ss)

n_ss=38/168;

a_ss=1;

k_ss=((n_ss^(1-alpha))*a_ss/((1/beta-(1-delta))/alpha))^(1/(alpha-1));

y_ss=a_ss*k_ss^alpha*n_ss^(1-alpha);

c_ss=y_ss-delta*k_ss;


//calculation of psi given steady state

Psi=((1-alpha)*y_ss)/(c_ss*n_ss^phi);


// Model equations

model;

// 1. Household's Euler equation
(1/c) = beta*(1/c(+1))*(1+r(+1)-delta);

// 2. Household's labor supply equation
Psi * n^phi= (1-alpha)*y/c;

// 3. Production function
y = a * k^alpha * n^(1 - alpha);

// 4. Capital law of motion
k = (1 - delta) * k(-1) + inv;

// 5. Resource constraint (social planner's constraint)
y=c+inv;

// 6. Firm's first-order condition for labor
w = (1 - alpha) * y/n;

// 7. Firm's first-order condition for capital
r = alpha * y/k;

// 8. Profit function (assuming zero profits in steady state)
pi = y - w * n - r * k;

// 9. Technology shock process
log(a) = rho * log(a(-1)) + e;

end;


// Initial values for steady state

initval;
    c = c_ss;        // Close to expected consumption
    n = n_ss;        // Adjusted guess for labor supply
    k = k_ss;         // Initial guess for capital
    inv=delta*k_ss; // Based on k and delta
    w = (1-alpha)*y_ss/n_ss;       // Initial guess for wage
    r = alpha*y_ss/k_ss;       // Initial guess for interest rate on bonds
    a = 1;          // Technology normalized
    y = y_ss;          // Initial guess for output
    pi = 0;         // Zero profits assumed
    e = 0;          // No shock in steady state
end;

//shocks

shocks;
var e; stderr sigma;
end;

// Question 2 - Computing the steady state

steady;


// Running stochastic simulation and generating impulse response functions

stoch_simul(order=1, irf=50, periods=200);