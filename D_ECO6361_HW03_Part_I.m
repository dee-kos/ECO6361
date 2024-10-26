% Loading data
data = readtable('Quarterly.csv');
time = datetime(data.DATE, 'InputFormat', 'yyyy-MM-dd');

% Defining variables based on data, other than to Smets and Wouters (2007), here I referred to,
% Smets, Frank, and Wouters, Rafael. Replication data for: Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach.
% Nashville, TN: American Economic Association [publisher], 2007.
% Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2019-12-07.
% https://doi.org/10.3886/E116269V1
 
consumption = log((data.PCEC ./ data.GDPDEF) ./ data.CLF16OV) * 100;
% PCEC: Personal Consumption Expenditures - Billions of Dollars, Seasonally Adjusted Annual Rate
% GDPDEF:Gross Domestic Product: Implicit Price Deflator, Index 2017=100, Seasonally Adjusted (GDPDEF)
% CLF160V: Civilian Labor Force Level, Thousands of Persons, Seasonally Adjusted (CLF16OV)

investment = log((data.FPI ./ data.GDPDEF) ./ data.CLF16OV) * 100;
% FPI: Civilian Labor Force Level, Thousands of Persons, Seasonally Adjusted (CLF16OV)

output = log(data.GDPC1 ./ data.CLF16OV) * 100;
% GDPC1: Real Gross Domestic Product, Billions of Chained 2017 Dollars, Seasonally Adjusted Annual Rate (GDPC1)

hours_worked = log((data.PRS85006023 .* data.CE16OV / 100) ./ data.CLF16OV) * 100;
% PRS85006023:Nonfarm Business, All Persons, Average Weekly Hours Duration : index, 2017 = 100, Seasonally Adjusted
% CE16OV:Civilian Employment: Sixteen Years & Over, Thousands, 2017=100, Seasonally Adjusted

inflation = [NaN; log(data.GDPDEF(2:end) ./ data.GDPDEF(1:end-1)) * 100];

real_wage = log(data.COMPRNFB ./ data.GDPDEF) * 100;
% COMPRNFB: Nonfarm Business Sector: Real Hourly Compensation for All Workers, Index 2017=100, Seasonally Adjusted (COMPRNFB)

interest_rate = data.FEDFUNDS;
% FEDFUNDS:Federal Funds Effective Rate, Percent, Not Seasonally Adjusted (FEDFUNDS)

% Setting parameters for HP filter quaterly data
lambdaa = 1600;

% Detrending and demeaning each series

[~, consumption_cycle] = hpfilter(consumption, lambdaa);
consumption_demeaned = consumption_cycle - mean(consumption_cycle, 'omitnan');

[~, investment_cycle] = hpfilter(investment, lambdaa);
investment_demeaned = investment_cycle - mean(investment_cycle, 'omitnan');

[~, output_cycle] = hpfilter(output, lambdaa);
output_demeaned = output_cycle - mean(output_cycle, 'omitnan');

[~, hours_cycle] = hpfilter(hours_worked, lambdaa);
hours_demeaned = hours_cycle - mean(hours_cycle, 'omitnan');

[~, inflation_cycle] = hpfilter(inflation, lambdaa);
inflation_demeaned = inflation_cycle - mean(inflation_cycle, 'omitnan');

[~, real_wage_cycle] = hpfilter(real_wage, lambdaa);
real_wage_demeaned = real_wage_cycle - mean(real_wage_cycle, 'omitnan');

[~, interest_rate_cycle] = hpfilter(interest_rate, lambdaa);
interest_rate_demeaned = interest_rate_cycle - mean(interest_rate_cycle, 'omitnan');

% Defining NBER defined recessions

NBER_recessions = [
    datetime(1960,1,1), datetime(1961,1,1);
    datetime(1969,12,1), datetime(1970,11,1);
    datetime(1973,11,1), datetime(1975,3,1);
    datetime(1980,1,1), datetime(1980,7,1);
    datetime(1981,7,1), datetime(1982,11,1);
    datetime(1990,7,1), datetime(1991,3,1);
    datetime(2001,3,1), datetime(2001,11,1);
    datetime(2007,12,1), datetime(2009,6,1);
    datetime(2020,2,1), datetime(2020,4,1)
];

% Plotting series in a single tile plot with recessions

figure;
variables = {consumption_demeaned, investment_demeaned, output_demeaned, ...
             hours_demeaned, inflation_demeaned, real_wage_demeaned, interest_rate_demeaned};

titles = {'Consumption', 'Investment', 'Output', 'Hours Worked', 'Inflation', 'Real Wage', 'Interest Rate'};

for i = 1:length(variables)
    if length(time) ~= length(variables{i})
        time_adj = time(1:length(variables{i}));
    else
        time_adj = time;
    end

    subplot(3, 3, i);
    
    plot(time_adj, variables{i}, 'LineWidth', 1.2);
    hold on;
    

    yl = ylim; % Get y-axis limits

    % Recession Shading
    for j = 1:size(NBER_recessions, 1)
        x = [NBER_recessions(j,1), NBER_recessions(j,2), NBER_recessions(j,2), NBER_recessions(j,1)];
        y = [yl(1), yl(1), yl(2), yl(2)];
        fill(x, y, [0.9 0.9 0.9], 'EdgeColor', [0.9 0.9 0.9]);
    end

    ylim(yl);

    hold on;

    plot(time_adj, variables{i}, 'LineWidth', 1.2);

    title(titles{i});
    xlabel('Year');
    ylabel('Detrended, Demeaned');
    grid on;
end

sgtitle('Detrended and Demeaned Series with NBER Recession Bars');