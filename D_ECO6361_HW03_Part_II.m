% Loading new data
data = readtable('Alternate.csv');
time = datetime(data.DATE, 'InputFormat', 'yyyy-MM-dd');
lambdaa = 1600;


% Original hours worked and inflation with same definitions as Part I
hours_worked = log((data.PRS85006023 .* data.CE16OV / 100) ./ data.CLF16OV) * 100;
inflation = [NaN; log(data.GDPDEF(2:end) ./ data.GDPDEF(1:end-1)) * 100];

% Alternative variables: average weekly hours of production and nonsupervisory employees across private industries for hours worked (AWHNONAG), 
% Consumer Price Index for All Urban Consumers for inflation
% Used simillar equations as above
alt_labor_hours = log((data.AWHNONAG .* data.CE16OV / 100) ./ data.CLF16OV) * 100;
alt_inflation = log(data.CPIAUCSL) * 100;

% Applying HP filter and demean each variable
[~, hours_demeaned] = hpfilter(hours_worked, lambdaa);
hours_demeaned = hours_demeaned - mean(hours_demeaned, 'omitnan');

[~, inflation_demeaned] = hpfilter(inflation, lambdaa);
inflation_demeaned = inflation_demeaned - mean(inflation_demeaned, 'omitnan');

[~, alt_labor_hours_demeaned] = hpfilter(alt_labor_hours, lambdaa);
alt_labor_hours_demeaned = alt_labor_hours_demeaned - mean(alt_labor_hours_demeaned, 'omitnan');

[~, alt_inflation_demeaned] = hpfilter(alt_inflation, lambdaa);
alt_inflation_demeaned = alt_inflation_demeaned - mean(alt_inflation_demeaned, 'omitnan');

% Ensuring all series have the same length as time, this is introduced
% becuase I recurrently came up with a dimension issue. Also, AWHNONAG
% starts from 1964, thus I have truncated

min_length = min([length(time), length(hours_demeaned), length(inflation_demeaned), ...
                  length(alt_labor_hours_demeaned), length(alt_inflation_demeaned)]);

% Truncate all series to min_length
time = time(1:min_length);
hours_demeaned = hours_demeaned(1:min_length);
inflation_demeaned = inflation_demeaned(1:min_length);
alt_labor_hours_demeaned = alt_labor_hours_demeaned(1:min_length);
alt_inflation_demeaned = alt_inflation_demeaned(1:min_length);

% Plotting comparison for original vs alternative hours worked and inflation
figure;

% Original vs Alternative Hours Worked
subplot(2, 1, 1);
plot(time, hours_demeaned, 'b', 'LineWidth', 1.5); hold on;
plot(time, alt_labor_hours_demeaned, 'r--', 'LineWidth', 1.5);
title('Original vs. Alternative Hours Worked');
xlabel('Year');
ylabel('Detrended, Demeaned');
legend('Original Hours Worked', 'Alternative (AWHNONAG)', 'Location', 'Best');
grid on;

% Original vs Alternative Inflation
subplot(2, 1, 2);
plot(time, inflation_demeaned, 'b', 'LineWidth', 1.5); hold on;
plot(time, alt_inflation_demeaned, 'r--', 'LineWidth', 1.5);
title('Original vs. Alternative Inflation');
xlabel('Year');
ylabel('Detrended, Demeaned');
legend('Original Inflation', 'Alternative (CPIAUCSL)', 'Location', 'Best');
grid on;

sgtitle('Comparison of Original and Alternative Variables for Hours Worked and Inflation');