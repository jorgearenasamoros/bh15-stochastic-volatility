% READ_DATA  Loads the data and builds the estimation sample (script).
% Adapted from readData1a.m in the Baumeister and Hamilton (2015) replication files.
%
% data/labor_data.csv, quarterly, 1947:Q1-2026:Q2:
%   col 1 = date (decimal year)
%   col 2 = real hourly compensation, nonfarm business sector (FRED: COMPRNFB)
%   col 3 = total nonfarm payroll employment, last month of the quarter (FRED: PAYEMS)
% Estimation sample: 1970:Q1-2019:Q4.

nlags = 8;                                           % lags in the VAR
rootdir = fileparts(mfilename('fullpath'));
labor_data = load(fullfile(rootdir,'data','labor_data.csv'));
wage       = 100*(log(labor_data(2:end,2)) - log(labor_data(1:end-1,2)));
employment = 100*(log(labor_data(2:end,3)) - log(labor_data(1:end-1,3)));

varnames   = {' wage'; ' employment'};
shocknames = {' demand'; ' supply'};
yall  = [wage employment];
dates = labor_data(2:end,1);                         % date of each row of yall
tstart = find(abs(dates-1970)<1e-6);                 % 1970:Q1
tend   = find(abs(dates-2019.75)<1e-6);              % 2019:Q4
YY = yall(tstart:tend,:);
time = (1970:0.25:dates(tend))';
