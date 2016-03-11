%% Read in the data

data = readtable('../data/all_condensed_v6.csv');
data = data(:,2:size(data,2));
samples = table2array(data(:,3:size(data,2)));
samples_year = table2array(data(:,2:size(data,2)));
headers = data.Properties.VariableNames;

%% Add calculated percent_TB as a column to the dataset

per_TB = table2array(data(:,3))./(table2array(data(:,24))*1000)*100;
new_data = [array2table(per_TB), data(:,4:size(data,2))];
new_data = [data(:,1:2), new_data];
headers_new = new_data.Properties.VariableNames;

%writetable(new_data, '../data/condensed_with_perTB.csv');

%% Get year data

% Here we extract numerical arrays using find_year_data
year_2007 = find_year_data(new_data,2007);
year_2008 = find_year_data(new_data,2008);
year_2009 = find_year_data(new_data,2009);
year_2010 = find_year_data(new_data,2010);
year_2011 = find_year_data(new_data,2011);
year_2012 = find_year_data(new_data,2012);
year_2013 = find_year_data(new_data,2013);

% Now we update each variable so that we are only looking at risk factors
xyear_2007 = year_2007(:, 2:size(year_2007,2));
xyear_2008 = year_2008(:, 2:size(year_2008,2));
xyear_2009 = year_2009(:, 2:size(year_2009,2));
xyear_2010 = year_2010(:, 2:size(year_2010,2));
xyear_2011 = year_2011(:, 2:size(year_2011,2));
xyear_2012 = year_2012(:, 2:size(year_2012,2));
xyear_2013 = year_2013(:, 2:size(year_2013,2));
xheaders_year = headers_new(4:numel(headers_new));

% Standardize predictors for PCA
xsyear_2007 = zscore_xnan(xyear_2007); 
xsyear_2008 = zscore_xnan(xyear_2008); 
xsyear_2009 = zscore_xnan(xyear_2009); 
xsyear_2010 = zscore_xnan(xyear_2010); 
xsyear_2011 = zscore_xnan(xyear_2011); 
xsyear_2012 = zscore_xnan(xyear_2012); 
xsyear_2013 = zscore_xnan(xyear_2013); 

% Get response vectors
response_2007 = year_2007(:, 1);
response_2008 = year_2008(:, 1);
response_2009 = year_2009(:, 1);
response_2010 = year_2010(:, 1);
response_2011 = year_2011(:, 1);
response_2012 = year_2012(:, 1);
response_2013 = year_2013(:, 1);

%% Run PCA-ALS

