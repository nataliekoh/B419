%% Plot predictor variables against TBprev for each year

% 2007
fig1 = figure; 
for i = 1:size(xsyear_2007,2),
    subplot(6,4,i);
    scatter(xsyear_2007(:,i), response_2007)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2007: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig1, '../img/2007_prelimplot', '-djpeg');

% 2008
fig2 = figure; 
for i = 1:size(xsyear_2008,2),
    subplot(6,4,i);
    scatter(xsyear_2008(:,i), response_2008)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2008: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig2, '../img/2008_prelimplot', '-djpeg');

% 2009
fig3 = figure; 
for i = 1:size(xsyear_2009,2),
    subplot(6,4,i);
    scatter(xsyear_2009(:,i), response_2009)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2009: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig3, '../img/2009_prelimplot', '-djpeg');

% 2010
fig4 = figure; 
for i = 1:size(xsyear_2010,2),
    subplot(6,4,i);
    scatter(xsyear_2010(:,i), response_2010)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2010: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig4, '../img/2010_prelimplot', '-djpeg');

% 2011
fig5 = figure; 
for i = 1:size(xsyear_2011,2),
    subplot(6,4,i);
    scatter(xsyear_2011(:,i), response_2011)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2011: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig5, '../img/2011_prelimplot', '-djpeg');

% 2012
fig6 = figure; 
for i = 1:size(xsyear_2012,2),
    subplot(6,4,i);
    scatter(xsyear_2012(:,i), response_2012)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2012: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig6, '../img/2012_prelimplot', '-djpeg');

% 2013
fig7 = figure; 
for i = 1:size(xsyear_2013,2),
    subplot(6,4,i);
    scatter(xsyear_2013(:,i), response_2013)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2013: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig7, '../img/2013_prelimplot', '-djpeg');

%% Identify risk factor indexes (non-empty plots) for data analysis

% SEE data-analyses/riskfactorindices.txt 

%% Access sparseness of dataset

data = readtable('data/condensed_with_perTB.csv');
% n = no. of observations
n = size(data, 1);
sum(~isnan(table2array(new_data(:,11))))/n
sum(~isnan(table2array(new_data(:,12))))/n
sum(~isnan(table2array(new_data(:,13))))/n
sum(~isnan(table2array(new_data(:,14))))/n
sum(~isnan(table2array(new_data(:,15))))/n
sum(~isnan(table2array(new_data(:,16))))/n
sum(~isnan(table2array(new_data(:,17))))/n
sum(~isnan(table2array(new_data(:,18))))/n

% total number of data points available
sum(sum(~isnan(table2array(new_data(:,3:size(new_data,2))))))

% total number of data points available as %
sum(sum(~isnan(table2array(new_data(:,3:size(new_data,2))))))/...
    (size(table2array(new_data(:,3:size(new_data,2))),1)*...
    size(table2array(new_data(:,3:size(new_data,2))),2))*100

% 27.06% of our dataset is filled;
% 72.94% of our dataset is made up of missing values!!!!!


