%% Load datasets and do some quick preliminary sorting

% To-do: store units in table properties? 
% d_tbhiv: leave as is 
% d_tbresis: extract col 1,2,4,6
% d_tbtreat: extract col 1,2,3,4
% d_tobacco: extract col 1,2,3,4
% d_alcohol: transpose and extract only data for 'All types' (across years)
% d_diabetes: extract col 1,2,3
% d_healthcare: extract col 1:5
% d_water: extract col 1,2,5,8
% d_life: extract col 1:3
% d_popden: transpose and extract data only from year 2010; col 1, 6
% d_chnutri: extract col 1,2,5

% Load as tables
d_water = readtable('data/WHO_water+sani+hyg.csv');
    d_water = d_water(2:size(d_water,1),[1,2,5,8]); 
    d_water.Properties.VariableNames = {'Country' 'Year' 'Per_cleanwater_access' 'Per_sanitation_access'};
d_tobacco = readtable('data/WHO_tobacco.csv');
    d_tobacco = d_tobacco(2:size(d_tobacco),1:4);
    d_tobacco.Properties.VariableNames = {'Country' 'Year' 'Prevalence_female_smoking' 'Prevalence_male_smoking'};
d_tbtreat = readtable('data/WHO_tbtreatments.csv');
    d_tbtreat = d_tbtreat(:,1:4);
    d_tbtreat.Properties.VariableNames = {'Country' 'Year' 'Per_success_treat_oldTBcases' 'Per_success_treat_newTBcases'};
d_tbresis = readtable('data/WHO_tbresistance.csv');
    d_tbresis = d_tbresis(:,[1,2,4,6]);
    d_tbresis.Properties.VariableNames = {'Country' 'Year' 'Per_new_RRorMDR_TBres' 'No_RRorMDR_TBres'};
d_tbhiv = readtable('data/WHO_tb+hiv.csv');
    d_tbhiv.Properties.VariableNames = {'Country' 'Year' 'No_TB' 'Per_TBwithHIV'...
        'Per_TBwithHIV_tested' 'Per_TBwithHIV_CPT' 'Per_TBwithHIV_ART'};
d_popden = readtable('data/WHO_popdensity.csv');
    d_popden = d_popden(2:size(d_popden,1),[1,6]); % Needs to be transposed later
    d_popden.Properties.VariableNames = {'Country' 'Per_urban_pop_2010'};
d_lifestats = readtable('data/WHO_life.csv');
    d_lifestats = d_lifestats(:,1:3);
    d_lifestats.Properties.VariableNames = {'Country' 'Year' 'Life_expectancy'};
d_healthcare = readtable('data/WHO_healthcare.csv');
    d_healthcare = d_healthcare(:,1:5);
    d_healthcare.Properties.VariableNames = {'Country' 'Year' 'Hospital_density' 'Healthposts_density' 'Heathcentres_density'};
d_diabetes = readtable('data/WHO_diabetes.csv','ReadVariableNames', false);
    d_diabetes = d_diabetes(4:size(d_diabetes,1),1:3);
    d_diabetes.Properties.VariableNames = {'Country' 'Year' 'Diabetes_mortality_rate'};
d_chnutri = readtable('data/WHO_childnutrition.csv');
    d_chnutri = d_chnutri(:,[1,2,5]);
    d_chnutri.Properties.VariableNames = {'Country' 'Year' 'Per_childbelow5_underweight'};
d_alcohol = readtable('data/WHO_alcohol.csv'); 
    d_alcohol = d_alcohol(2:size(d_alcohol,1),[1,3:size(d_alcohol,2)]); % Needs to be transposed later
    d_alcohol.Properties.VariableNames = {'Country' 'Type' 'al2013' 'al2012' 'al2011' 'al2010'...
    'al2009' 'al2008' 'al2007' 'al2006' 'al2005' 'al2004' 'al2003' 'al2002' 'al2001' 'al2000'};

% Load as arrays (for easier manipulation later) 
da_water = table2array(d_water);
da_tobacco = table2array(d_tobacco);
da_tbtreat = table2array(d_tbtreat);
da_tbresis = table2array(d_tbresis);
da_tbhiv = table2array(d_tbhiv);
da_popden = table2array(d_popden);
da_lifestats = table2array(d_lifestats);
da_healthcare = table2array(d_healthcare);
da_diabetes = table2array(d_diabetes);
da_chnutri = table2array(d_chnutri);
da_alcohol = table2array(d_alcohol);

%% Clean up and subset datasets

% d_popden cleanup here
year_2010 = repmat(2010,size(d_popden,1),1);
popden_table = table(table2array(d_popden(:,1)), year_2010, table2array(d_popden(:,2)));
popden_table.Properties.VariableNames = {'Country' 'Year' 'Per_urban_pop_2010'};

% d_alcohol cleanup here
tmp = zeros(size(da_alcohol, 1), 1); 
for i = 1:size(da_alcohol, 1); tmp(i) = strncmp(da_alcohol(i,2), ' All types', 10); end
new_dal = da_alcohol(find(tmp==1),:);
new_d_alcohol = cell2table(new_dal);
tmp2 = new_dal'; 
tmp3 = tmp2(3:size(tmp2,1),:); 

countries = tmp2(1,:)';
years = (2013:-1:2000)';
repyears = repmat(years,size(countries,1),1);
repcountries = cell(size(years,1)*size(countries,1),1);
a = 1:14:size(repcountries, 1)+14; 

for i = 1:size(countries, 1); 
    repcountries(a(i):a(i+1)-1,1) = repmat(countries(i,1),14,1);
end;

repyears = repmat(years,size(countries,1),1);
al_table = table(char(repcountries), repyears, reshape(test, [], 1));
al_table.Properties.VariableNames = {'Country' 'Year' 'Alcohol_consumption_liters_per_cap'};

%% Write edited tables out to csv files

writetable(d_tbhiv, 'd_tbhiv.csv');
writetable(d_tbresis, 'd_tbresis.csv');
writetable(d_tbtreat, 'd_tbtreat.csv');
writetable(d_tobacco, 'd_tobacco.csv');
writetable(al_table, 'd_alcohol.csv');
writetable(d_diabetes, 'd_diabetes.csv');
writetable(d_healthcare, 'd_healthcare.csv');
writetable(d_water, 'd_water.csv');
writetable(d_lifestats, 'd_lifestats.csv');
writetable(popden_table, 'd_popden.csv');
writetable(d_chnutri, 'd_chnutri.csv');

%% Merge - It's a pain to do this in matlab so I'll be using R instead.
% -------Ignore this section ---------

% Sort order of datasets (for merging): 
% tb_hiv, tbresistance, tbtreatments, tobacco, alcohol, ...
% diabetes, healthcare, water+sani+hyg, life, popdensity, childnutrition

headers = {'Country', 'Year', 'No_TB', 'Per_TBwithHIV',...
    'Per_TBwithHIV_tested', 'Per_TBwithHIV_CPT', 'Per_TBwithHIV_ART',...
    'Per_new_RRorMDR_TBres', 'No_RRorMDR_TBres', 'Per_success_treat_oldTBcases',...
    'Per_success_treat_newTBcases', 'Prevalence_female_smoking', 'Prevalence_male_smoking',...
    'Alcohol_consumption_liters_per_cap', 'Diabetes_mortality_rate',...
    'Hospital_density', 'Healthposts_density', 'Heathcentres_density',...
    'Per_cleanwater_access', 'Per_sanitation_access',...
    'Life_expectancy', 'Per_urban_pop_2010', 'Per_childbelow5_underweight'}

merged_data1 = outerjoin(d_tbhiv, d_tbresis, 'MergeKeys',true);
merged_data2 = outerjoin(merged_data1, d_tbtreat, 'MergeKeys',true);
merged_data3 = outerjoin(merged_data2, d_tobacco, 'MergeKeys',true);
merged_data4 = outerjoin(merged_data3, al_table, 'MergeKeys',true);
% merged_data5 = outerjoin(merged_data4, d_diabetes, 'MergeKeys',true);
% merged_data6 = outerjoin(merged_data5, d_healthcare, 'MergeKeys',true);
% merged_data7 = outerjoin(merged_data6, d_water, 'MergeKeys',true);
% merged_data8 = outerjoin(merged_data7, d_life, 'MergeKeys',true);
% merged_data9 = outerjoin(merged_data8, popden_table, 'MergeKeys',true);
% merged_data10 = outerjoin(merged_data9, d_chnutri, 'MergeKeys',true);
% 
% 
