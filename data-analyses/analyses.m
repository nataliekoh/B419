%% load the data

% Make sure working directory is set as B419 folder, NOT the same folder 
% as where this script is.
file = readtable('data/all_condensed_v6.csv');

cd 'data-analyses'; % functions in data-analyses folder.

%% basic descriptions
% col = desc
% 3 = year
% 4 = # of TB
% 5 = TB patients with known HIV status (%)
% 6 = TB patients (tested) that are HIV-positive (%)
% 7 = HIV + TB patients on CPT (%)
% 8 = HIV + TB patients on ART (%)
% 5 = Pop (thousands)
% 26 = Pop under 15 (%)
% 27 = Pop over 60 (%)

unique_countries = unique(file(:, 2));

%% TB prevalence calculator
% will probably move into function 

figure;
title('TB prevalence (# of TB patients/population)');
xlabel('Year');
ylabel('% of population with TB');
hold on;

n_countries = numel(unique_countries);
for i = 1:n_countries
    country_name = unique_countries{i, 1};
    
    country_indices = find(strcmp(file{:, 2}, country_name) == 1);
    years = file{country_indices, 3};
    nTB_per_pop = file{country_indices, 4} ./ (file{country_indices, 25} * 1000);
    
    plot(years, nTB_per_pop, '.-');
end;

% legend(unique_countries{1:n_countries, 1});

%% Data pull

data_2005 = getAllDataForYear(file, 2005);
TBprev_2005 = data_2005{:, 4} ./ (data_2005{:, 25} * 1000) * 100;

data_2010 = getAllDataForYear(file, 2010);
TBprev_2010 = data_2010{:, 4} ./ (data_2010{:, 25} * 1000) * 100;

%% Histogram of TB prevalence in 2010
h_binwidth = 0.025;
figure;
hold on;
histogram(TBprev_2005, 'BinWidth', h_binwidth);
histogram(TBprev_2010, 'BinWidth', h_binwidth);
legend('2005', '2010', '2013');
title('TB prevalence in 2005 vs 2010');
xlabel('Percentage of population with TB');
ylabel('Number of countries');

%% plot testing
figure;
suptitle('TB prevalence per population versus different risk factors (TB prevalence/pop in % on y-axis)');
labels = fieldnames(file);
for i = 5:27,
    subplot(4, 6, i-4);
    plot(data_2010{:, i}, TBprev_2010, '.');
    title(['col #', num2str(i)]);
end;

%%
[coeff,score,pcvar, mu] = ppca(data_2010{:, [5:14,18,23:27]}, 1);
% figure;
% scatter(score(:,1), score(:,2), TBprev_2010);