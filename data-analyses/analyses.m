%% load the data

% Make sure working directory is set as B419 folder.
file = readtable('data/all_condensed_v3.csv');

%% basic descriptions
% col = desc
% 1 = year
% 2 = # of TB
% 3 = % of TB patients with known HIV status
% 4 = % of TB patients (tested) that are HIV-positive
% 5 = % of HIV + TB patients on CPT
% 6 = % of HIV + TB patients on ART


% todo - merge cote d'ivoire
unique_countries = numel(unique(file(:, 2)));

% view afghanistan # of TB
afghan_i = find(strcmp(file{:, 1}, 'Afghanistan') == 1);

afghan_years = file{afghan_i, 2};
afghan_nTB = file{afghan_i, 3};

figure;
plot(afghan_years, afghan_nTB);