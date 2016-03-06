%% load the data

% Make sure working directory is set as B419 folder, NOT the same folder 
% as where this script is.
cd '..';
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
legend('2005', '2010');
title('TB prevalence in 2005 vs 2010');
xlabel('Percentage of population with TB');
ylabel('Number of countries');

%% plot testing
count = 1;
xlabels = readtable('subplot_labels.csv');
for i = 5:27,
    if mod(count, 7) == 0 || count == 1
        figure;
        count = 1;
    end;
    subplot(3, 2, count);
    count = count + 1;
    plot(data_2010{:, i}, TBprev_2010, '.');
    ylabel('% population with TB');
    xlabel(xlabels{i, 1});
end;

%% pca
% [coeff,score,pcvar, mu] = ppca(data_2010{:, [5:14,18,23, 26,27]}, 3);

[coeff, score, var] = pca(data_2010{:, [5:14,18,23, 26,27]}, 'algorithm', 'als');
stairs(cumsum(var) / sum(var));

%% plot projection onto pc
qual_prev = cell(numel(TBprev_2010), 1);

% http://www.tbfacts.org/tb-statistics/
for i = 1:numel(TBprev_2010)
    if TBprev_2010(i, 1) < 0.1
        qual_prev{i, 1} = 'Low';
    elseif TBprev_2010(i, 1) < 0.3
        qual_prev{i, 1} = 'Medium';
    elseif TBprev_2010(i, 1) >= 0.3
        qual_prev{i, 1} = 'High';
    else
        qual_prev{i, 1} = 'N/A';
    end
end;

figure;
gscatter(score(:,1), score(:,2), qual_prev);
