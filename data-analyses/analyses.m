%% load the data

% Make sure working directory is set in same folder as this script
file = readtable('../data/all_condensed_v7.csv');

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

%% Data pull
unique_countries = unique(file(:, 2));

% 2005 data
data_2005 = getAllDataForYear(file, 2005);
TBprev_2005 = data_2005{:, 4} ./ (data_2005{:, 25} * 1000) * 100;

% 2010 data
data_2010 = getAllDataForYear(file, 2010);
TBprev_2010 = data_2010{:, 4} ./ (data_2010{:, 25} * 1000) * 100;

%% Histogram of TB prevalence in 2005 and 2010
h_binwidth = 0.025;
figure;
hold on;
histogram(TBprev_2005, 'BinWidth', h_binwidth);
histogram(TBprev_2010, 'BinWidth', h_binwidth);
legend('2005', '2010');
title('TB prevalence in 2005 vs 2010');
xlabel('Percentage of population with TB');
ylabel('Number of countries');

%% Individual plots of all risk factors vs TB prevalence in 2010
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

%% classifications test
cv = zeros(3, 3);
[p, l, q, t] = classifications(2007, [1:8, 20:23]);
cv(1, :) = [l, q, t];

[p, l, q, t] = classifications(2010, [1:8, 14:15, 20:23]);
cv(2, :) = [l, q, t];

[p, l, q, t] = classifications(2012, [1:8, 12, 18, 20:23]);
cv(3, :) = [l, q, t];

%% Plot cross validated accuracies
figure;
bar([2007,2010,2012], cv);
legend('LDA', 'QDA', 'fitctree', 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Different classification methods on classifying high prevalence');
xlabel('Year');
ylabel('% cross-validated accuracy');

%% OLD METHODS - Further cleaning of 2010 data in order to be used with PCA
% clean_indices = find(~isnan(TBprev_2010) == 1);
% clean_2010 = TBprev_2010(clean_indices);
% clean_2010_countries = unique_countries{clean_indices, 1};
% clean_2010_data = data_2010{clean_indices, [5:14,18,19,23,24,26,27]};
% 
% [coeff, score, var] = pca(clean_2010_data, 'algorithm', 'als');
% stairs(cumsum(var) / sum(var));
% 
% %% Plotting the projections with prevalence categories
% qual_prev = cell(numel(clean_2010), 1);
% 
% % Using the below link on page 25. The prevalence rate of high-burden
% % countries in 2014 were averaged to mark the threshold between high and
% % low prevalence. The average rate was rounded to 0.2% of total population
% % having TB.
% % http://apps.who.int/iris/bitstream/10665/191102/1/9789241565059_eng.pdf?ua=1
% for i = 1:numel(clean_2010)
%     if clean_2010(i, 1) < 0.2
%         qual_prev{i, 1} = 'Not-High';
%     elseif clean_2010(i, 1) >= 0.2
%         qual_prev{i, 1} = 'High';
%     end
% end;
% 
% % Quick plot to see what the first 2 PC's look like.
% figure;
% gscatter(score(:,1), score(:,2), qual_prev);
% title('PCA of Risk Factors with Raw High Prevalence Values');
% xlabel('PC 1');
% ylabel('PC 2');
% 
% %% LDA classifier
% 
% % lda classifier cross validation runs
% test_fraction = 0.3;
% n_trials = 100;
% all_cv_lda = zeros(n_trials, 1);
% 
% for i = 1:n_trials
%     permuted = randperm(numel(qual_prev));
%     test = permuted(1 : floor(numel(qual_prev) * test_fraction));
%     train = permuted(ceil((numel(qual_prev) * test_fraction)) : end);
% 
%     % lda
%     lda = fitcdiscr(score(train, 1:6), qual_prev(train, :));
%     label_lda = predict(lda, score(test, 1:6));
%     all_cv_lda(i) = mean(strcmp(qual_prev(test, :), label_lda));
% end;
% 
% disp(mean(all_cv_lda));
% 
% %% Classification tree classifier
% 
% test_fraction = 0.3;
% n_trials = 100;
% all_cv_tree = zeros(n_trials, 1);
% 
% for i = 1:n_trials
%     permuted = randperm(numel(qual_prev));
%     test = permuted(1 : floor(numel(qual_prev) * test_fraction));
%     train = permuted(ceil((numel(qual_prev) * test_fraction)) : end);
% 
%     % classification trees
%     tree = fitctree(score(train, 1:6), qual_prev(train, :));
%     label_tree = predict(tree, score(test, 1:6));
%     all_cv_tree(i) = mean(strcmp(qual_prev(test, :), label_tree));
% end;
% 
% disp(mean(all_cv_tree));
% %% QDA classifier
% 
% % qda classifier cross validation runs
% test_fraction = 0.3;
% n_trials = 100;
% all_cv_qda = zeros(n_trials, 1);
% 
% for i = 1:n_trials
%     permuted = randperm(numel(qual_prev));
%     test = permuted(1 : floor(numel(qual_prev) * test_fraction));
%     train = permuted(ceil((numel(qual_prev) * test_fraction)) : end);
% 
%     % qda
%     qda = fitcdiscr(score(train, 1:6), qual_prev(train, :), 'DiscrimType', 'quadratic');
%     label_qda = predict(lda, score(test, 1:6));
%     all_cv_qda(i) = mean(strcmp(qual_prev(test, :), label_qda));
% end;
% 
% disp(mean(all_cv_qda));
% 
