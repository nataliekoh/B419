function [n_pc, l, q, t] = classifications(y, columns)
%classifications - description

%% data load

file = readtable('../data/condensed_with_perTB.csv');

%% pull specific year and modify
year = find_year_data(file, y);
usable_rows = find(~isnan(year(:,1)) == 1);
xyear = year(usable_rows, 2:size(year,2));

zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')),...
    std(x, 'omitnan'));

xsyear = zscore_xnan(xyear);
xssyear = xsyear(:, columns);

[coeff,score,latent,~,~,mu] = pca(xssyear, 'Algorithm', 'als',...
    'Centered',false);

filtered_var_norm = find(cumsum(latent)/sum(latent) >= .95, 1);
n_pc = filtered_var_norm(1);

%% Plotting the projections with prevalence categories
qual_prev = cell(numel(usable_rows), 1);
filtered_rows = year(usable_rows, 1);

% Using the below link on page 25. The prevalence rate of high-burden
% countries in 2014 were averaged to mark the threshold between high and
% low prevalence. The average rate was rounded to 0.2% of total population
% having TB.
% http://apps.who.int/iris/bitstream/10665/191102/1/9789241565059_eng.pdf?ua=1
for i = 1:numel(filtered_rows)
    if filtered_rows(i, 1) < 0.2
        qual_prev{i, 1} = 'Not-High';
    elseif filtered_rows(i, 1) >= 0.2
        qual_prev{i, 1} = 'High';
    else
        qual_prev{i, 1} = 'No data available';
    end
end;

%% lda
% lda classifier cross validation runs
test_fraction = 0.3;
n_trials = 100;
all_cv_lda = zeros(n_trials, 1);

for i = 1:n_trials
    permuted = randperm(numel(qual_prev));
    test = permuted(1 : floor(numel(qual_prev) * test_fraction));
    train = permuted(ceil((numel(qual_prev) * test_fraction)) : end);

    % lda
    lda = fitcdiscr(score(train, 1:n_pc), qual_prev(train, :));
    label_lda = predict(lda, score(test, 1:n_pc));
    all_cv_lda(i) = mean(strcmp(qual_prev(test, :), label_lda));
end;

l = mean(all_cv_lda);

%% Classification tree classifier

test_fraction = 0.3;
n_trials = 100;
all_cv_tree = zeros(n_trials, 1);

for i = 1:n_trials
    permuted = randperm(numel(qual_prev));
    test = permuted(1 : floor(numel(qual_prev) * test_fraction));
    train = permuted(ceil((numel(qual_prev) * test_fraction)) : end);

    % classification trees
    tree = fitctree(score(train, 1:n_pc), qual_prev(train, :));
    label_tree = predict(tree, score(test, 1:n_pc));
    all_cv_tree(i) = mean(strcmp(qual_prev(test, :), label_tree));
end;

t = mean(all_cv_tree);
%% QDA classifier

% qda classifier cross validation runs
test_fraction = 0.3;
n_trials = 100;
cv_qda = 0;
count = 0;

for i = 1:n_trials
    try
        permuted = randperm(numel(qual_prev));
        test = permuted(1 : floor(numel(qual_prev) * test_fraction));
        train = permuted(ceil((numel(qual_prev) * test_fraction)) : end);

        % qda
        qda = fitcdiscr(score(train, 1:n_pc), qual_prev(train, :), 'DiscrimType', 'quadratic');
        label_qda = predict(lda, score(test, 1:n_pc));
        
        cv_qda = cv_qda + mean(strcmp(qual_prev(test, :), label_qda));
        count = count + 1;
    catch ME
        warning('Error with covariance matrices for QDA');
    end
end;

q = cv_qda / count;

end