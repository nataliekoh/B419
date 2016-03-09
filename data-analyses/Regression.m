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

%% Quick test for VBPCA vs ALSPCA using subset of data from 2010
% Subset the data so that there are fewer NaN values to impute
% col 5-14, 18, 23, 26, 27 in data_2010 from analyses.m
% col 1-10, 14, 19, 22, 23 in year_2010 from this script

% Standardize the data 
zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')),...
    std(x, 'omitnan'));
Z2010 = zscore_xnan(year_2010); 
 
% Now normalize data 
normcol_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, nanmin(x)),...
    bsxfun(@minus, nanmax(x), nanmin(x)));
norm2010 = normcol_xnan(Z2010); 

% Subset the data to match with Vincent's 
subset_2010 = norm2010(:, [2:11, 15, 20, 23:24]); 

% Run the reconstruction using VB
opts = struct( 'maxiters', 70,...
               'algorithm', 'vb',...
               'uniquesv', 0,...
               'cfstop', [ 100 0 0 ],...
               'rotate2pca', 1,...
               'minangle', 0 );
[A10, S10, Mu10, V10, cv10, hp10, lc10] = pca_full(subset_2010, 7, opts);

% Reconstruct the data
Xrec2010 = repmat(Mu10, 1, size(subset_2010,2)) + A10*S10; 
[C, S, L] = pca(Xrec2010);

%% Use the ALS on 2010 (start with year 2010) 

% Subsetting the data to match Vincent's
raw_2010 = xyear_2010(:, [1:10, 14, 19, 22:23]);
headers_subset = xheaders_year([1:10, 14, 19, 22:23]);
[coeff10,score10,latent10,~] = pca(raw_2010, 'Algorithm', 'als');

fig1 = figure; 
stairs(cumsum(latent10)/sum(latent10));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - ALS algorithm on Year 2010');
print(fig1, 'img/PCvariance_2010', '-depsc')

fig2 = figure;
biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',headers_subset);

fig3 = figure;
plot(1:14,coeff10(:,1:6),'-');
xlabel('Variable');
ylabel('PCA Loading');
title('PCA Loadings on Risk Factors for 2010');
legend({'1st PC' '2nd PC' '3rd PC'  ...
	'4th PC', '5th PC', '6th PC'},'location','NW');
print(fig3, 'img/PCloadings_2010', '-depsc');

%% Use the ALS on 2007

raw_2007 = xyear_2007(:, [1:10, 14, 19, 22:23]);
[coeff7,score7,latent7,~] = pca(raw_2007, ...
    'Algorithm', 'als');

fig1 = figure; 
stairs(cumsum(latent7)/sum(latent7));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - ALS algorithm on Year 2007'); 
print(fig1, 'img/PCvariance_2007', '-depsc');

fig2 = figure;
biplot(coeff7(:,1:2),'scores',score7(:,1:2),'varlabels',headers_subset);

fig3 = figure;
plot(1:14,coeff7(:,1:6),'-');
xlabel('Variable');
ylabel('PCA Loading');
title('PCA Loadings on Risk Factors for 2007');
legend({'1st PC' '2nd PC' '3rd PC'  ...
	'4th PC', '5th PC', '6th PC'},'location','NW');
print(fig3, 'img/PCloadings_2007', '-depsc');

%% Use the ALS on 2008

raw_2008 = xyear_2008(:, [1:10, 14, 19, 22:23]);
[coeff8,score8,latent8,~] = pca(raw_2008, ...
    'Algorithm', 'als');

fig1 = figure; 
stairs(cumsum(latent8)/sum(latent8));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - ALS algorithm on Year 2008');
print(fig1, 'img/PCvariance_2008', '-depsc');

fig2 = figure;
biplot(coeff8(:,1:2),'scores',score8(:,1:2),'varlabels',headers_subset);

fig3 = figure;
plot(1:14,coeff8(:,1:6),'-');
xlabel('Variable');
ylabel('PCA Loading');
title('PCA Loadings on Risk Factors for 2008');
legend({'1st PC' '2nd PC' '3rd PC'  ...
	'4th PC', '5th PC', '6th PC'},'location','NW');
print(fig3, 'img/PCloadings_2008', '-depsc');

%% Use the ALS on 2009

raw_2009 = xyear_2009(:, [1:10, 14, 19, 22:23]);
[coeff9,score9,latent9,~] = pca(raw_2009, ...
    'Algorithm', 'als');

fig1 = figure; 
stairs(cumsum(latent9)/sum(latent9));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - ALS algorithm on Year 2009'); 
print(fig1, 'img/PCvariance_2009', '-depsc');

fig2 = figure;
biplot(coeff9(:,1:2),'scores',score9(:,1:2),'varlabels',headers_subset);

fig3 = figure;
plot(1:14,coeff9(:,1:6),'-');
xlabel('Variable');
ylabel('PCA Loading');
title('PCA Loadings on Risk Factors for 2009');
legend({'1st PC' '2nd PC' '3rd PC'  ...
	'4th PC', '5th PC', '6th PC'},'location','NW');
print(fig3, 'img/PCloadings_2009', '-depsc');

%% Use the ALS on 2010

raw_2010 = xyear_2010(:, [1:10, 14, 19, 22:23]);
[coeff10,score10,latent10,~] = pca(raw_2010, ...
    'Algorithm', 'als');

fig1 = figure; 
stairs(cumsum(latent10)/sum(latent10));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - ALS algorithm on Year 2010'); 
print(fig1, 'img/PCvariance_2010', '-depsc');

fig2 = figure;
biplot(coeff10(:,1:2),'scores',score10(:,1:2),'varlabels',headers_subset);

fig3 = figure;
plot(1:14,coeff10(:,1:6),'-');
xlabel('Variable');
ylabel('PCA Loading');
title('PCA Loadings on Risk Factors for 2010');
legend({'1st PC' '2nd PC' '3rd PC'  ...
	'4th PC', '5th PC', '6th PC'},'location','NW');
print(fig3, 'img/PCloadings_2010', '-depsc');

%% Use the ALS on 2011

raw_2011 = xyear_2011(:, [1:10, 14, 19, 22:23]);
[coeff11,score11,latent11,~] = pca(raw_2011, ...
    'Algorithm', 'als');

fig1 = figure; 
stairs(cumsum(latent11)/sum(latent11));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - ALS algorithm on Year 2011'); 
print(fig1, 'img/PCvariance_2011', '-depsc');

fig2 = figure;
biplot(coeff11(:,1:2),'scores',score11(:,1:2),'varlabels',headers_subset);

fig3 = figure;
plot(1:14,coeff11(:,1:6),'-');
xlabel('Variable');
ylabel('PCA Loading');
title('PCA Loadings on Risk Factors for 2011');
legend({'1st PC' '2nd PC' '3rd PC'  ...
	'4th PC', '5th PC', '6th PC'},'location','NW');
print(fig3, 'img/PCloadings_2011', '-depsc');

%% Use the ALS on 2012

raw_2012 = xyear_2012(:, [1:10, 14, 19, 22:23]);
[coeff12,score12,latent12,~] = pca(raw_2012, ...
    'Algorithm', 'als');

fig1 = figure; 
stairs(cumsum(latent12)/sum(latent12));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - ALS algorithm on Year 2012'); 
print(fig1, 'img/PCvariance_2012', '-depsc');

fig2 = figure;
biplot(coeff12(:,1:2),'scores',score12(:,1:2),'varlabels',headers_subset);

fig3 = figure;
plot(1:14,coeff12(:,1:6),'-');
xlabel('Variable');
ylabel('PCA Loading');
title('PCA Loadings on Risk Factors for 2012');
legend({'1st PC' '2nd PC' '3rd PC'  ...
	'4th PC', '5th PC', '6th PC'},'location','NW');
print(fig3, 'img/PCloadings_2012', '-depsc');

%% Extract the most important risk factors from PCA output

risk7 = getPCvariables(headers_subset, coeff7, 6)
risk8 = getPCvariables(headers_subset, coeff8, 6)
risk9 = getPCvariables(headers_subset, coeff9, 6)
risk10 = getPCvariables(headers_subset, coeff10, 6)
risk11 = getPCvariables(headers_subset, coeff11, 6)
risk12 = getPCvariables(headers_subset, coeff12, 6)

%% Regression Part I: Test using 2010 raw data
% Here we will be using the 2010 PCA output just to extract the most important factors

% Find the factors that are contributing to the first 6 PCs
trans_header = headers_subset';
factor10_tmp = table(trans_header,coeff10);
factor12_tmp = table(trans_header,coeff12);

% If we look at factor_tmp, we can see the risk factors that are 
% weighted most heavily in each progressive PC:
% 1) {'No_RRorMDR_TBres'} - col 6
% 2) {'Per_TBwithHIV_CPT'} - col 3
% 3) {'Per_new_RRorMDR_TBres'} - col 5
% 4) {'Healthposts_density'} - col 11
% 5) {'Per_TBwithHIV_ART'} - col 4
% 6) {'Per_TBwithHIV'} - col 1

% Now we are going to regress these variables against per_TB
    % First, create the predictor matrix
    predictors_2010 = raw_2010(:, [6 3 5 11 4 1]);
    % Second, create the response matrix, i.e. TB prev %
    response_2010 = year_2010(:, 1);

% Center and scale data before feeding into regression
predictors_2010 = zscore_xnan(predictors_2010); 
predictors_2010 = normcol_xnan(predictors_2010);
response_2010 = zscore_xnan(response_2010); 
response_2010 = normcol_xnan(response_2010); 
    
% Subset data for testing and training
test_frac = 0.2;
perm_2010 = randperm(size(predictors_2010, 1));
test_predictors_2010 = predictors_2010(1:floor(test_frac*size(predictors_2010,1)),:);
train_predictors_2010 = predictors_2010(ceil(test_frac*size(predictors_2010,1)):end,:);
test_response_2010 = response_2010(1:floor(test_frac*size(predictors_2010,1)),:);
train_response_2010 = response_2010(ceil(test_frac*size(predictors_2010,1)):end,:);

% linear fit
mdl_2010 = fitlm(train_predictors_2010, train_response_2010);
anova(mdl_2010)
anova(mdl_2010, 'summary')

% polynomial fit
mdlquad_2010 = fitlm(train_predictors_2010, train_response_2010,...
    'purequadratic');
anova(mdlquad_2010)

% Check residuals
plotResiduals(mdl_2010);
outl = find(mdl_2010.Residuals.Raw > 0.5);
mdl2_2010 = fitlm(train_predictors_2010, train_response_2010,...
    'Exclude',outl);
plotResiduals(mdl2_2010);
plotResiduals(mdl2_2010,'fitted');

% Test our fit on training data - does not work with NaN values!!
ypred = predict(mdl_2010, test_predictors_2010);
plot(train_predictors_2010, train_response_2010, 'o',...
    test_predictors_2010, ypred, 'x')
legend('Data','Predictions'); 

%% Regression Part II: Use PCA scores instead for linear regression

% Create response matrices for all years (i.e. TB prev %) 
response_2007 = year_2007(:, 1);
response_2008 = year_2008(:, 1);
response_2009 = year_2009(:, 1);
response_2010 = year_2010(:, 1);
response_2011 = year_2011(:, 1);
response_2012 = year_2012(:, 1);

% extract scores from the first 6 PCs
PC_6_2007 = score7(:,1:6);
PC_6_2008 = score8(:,1:6);
PC_6_2009 = score9(:,1:6);
PC_6_2010 = score10(:,1:6);
PC_6_2011 = score11(:,1:6);
PC_6_2012 = score12(:,1:6);

% subset data for testing and training
test_frac = 0.2;
perm_PC2010 = randperm(size(PC_6_2010, 1));
test_pred_PC2010 = PC_6_2010(1:floor(test_frac*size(PC_6_2010,1)),:);
train_pred_PC2010 = PC_6_2010(ceil(test_frac*size(PC_6_2010,1)):end,:);
test_resp_PC2010 = response_2010(1:floor(test_frac*size(response_2010,1)),:);
train_resp_PC2010 = response_2010(ceil(test_frac*size(response_2010,1)):end,:);

% fit linear models
mdlPC_2007 = fitlm(PC_6_2007, response_2007)
mdlPC_2008 = fitlm(PC_6_2008, response_2008)
mdlPC_2009 = fitlm(PC_6_2009, response_2009)
mdlPC_2010 = fitlm(PC_6_2010, response_2010)
mdlPC_2011 = fitlm(PC_6_2011, response_2011)
mdlPC_2012 = fitlm(PC_6_2011, response_2012)
mdl = fitlm(train_pred_PC2010, train_resp_PC2010)

% print anova results
anova(mdlPC_2007, 'summary')
anova(mdlPC_2008, 'summary')
anova(mdlPC_2009, 'summary')
anova(mdlPC_2010, 'summary')
anova(mdlPC_2011, 'summary')
anova(mdlPC_2012, 'summary')

%% Crossval Prediction

ypred = predict(mdl, test_pred_PC2010);
nanmean(ypred) - nanmean(test_resp_PC2010)

%% Exclude outliers and refit models
plotResiduals(mdlPC_2010);
outlPC = find(mdlPC_2010.Residuals.Raw > 0.6);
mdlPC2_2010 = fitlm(PC_6, response_2010,...
    'Exclude',outlPC);
figure;
plotResiduals(mdlPC2_2010);
figure;
plotResiduals(mdlPC2_2010,'fitted');

%% Try quadratic model fitting

mdlPCquad_2010 = fitlm(PC_6_2010, response_2010,...
    'purequadratic')
mdlPCquad_2012 = fitlm(PC_6_2012, response_2012,...
    'purequadratic')

%% Stepwise regression

mdlPCstep_2010 = stepwiselm(PC_6_2010,response_2010,'linear')

%% --- Working section ---
% betaPCR = regress(response_2010, score(:,1:6));
% betaPCR = coeff(:,1:6)*betaPCR;
% betaPCR = [mean(y) - mean(X)*betaPCR; betaPCR];
% yfitPCR = [ones(n,1) X]*betaPCR;
% RSS_PCR = sum((y-yfitPCR).^2);
% rsquaredPCR = 1 - RSS_PCR/TSS






