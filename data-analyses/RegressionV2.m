%% !!!THIS VERSION IS MUCH CLEANER THAN Regression.m

%% Read in the data

new_data = readtable('../data/condensed_with_perTB.csv');
headers_new = new_data.Properties.VariableNames;

%% Get year data (pardon the long variable names)

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
zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')),...
    std(x, 'omitnan'));
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

% Subset predictor data to reduce sparseness - see data-analyses/riskfactorindices.txt
xssyear_2007 = xsyear_2007(:, [1:8, 20:23]);
xssyear_2008 = xsyear_2008(:, [1:8, 20:23]);
xssyear_2009 = xsyear_2009(:, [1:8, 20:23]);
xssyear_2010 = xsyear_2010(:, [1:8, 14:15, 20:23]);
xssyear_2011 = xsyear_2011(:, [1:8, 20:23]);
xssyear_2012 = xsyear_2012(:, [1:8, 12, 18, 20:23]);
xssyear_2013 = xsyear_2013(:, [1:8, 13:15, 18, 20:23]);

%% Grab headers for xssyear____

headers_A = xheaders_year([1:8, 20:23]); % for 2007,2008,2009,2011
headers_B = xheaders_year([1:8, 14:15, 20:23]); % for 2010
headers_C = xheaders_year([1:8, 12, 18, 20:23]); % for 2012
headers_D = xheaders_year([1:8, 13:15, 18, 20:23]); % for 2013

%% Run PCA-ALS

[coeff7,score7,latent7,~,~,mu7] = pca(xssyear_2007, 'Algorithm', 'als',...
    'Centered',false);
[coeff8,score8,latent8,~,~,mu8] = pca(xssyear_2008, 'Algorithm', 'als',...
    'Centered',false);
[coeff9,score9,latent9,~,~,mu9] = pca(xssyear_2009, 'Algorithm', 'als',...
    'Centered',false);
[coeff10,score10,latent10,~,~,mu10] = pca(xssyear_2010, 'Algorithm', 'als',...
    'Centered',false);
[coeff11,score11,latent11,~,~,mu12] = pca(xssyear_2011, 'Algorithm', 'als',...
    'Centered',false);
[coeff12,score12,latent12,~,~,mu11] = pca(xssyear_2012, 'Algorithm', 'als',...
    'Centered',false);
[coeff13,score13,latent13,~,~,mu13] = pca(xssyear_2013, 'Algorithm', 'als',...
    'Centered',false);

%% Plot explained variance

figure;
subplot(4,2,1);
stairs(cumsum(latent7)/sum(latent7));
xlabel('No. of Principal Components');
ylabel('% Var Explained')
title('2007');
refline(0,0.95);
subplot(4,2,2);
stairs(cumsum(latent8)/sum(latent8));
xlabel('No. of Principal Components');
ylabel('% Var Explained')
title('2008');
refline(0,0.95);
subplot(4,2,3);
stairs(cumsum(latent9)/sum(latent9));
xlabel('No. of Principal Components');
ylabel('% Variance Explained')
title('2009');
refline(0,0.95);
subplot(4,2,4);
stairs(cumsum(latent10)/sum(latent10));
xlabel('No. of Principal Components');
ylabel('% Var Explained')
title('2010');
refline(0,0.95);
subplot(4,2,5);
stairs(cumsum(latent11)/sum(latent11));
xlabel('No. of Principal Components');
ylabel('% Var Explained')
title('2011');
refline(0,0.95);
subplot(4,2,6);
stairs(cumsum(latent12)/sum(latent12));
xlabel('No. of Principal Components');
ylabel('% Var Explained')
title('2012');
refline(0,0.95);
subplot(4,2,7);
stairs(cumsum(latent13)/sum(latent13));
xlabel('No. of Principal Components');
ylabel('% Var Explained')
title('2013');
refline(0,0.95);
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'PCA-ALS on years 2007-2013 by row', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print('img/PCvariance_ALLyears', '-dpdf');

%% Plot PC1 vs PC2 for all years

figure;
subplot(4,2,1);
scatter(score7(:,1),score7(:,2))
xlabel('PC1');
ylabel('PC2')
title('2007');
subplot(4,2,2);
scatter(score8(:,1),score8(:,2))
xlabel('PC1');
ylabel('PC2')
title('2008');
subplot(4,2,3);
scatter(score9(:,1),score9(:,2))
xlabel('PC1');
ylabel('PC2')
title('2009');
subplot(4,2,4);
scatter(score10(:,1),score10(:,2))
xlabel('PC1');
ylabel('PC2')
title('2010');
subplot(4,2,5);
scatter(score11(:,1),score11(:,2))
xlabel('PC1');
ylabel('PC2')
title('2011');
subplot(4,2,6);
scatter(score12(:,1),score12(:,2))
xlabel('PC1');
ylabel('PC2')
title('2012');
subplot(4,2,7);
scatter(score13(:,1),score13(:,2))
xlabel('PC1');
ylabel('PC2')
title('2013');
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'PCA-ALS on years 2007-2013 by row', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print('img/PC1PC2plot_ALLyears', '-dpdf');

%% Number of PCs to extract from each year 
% see figure from above

% 2007: 8
% 2008: 6
% 2009: 4
% 2010: 6
% 2011: 6
% 2012: 4
% 2013: 9

% extract scores for the number of PCs listed above
PCs_2007 = score7(:,1:8);
PCs_2008 = score8(:,1:6);
PCs_2009 = score9(:,1:4);
PCs_2010 = score10(:,1:6);
PCs_2011 = score11(:,1:6);
PCs_2012 = score12(:,1:4);
PCs_2013 = score13(:,1:9);

%% Extract the most important risk factors from PCA output

risk7 = getPCvariables(headers_A, coeff7, 8)
risk8 = getPCvariables(headers_A, coeff8, 6)
risk9 = getPCvariables(headers_A, coeff9, 4)
risk10 = getPCvariables(headers_B, coeff10, 6)
risk11 = getPCvariables(headers_A, coeff11, 6)
risk12 = getPCvariables(headers_C, coeff12, 4)
risk13 = getPCvariables(headers_D, coeff13, 9)

%% Reconstruct the data for each year

recon_2007 = score7*coeff7';
recon_2008 = score8*coeff8';
recon_2009 = score9*coeff9';
recon_2010 = score10*coeff10';
recon_2011 = score11*coeff11';
recon_2012 = score12*coeff12';
recon_2013 = score13*coeff13';

%% Linear regression using PC scores

% fit linear models
LMPCs_2007 = fitlm(PCs_2007, response_2007)
LMPCs_2008 = fitlm(PCs_2008, response_2008)
LMPCs_2009 = fitlm(PCs_2009, response_2009)
LMPCs_2010 = fitlm(PCs_2010, response_2010)
LMPCs_2011 = fitlm(PCs_2011, response_2011)
LMPCs_2012 = fitlm(PCs_2012, response_2012)
LMPCs_2013 = fitlm(PCs_2013, response_2013)

%% ANOVA on linear regression models

anova(LMPCs_2007, 'summary')
anova(LMPCs_2008, 'summary')
anova(LMPCs_2009, 'summary')
anova(LMPCs_2010, 'summary') % MODEL NOT SIGNIFICANT?
anova(LMPCs_2011, 'summary')
anova(LMPCs_2012, 'summary')
anova(LMPCs_2013, 'summary')

%% Lasso regression using PC scores
% NOTE: THIS TAKES QUITE A WHILE TO RUN

niter = 100;
for i = 1:niter,
    [Blasso_2007 Fitlasso_2007] = lasso(PCs_2007, response_2007, 'CV', 10);
    [Blasso_2008 Fitlasso_2008] = lasso(PCs_2008, response_2008, 'CV', 10);
    [Blasso_2009 Fitlasso_2009] = lasso(PCs_2009, response_2009, 'CV', 10);
    [Blasso_2010 Fitlasso_2010] = lasso(PCs_2010, response_2010, 'CV', 10);
    [Blasso_2011 Fitlasso_2011] = lasso(PCs_2011, response_2011, 'CV', 10);
    [Blasso_2012 Fitlasso_2012] = lasso(PCs_2012, response_2012, 'CV', 10);
    [Blasso_2013 Fitlasso_2013] = lasso(PCs_2013, response_2013, 'CV', 10);
end

%% Find lambda which gives the smallest MSE

Fitlasso_2007.LambdaMinMSE % 0.0026
Fitlasso_2007.IndexMinMSE % 73
Fitlasso_2008.LambdaMinMSE % 0.0029
Fitlasso_2008.IndexMinMSE % 70
Fitlasso_2009.LambdaMinMSE % 0.016
Fitlasso_2009.IndexMinMSE % 67
Fitlasso_2010.LambdaMinMSE % 0.0083
Fitlasso_2010.IndexMinMSE % 92
Fitlasso_2011.LambdaMinMSE % 3.8458e-06
Fitlasso_2011.IndexMinMSE % 1
Fitlasso_2012.LambdaMinMSE % 3.5352e-04
Fitlasso_2012.IndexMinMSE % 51
Fitlasso_2013.LambdaMinMSE % 0.0060
Fitlasso_2013.IndexMinMSE % 76

%% Find betas for each year

Blasso_2007(:,Fitlasso_2007.IndexMinMSE)
Blasso_2008(:,Fitlasso_2008.IndexMinMSE)
Blasso_2009(:,Fitlasso_2009.IndexMinMSE)
Blasso_2010(:,Fitlasso_2010.IndexMinMSE)
Blasso_2011(:,Fitlasso_2011.IndexMinMSE)
Blasso_2012(:,Fitlasso_2012.IndexMinMSE)
Blasso_2013(:,Fitlasso_2013.IndexMinMSE)

%% Lasso plots for PC scores

lassoPlot(Blasso_2010, Fitlasso_2010, 'PlotType', 'CV');

%% Lasso regression using reconstructed data
% NOTE: THIS TAKES A LONGGG TIME TO RUN

niter = 100;
for i = 1:niter,
    [BRlasso_2007 FRlasso_2007] = lasso(recon_2007, response_2007, 'CV', 10);
    [BRlasso_2008 FRlasso_2008] = lasso(recon_2008, response_2008, 'CV', 10);
    [BRlasso_2009 FRlasso_2009] = lasso(recon_2009, response_2009, 'CV', 10);
    [BRlasso_2010 FRlasso_2010] = lasso(recon_2010, response_2010, 'CV', 10);
    [BRlasso_2011 FRlasso_2011] = lasso(recon_2011, response_2011, 'CV', 10);
    [BRlasso_2012 FRlasso_2012] = lasso(recon_2012, response_2012, 'CV', 10);
    [BRlasso_2013 FRlasso_2013] = lasso(recon_2013, response_2013, 'CV', 10);
end

%% Find lambda which gives the smallest MSE

FRlasso_2007.LambdaMinMSE % 0.0069
FRlasso_2007.IndexMinMSE % 79
FRlasso_2008.LambdaMinMSE % 0.0095
FRlasso_2008.IndexMinMSE % 83
FRlasso_2009.LambdaMinMSE % 1.3783e-04
FRlasso_2009.IndexMinMSE % 36
FRlasso_2010.LambdaMinMSE % 0.0013
FRlasso_2010.IndexMinMSE % 60
FRlasso_2011.LambdaMinMSE % 0.0026
FRlasso_2011.IndexMinMSE % 69
FRlasso_2012.LambdaMinMSE % 0.0037
FRlasso_2012.IndexMinMSE % 69
FRlasso_2013.LambdaMinMSE % 0.0110
FRlasso_2013.IndexMinMSE % 81

%% Find betas for each year

BRlasso_2007(:,FRlasso_2007.IndexMinMSE)
BRlasso_2008(:,FRlasso_2008.IndexMinMSE)
BRlasso_2009(:,FRlasso_2009.IndexMinMSE)
BRlasso_2010(:,FRlasso_2010.IndexMinMSE)
BRlasso_2011(:,FRlasso_2011.IndexMinMSE)
BRlasso_2012(:,FRlasso_2012.IndexMinMSE)
BRlasso_2013(:,FRlasso_2013.IndexMinMSE)

%% Matrix visualizations

figure;
imagesc(coeff7); colorbar;
xlabel('Principal Component'); 
ylabel('Variable');
title('Variable Loadings for 2007');
set(gca,'Xtick',1:12, 'Ytick',1:12,...
    'YTickLabel',{'%TBwithHIV',...
    '%TBwithHIVtested',...
    '%TBwithHIV-CPT',...
    '%TBwithHIV-ART',...
    '%newRRorMDRTBres',...
    'NoRRorMDRTBres',...
    '%successtreatoldTBcases',...
    '%successtreatnewTBcases',...
    '%childbelow5underweight',...
    'Pop',...
    'Popunder15',...
    'Popover60'});
print('img/PCLoadingmatrix_2007', '-djpeg');

figure;
imagesc(coeff8); colorbar;
xlabel('Principal Component'); 
ylabel('Variable');
title('Variable Loadings for 2008');
set(gca,'Xtick',1:12, 'Ytick',1:12,...
    'YTickLabel',{'%TBwithHIV',...
    '%TBwithHIVtested',...
    '%TBwithHIV-CPT',...
    '%TBwithHIV-ART',...
    '%newRRorMDRTBres',...
    'NoRRorMDRTBres',...
    '%successtreatoldTBcases',...
    '%successtreatnewTBcases',...
    '%childbelow5underweight',...
    'Pop',...
    'Popunder15',...
    'Popover60'});
print('img/PCLoadingmatrix_2008', '-djpeg');

figure;
imagesc(coeff9); colorbar;
xlabel('Principal Component'); 
ylabel('Variable');
title('Variable Loadings for 2009');
set(gca,'Xtick',1:12, 'Ytick',1:12,...
    'YTickLabel',{'%TBwithHIV',...
    '%TBwithHIVtested',...
    '%TBwithHIV-CPT',...
    '%TBwithHIV-ART',...
    '%newRRorMDRTBres',...
    'NoRRorMDRTBres',...
    '%successtreatoldTBcases',...
    '%successtreatnewTBcases',...
    '%childbelow5underweight',...
    'Pop',...
    'Popunder15',...
    'Popover60'});
print('img/PCLoadingmatrix_2009', '-djpeg');

figure;
imagesc(coeff10); colorbar;
xlabel('Principal Component'); 
ylabel('Variable');
title('Variable Loadings for 2010');
set(gca,'Xtick',1:14, 'Ytick',1:14,...
    'YTickLabel',{'%TBwithHIV',...
    '%TBwithHIVtested',...
    '%TBwithHIV-CPT',...
    '%TBwithHIV-ART',...
    '%newRRorMDRTBres',...
    'NoRRorMDRTBres',...
    '%successtreatoldTBcases',...
    '%successtreatnewTBcases',...
    'Healthpostsdensity',...
    'Heathcentresdensity',...
    '%childbelow5underweight',...
    'Pop',...
    'Popunder15',...
    'Popover60'});
print('img/PCLoadingmatrix_2010', '-djpeg');

figure;
imagesc(coeff11); colorbar;
xlabel('Principal Component'); 
ylabel('Variable');
title('Variable Loadings for 2011');
set(gca,'Xtick',1:12, 'Ytick',1:12,...
    'YTickLabel',{'%TBwithHIV',...
    '%TBwithHIVtested',...
    '%TBwithHIV-CPT',...
    '%TBwithHIV-ART',...
    '%newRRorMDRTBres',...
    'NoRRorMDRTBres',...
    '%successtreatoldTBcases',...
    '%successtreatnewTBcases',...
    '%childbelow5underweight',...
    'Pop',...
    'Popunder15',...
    'Popover60'});
print('img/PCLoadingmatrix_2011', '-djpeg');

figure;
imagesc(coeff12); colorbar;
xlabel('Principal Component'); 
ylabel('Variable');
title('Variable Loadings for 2012');
set(gca,'Xtick',1:14, 'Ytick',1:14,...
    'YTickLabel',{'%TBwithHIV',...
    '%TBwithHIVtested',...
    '%TBwithHIV-CPT',...
    '%TBwithHIV-ART',...
    '%newRRorMDRTBres',...
    'NoRRorMDRTBres',...
    '%successtreatoldTBcases',...
    '%successtreatnewTBcases',...
    'Diabetesmortalityrate',...
    'Lifeexpectancy',...
    '%childbelow5underweight',...
    'Pop',...
    'Popunder15',...
    'Popover60'});
print('img/PCLoadingmatrix_2012', '-djpeg');

figure;
imagesc(coeff13); colorbar;
xlabel('Principal Component'); 
ylabel('Variable');
title('Variable Loadings for 2013');
set(gca,'Xtick',1:16, 'Ytick',1:16,...
    'YTickLabel',{'%TBwithHIV',...
    '%TBwithHIVtested',...
    '%TBwithHIV-CPT',...
    '%TBwithHIV-ART',...
    '%newRRorMDRTBres',...
    'NoRRorMDRTBres',...
    '%successtreatoldTBcases',...
    '%successtreatnewTBcases',...
    'Hospitaldensity',...
    'Healthpostsdensity',...
    'Heathcentresdensity',...
    'Lifeexpectancy',...
    '%childbelow5underweight',...
    'Pop',...
    'Popunder15',...
    'Popover60'});
print('img/PCLoadingmatrix_2013', '-djpeg');

%% Concatenate results from lasso regression for visualization

header_tmp1 = outerjoin(cell2table(headers_A'), cell2table(headers_B'), 'MergeKeys', true);
header_tmp2 = outerjoin(header_tmp1, cell2table(headers_C'), 'MergeKeys', true);
headers_all = outerjoin(header_tmp2, cell2table(headers_D'), 'MergeKeys', true);

T = table; 
T.varcode = [1:size(headers_all,1)]';
T.varname = table2cell(headers_all);

t2007 = coeff7(:,1:7);
h2007 = headers_A(1:7)';
joined2007 = [headers_A' array2table(t2007)];

t2008 = coeff8(:,[1:2 4:6]);
h2008 = headers_A([1:2 4:6])';
joined2008 = [headers_A' array2table(t2008)];

t2009 = coeff9(:,1:4);
h2009 = headers_A([1:4])';
joined2009 = [headers_A' array2table(t2009)];

t2010 = coeff10(:,[2 5]);
h2010 = headers_B([2 5])';
joined2010 = [headers_B' array2table(t2010)];

t2011 = coeff11(:,1:6);
h2011 = headers_A(1:6)';
joined2011 = [headers_A' array2table(t2011)];

t2012 = coeff12(:,1:4);
h2012 = headers_C(1:4)';
joined2012 = [headers_C' array2table(t2012)];

t2013 = coeff13(:,[1 3 5 9]);
h2013 = headers_D([1 3 5 9])';
joined2013 = [headers_D' array2table(t2013)];

tmp1 = outerjoin(T, joined2007, 'LeftKeys', 2, 'RightKeys', 1, 'MergeKeys', true);
tmp2 = outerjoin(tmp1, joined2008, 'LeftKeys', 2, 'RightKeys', 1, 'MergeKeys', true);
tmp3 = outerjoin(tmp2, joined2009, 'LeftKeys', 2, 'RightKeys', 1, 'MergeKeys', true);
tmp4 = outerjoin(tmp3, joined2010, 'LeftKeys', 2, 'RightKeys', 1, 'MergeKeys', true);
tmp5 = outerjoin(tmp4, joined2011, 'LeftKeys', 2, 'RightKeys', 1, 'MergeKeys', true);
tmp6 = outerjoin(tmp5, joined2012, 'LeftKeys', 2, 'RightKeys', 1, 'MergeKeys', true);
tableall = outerjoin(tmp6, joined2013, 'LeftKeys', 2, 'RightKeys', 1, 'MergeKeys', true);

tt = table2array(tableall(:,3:size(tableall,2)));

%%
figure;
imagesc(tt); colorbar;
%xlabel('Principal Component'); 
%ylabel('Variable');
title('Factor Loadings from Significant PCs across Years');
set(gca,'Xtick',1:32, 'Ytick',1:17,...
    'YTickLabel',{'Diabetes.mortality.rate',...     
    'Healthposts.density',...         
    'Heathcentres.density',...        
    'Hospital.density',...            
    'Life.expectancy',...             
    'No.RRorMDR.TBres',...            
    '%TBwithHIV',...               
    '%TBwithHIV-ART',...           
    '%TBwithHIV-CPT',...           
    '%TBwithHIV.tested',...        
    '%childbelow5.underweight',... 
    '%new.RRorMDR.TBres',...       
    '%success.treat.newTBcases',...
    '%success.treat.oldTBcases',...
    'Pop',...                         
    'Pop.over.60',...                 
    'Pop.under.15'},...
    'XTickLabel',{'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7',...
    'PC1', 'PC2', 'PC4', 'PC5', 'PC6',...
    'PC1', 'PC2', 'PC3', 'PC4',...
    'PC2', 'PC5',...
    'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6',...
    'PC1', 'PC2', 'PC3', 'PC4',...
    'PC1', 'PC3', 'PC5', 'PC9'},...
    'XTickLabelRotation', 70);
print('img/PCLoadingmatrix_LASSO', '-dpdf');
