%% Compare different methods for dimensionality reduction
% Encompasses: ALS-PCA, PPCA, VBPCA

%% Load data

data = readtable('../data/all_condensed_v6.csv');
data = data(:,2:size(data,2));
samples = table2array(data(:,3:size(data,2)));
samples_year = table2array(data(:,2:size(data,2)));
headers = data.Properties.VariableNames;

% Here we extract numerical arrays using find_year_data
year_2007 = find_year_data(data,2007);
year_2008 = find_year_data(data,2008);
year_2009 = find_year_data(data,2009);
year_2010 = find_year_data(data,2010);
year_2011 = find_year_data(data,2011);
year_2012 = find_year_data(data,2012);
year_2013 = find_year_data(data,2013);

% Here we extract data as tables using getAllDataForYear
data_2007 = getAllDataForYear(data,2007);
data_2008 = getAllDataForYear(data,2008);
data_2009 = getAllDataForYear(data,2009);
data_2010 = getAllDataForYear(data,2010);
data_2011 = getAllDataForYear(data,2011);
data_2012 = getAllDataForYear(data,2012);
data_2013 = getAllDataForYear(data,2013);

%% ALS-PCA: 2007-2011

[co_2007,score_2007,latent_2007,~,explained_2007,~] = pca(year_2007, 'Algorithm', 'als');
[co_2008,score_2008,latent_2008,~,explained_2008,~] = pca(year_2008, 'Algorithm', 'als');
[co_2009,score_2009,latent_2009,~,explained_2009,~] = pca(year_2009, 'Algorithm', 'als');
[co_2010,score_2010,latent_2010,~,explained_2010,~] = pca(year_2010, 'Algorithm', 'als');
[co_2011,score_2011,latent_2011,~,explained_2011,~] = pca(year_2011, 'Algorithm', 'als');

figure(1); 
scatter(score_2007(:,1), score_2007(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - ALS algorithm on Year 2007'); 

figure(2); 
scatter(score_2008(:,1), score_2008(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - ALS algorithm on Year 2008'); 

figure(3); 
scatter(score_2009(:,1), score_2009(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - ALS algorithm on Year 2009'); 

figure(4); 
scatter(score_2010(:,1), score_2010(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - ALS algorithm on Year 2010'); 

figure(5);  
stairs(cumsum(latent_2007)/sum(latent_2007));

figure(6);  
stairs(cumsum(latent_2008)/sum(latent_2008));

figure(7);  
stairs(cumsum(latent_2009)/sum(latent_2009));

figure(8);  
stairs(cumsum(latent_2010)/sum(latent_2010));

figure(9);
stairs(cumsum(latent_2011)/sum(latent_2011));

%% PPCA: This will not run without subsetting data further
% only 1 component allowed because of rank deficient matrix
% cannot handle very sparse data

[coeff_2007, score_2007, pcvar_2007, mu_2007] = ppca(year_2007, 1);
[coeff_2008, score_2008, pcvar_2008, mu_2008] = ppca(year_2008, 1);
[coeff_2010, score_2010, pcvar_2010, mu_2010] = ppca(year_2010, 1);

%% VBPCA: 2007-2012

% Standardize the data 
zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')),...
    std(x, 'omitnan'));
Z2007 = zscore_xnan(year_2007); 
Z2008 = zscore_xnan(year_2008); 
Z2009 = zscore_xnan(year_2009); 
Z2010 = zscore_xnan(year_2010); 
Z2011 = zscore_xnan(year_2011); 
Z2012 = zscore_xnan(year_2012); 

% Now normalize data 
normcol_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, nanmin(x)),...
    bsxfun(@minus, nanmax(x), nanmin(x)));
norm2007 = normcol_xnan(Z2007); 
norm2008 = normcol_xnan(Z2008); 
norm2009 = normcol_xnan(Z2009); 
norm2010 = normcol_xnan(Z2010); 
norm2011 = normcol_xnan(Z2011);
norm2012 = normcol_xnan(Z2012); 

% Specify the options for VB imputation
opts = struct( 'maxiters', 50,...
               'algorithm', 'vb',...
               'uniquesv', 0,...
               'cfstop', [ 100 0 0 ],...
               'rotate2pca', 1,...
               'minangle', 0 );
[A07, S07, Mu07, V07, cv07, hp07, lc07] = pca_full(norm2007, 2, opts);
[A08, S08, Mu08, V08, cv08, hp08, lc08] = pca_full(norm2008, 2, opts);
[A09, S09, Mu09, V09, cv09, hp09, lc09] = pca_full(norm2009, 2, opts);
[A10, S10, Mu10, V10, cv10, hp10, lc10] = pca_full(norm2010, 2, opts);
[A11, S11, Mu11, V11, cv11, hp11, lc11] = pca_full(norm2011, 2, opts);
[A12, S12, Mu12, V12, cv12, hp12, lc12] = pca_full(norm2012, 2, opts);

% Reconstruction of missing data
Xrec2007 = repmat(Mu07, 1, 24) + A07*S07; 
Xrec2008 = repmat(Mu08, 1, 24) + A08*S08; 
Xrec2009 = repmat(Mu09, 1, 24) + A09*S09; 
Xrec2010 = repmat(Mu10, 1, 24) + A10*S10; 
Xrec2011 = repmat(Mu11, 1, 24) + A11*S11; 
Xrec2012 = repmat(Mu12, 1, 24) + A12*S12; 

% Feed in the new model into MATLAB's PCA function
[coeff07, score07, lat07, tsq07, ex07, mu07] = pca(Xrec2007);
[coeff08, score08, lat08, tsq08, ex08, mu08] = pca(Xrec2008);
[coeff09, score09, lat09, tsq09, ex09, mu09] = pca(Xrec2009);
[coeff10, score10, lat10, tsq10, ex10, mu10] = pca(Xrec2010);
[coeff11, score11, lat11, tsq11, ex11, mu11] = pca(Xrec2011);
[coeff12, score12, lat12, tsq12, ex12, mu12] = pca(Xrec2012);

% Plot some figures
figure(1); 
scatter(score07(:,1), score07(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2007'); 

figure(2); 
scatter(score08(:,1), score08(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2008'); 

figure(3); 
scatter(score09(:,1), score09(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2009'); 

figure(4); 
scatter(score10(:,1), score10(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2010'); 

figure(5); 
scatter(score11(:,1), score11(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2011'); 

figure(6); 
scatter(score12(:,1), score12(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2012'); 

figure(7);  
stairs(cumsum(lat07)/sum(lat07));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2007');

figure(8);  
stairs(cumsum(lat08)/sum(lat08));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2008');

figure(9);  
stairs(cumsum(lat09)/sum(lat09));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2009');

figure(10);  
stairs(cumsum(lat10)/sum(lat10));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2010');

figure(11);
stairs(cumsum(lat11)/sum(lat11));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2011');

figure(12);
stairs(cumsum(lat12)/sum(lat12));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2012');


%% Figures for interim report from VBPCA results

close all;

fig1 = figure; 
biplot(coeff10(:,1:2), 'scores', score10(:,1:2), 'VarLabels', headers(1,3:26));
title('Biplot of 2010 data using VBPCA')
    % Each point represents a country;
    % Each vector line represents a variable;
    % Length of each vector represents variance; longer = higher variance;
    % Angle between vectors represents degree of correlation;
%print(fig1, 'img/2010_Biplot', '-depsc2');

fig2 = figure; 
biplot(coeff07(:,1:2), 'scores', score07(:,1:2), 'VarLabels', headers(1,3:26));
title('Biplot of 2007 data using VBPCA')
%print(fig2, 'img/2007_Biplot', '-depsc2');

fig3 = figure; 
subplot(1,2,1);
scatter(score07(:,1), score07(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2007'); 
subplot(1,2,2);
scatter(score10(:,1), score10(:,2), 15, 'b', 'filled');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2010');
%print(fig3, 'img/PCplot_0710', '-depsc2');

fig4 = figure;
subplot(1,2,1);
stairs(cumsum(lat07)/sum(lat07));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2007');
refline(0,0.95)
subplot(1,2,2);
stairs(cumsum(lat10)/sum(lat10));
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2010');
refline(0,0.95)
%print(fig4, 'img/PCvarplot_0710', '-dpsc');

X = 1:1:24; 
Y1 = cumsum(lat07)/sum(lat07);
Y2 = cumsum(lat10)/sum(lat10);
[xb,yb] = stairs(X,Y1);
[xb2,yb2] = stairs(X,Y2);

fig5 = figure;
subplot(1,2,1);
plot(xb, yb);
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2007');
ax = gca;
ax.XTick = [1:2:24];
subplot(1,2,2);
plot(xb2, yb2);
xlabel('No. of Principal Components');
ylabel('% Variance Explained');
title('PCA - VB algorithm on Year 2010');
ax = gca;
ax.XTick = [1:2:24];
%print(fig5, 'img/xPCvarplot_0710', '-depsc');

%% Save these figures

% print(fig1, 'img/2010_Biplot', '-depsc2');
% print(fig2, 'img/2007_Biplot', '-depsc2');
% print(fig3, 'img/PCplot_0710', '-depsc2');
% print(fig4, 'img/PCvarplot_0710', '-depsc2');
% print(fig5, 'img/xPCvarplot_0710', '-depsc2');











