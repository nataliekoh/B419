
%% Load data

data = readtable('data/all_condensed_v4.csv');
data = data(:,2:size(data,2));
samples = table2array(data(:,3:size(data,2)));
samples_year = table2array(data(:,2:size(data,2)));
headers = data.Properties.VariableNames;

%% Which years to choose?

no_elements = sum(~isnan(samples),2);
see_years = horzcat(no_elements, data{:,2});

years = [1990:1:2015]'; 
years_sum = zeros(numel(years), 1);

for i = 1:numel(years),
    tmp = zeros(numel(years), 1);
    tmp = see_years(:,2) == years(i);
    [row, col] = find(tmp);
    for j = 1:size(row, 1);
        years_sum(i, 1) = sum(see_years(row(j,1),1));
    end
end

see_years_no = horzcat(years, years_sum);
    % Seems like we have the most data points for years 2007 - 2013.
    
%% Sort data by years to create a 3D matrix

% Note that the matrices do not have the same size
year_1990 = find_year_data(data,1990);
year_1991 = find_year_data(data,1991);
year_1992 = find_year_data(data,1992);
year_1993 = find_year_data(data,1993);
year_1994 = find_year_data(data,1994);
year_1995 = find_year_data(data,1995);
year_1996 = find_year_data(data,1996);
year_1997 = find_year_data(data,1997);
year_1998 = find_year_data(data,1998);
year_1999 = find_year_data(data,1999);
year_2000 = find_year_data(data,2000);
year_2001 = find_year_data(data,2001);
year_2002 = find_year_data(data,2002);
year_2003 = find_year_data(data,2003);
year_2004 = find_year_data(data,2004);
year_2005 = find_year_data(data,2005);
year_2006 = find_year_data(data,2006);
year_2007 = find_year_data(data,2007);
year_2008 = find_year_data(data,2008);
year_2009 = find_year_data(data,2009);
year_2010 = find_year_data(data,2010);
year_2011 = find_year_data(data,2011);
year_2012 = find_year_data(data,2012);
year_2013 = find_year_data(data,2013);
year_2014 = find_year_data(data,2014);
year_2015 = find_year_data(data,2015);

%% Now create the 3D matrix

% Largest year matrix is 195x24.
% We will do PCA on data only from years which meet this matrix size
% because we cannot concatenate matrices of different sizes

% samples_3d = cat(3, year_1990, year_1991, year_1992, year_1993, year_1994,...
%     year_1995, year_1996, year_1997, year_1998, year_1999, year_2000,...
%     year_2001, year_2002, year_2003, year_2004, year_2005, year_2006,...
%     year_2007, year_2008, year_2009, year_2010, year_2011,...
%     year_2012, year_2013, year_2014, year_2015);

samples_3d = cat(3, year_2003, year_2004, year_2006, year_2007, year_2010,...
     year_2011, year_2012, year_2013);

size(samples_3d, 3); % We can do PCA on 8 years' worth of data

%% PCA - ALS Algorithm: 
% 1st 2 PCs explain most of the variance of the data
% Only 1 cluster observable

% Okay... seems like we can only do it by year.
[co_2007,score_2007,latent_2007,~,explained_2007,~] = pca(year_2007, 'Algorithm', 'als');
[co_2008,score_2008,latent_2008,~,explained_2008,~] = pca(year_2008, 'Algorithm', 'als');
[co_2009,score_2009,latent_2009,~,explained_2009,~] = pca(year_2009, 'Algorithm', 'als');
[co_2010,score_2010,latent_2010,~,explained_2010,~] = pca(year_2010, 'Algorithm', 'als');
%t = score*coeff' + repmat(mu,size(samples_3d,1),1);
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

%% PCA - Variational Bayesian PCA [Ilin & Raiko 2010]
% Test first on 2007 data
% 1st 3 PCs explain most variance in data
% 2 clusters, perhaps? 

% Standardize the data 
zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')),...
    std(x, 'omitnan'));
Z2007 = zscore_xnan(year_2007); 

% Now normalize data 
normcol_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, nanmin(x)),...
    bsxfun(@minus, nanmax(x), nanmin(x)));
norm2007 = normcol_xnan(Z2007); 

opts = struct( 'maxiters', 30,...
               'algorithm', 'vb',...
               'uniquesv', 0,...
               'cfstop', [ 100 0 0 ],...
               'rotate2pca', 1,...
               'minangle', 0 );
[ A, S, Mu, V, cv, hp, lc ] = pca_full(norm2007, 2, opts);

Sv = cov_f2d(cv.S, cv.Isv);
% subspace2d(S, Sv);

% hax = tsplot( S, 'g' );
% Sbnd = sqrt(Sv)*3;
% addebars( Sbnd );

Xrec = repmat(Mu, 1, 24) + A*S; % Reconstruction of missing data

[coeff,score,latent,tsquared,explained,mu] = pca(Xrec)

%%
figure(1);
scatter(score(:,1),score(:,2));
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('PCA - VB algorithm on Year 2007'); 

figure(2);
stairs(cumsum(latent)/sum(latent));

figure(3); 
biplot(coeff(:,1:2), 'scores', score(:,1:2), 'VarLabels', headers(1,3:26));

%%

% Hotelling's T2 values represent a measure of the variation in each 
% sample within the model. It indicates how far each sample is from 
% the center (scores = 0) of the model. 
% Low T-Squared statistics indicate a well-fit model.

subspace(coeff, co_2007) 
% Seems like ALS and VBPCA perform similarly? ALS outputs are pretty sketchy though. 

%% PPCA for quick comparison -- errors

% --- incomplete ---

[c2_2007,s2_2007,pcvar_2007,mu2_2007] = ppca(year_2007, 3);

%% PCA - Least Squares PCA [Ilin & Raiko 2010]

% --- may be explored if time allows... ---






