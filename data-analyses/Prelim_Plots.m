%% Plot predictor variables against TBprev for each year

% 2007
fig1 = figure; 
for i = 1:size(xsyear_2007,2),
    subplot(6,4,i);
    scatter(xsyear_2007(:,i), response_2007)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2007: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig1, '../img/2007_prelimplot', '-djpeg');

% 2008
fig2 = figure; 
for i = 1:size(xsyear_2008,2),
    subplot(6,4,i);
    scatter(xsyear_2008(:,i), response_2008)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2008: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig2, '../img/2008_prelimplot', '-djpeg');

% 2009
fig3 = figure; 
for i = 1:size(xsyear_2009,2),
    subplot(6,4,i);
    scatter(xsyear_2009(:,i), response_2009)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2009: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig3, '../img/2009_prelimplot', '-djpeg');

% 2010
fig4 = figure; 
for i = 1:size(xsyear_2010,2),
    subplot(6,4,i);
    scatter(xsyear_2010(:,i), response_2010)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2010: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig4, '../img/2010_prelimplot', '-djpeg');

% 2011
fig5 = figure; 
for i = 1:size(xsyear_2011,2),
    subplot(6,4,i);
    scatter(xsyear_2011(:,i), response_2011)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2011: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig5, '../img/2011_prelimplot', '-djpeg');

% 2012
fig6 = figure; 
for i = 1:size(xsyear_2012,2),
    subplot(6,4,i);
    scatter(xsyear_2012(:,i), response_2012)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2012: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig6, '../img/2012_prelimplot', '-djpeg');

% 2013
fig7 = figure; 
for i = 1:size(xsyear_2013,2),
    subplot(6,4,i);
    scatter(xsyear_2013(:,i), response_2013)
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', '2013: Risk factors (x) by row against TB prevalence (y)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
print(fig7, '../img/2013_prelimplot', '-djpeg');

%% Identify risk factor indexes (non-empty plots) for data analysis
% SEE data-analyses/riskfactorindices.txt 
