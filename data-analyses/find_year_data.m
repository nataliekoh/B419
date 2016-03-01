function [yoursample] = find_year_data(data, year);
    samples = table2array(data(:,3:size(data,2)));
    samples_year = table2array(data(:,2:size(data,2)));
    [row_tmp c_tmp] = find(samples_year(:,1) == year);
    yoursample = samples(row_tmp, :);
end
