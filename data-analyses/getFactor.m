function [filtered] = getFactor(file, col, year)
%getFactor takes in a specific factor (by column index) and a specific year
%and returns a matrix consisting of values for every country.

year_indices = find(file{:, 3} == year);
filtered = file{year_indices, col};
end

