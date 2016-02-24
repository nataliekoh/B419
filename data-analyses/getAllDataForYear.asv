function [filtered] = getAllDataForYear(file, year)
%GETALLDATAFORYEAR Takes in original dataset in table form and a specific
% year. Returns a table with data for every country in that year.

year_indices = find(file{:, 2} == year);
filtered = file(year_indices, :);

end

