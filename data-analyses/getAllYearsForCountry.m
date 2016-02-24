function [filtered] = getAllYearsForCountry(file, country)
%getAllYearsForCountry Takes in original dataset in table form and a specific
% country (as string). Returns a table with all data for that country.

country_indices = find(strcmp(file{:, 1}, country) == 1);
filtered = file(country_indices, :);
end

