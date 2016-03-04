import csv
import pycountry

# read in csv file
data = []
with open('../data/all_condensed_v6.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)

    for row in reader:
        data.append(row)


countries = {}
for country in pycountry.countries:
    countries[country.name] = country.alpha3


##for row in data:
##    row[1] = countries.get(row[1], row[1])
##
##with open('../data/all_condensed_v7.csv', 'w', newline='') as csvfile:
##    writer = csv.writer(csvfile)
##    writer.writerows(data)
        
print(countries.get('Bolivia'))
