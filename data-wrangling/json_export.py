import csv
import json

# read in csv file
data = []
with open('../data/all_condensed_v7.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)

    for row in reader:
        data.append(row)

data_2010 = {}
for i in range(1, len(data)):
    if data[i][2] == '2010':
        try:
            tb_prev = float(data[i][3]) / (float(data[i][24]) * 1000) * 100
            data_2010[data[i][1]] = tb_prev
        except ValueError as e:
            print(str(e) + ' at row ' + str(i))

with open('tb_prev_2010.json', 'w') as output:
    json.dump(data_2010, output)
