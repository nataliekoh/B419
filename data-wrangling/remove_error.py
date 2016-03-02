import csv

# read in csv file
data = []
with open('../data/all_condensed_v5.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)

    for row in reader:
        data.append(row)

percents = [7, 8]
for row in range(1, len(data)):
    for col in percents:
        value = data[row][col - 1]
        if value:
            if float(value) > 100:
                print(row, col, value)
                data[row][col - 1] = ''
        
with open('../data/all_condensed_v6.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)
        

