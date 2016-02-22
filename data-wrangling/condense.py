import csv

# read in csv file
data = []
with open('../data/merged_all_datasets.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)

    for row in reader:
        data.append(row)
   
# trim extra whitespace from countries
for row in data:
    row[1] = row[1].strip()
    row[2] = row[2].strip()
    row[3] = row[3].replace(' ', '')

# data_b will need to be condensed into data_a
data_a = data[0:4929]
data_b = data[4929:4957]


# need to merge column 21, 14
def condense_life_expect():
    # merging column 21
    processed = sorted(data_a, key=lambda c: c[1] + c[2])
    i_remove = []
    
    for i in range(1, len(processed)):
        if processed[i - 1][2] == processed[i][2]:
            processed[i][21] = processed[i-1][21]
            i_remove.append(i-1)
    
    count = 0
    for i in i_remove:
        processed.pop(i - count)
        count += 1
    
    # merging column 14
    processed += data_b
    processed = sorted(processed, key=lambda c: c[1] + c[2])
    i_remove = []
    
    for i in range(1, len(processed)):
        if processed[i - 1][2] == processed[i][2]:
            processed[i-1][14] = processed[i][14]
            i_remove.append(i)
    
    count = 0
    for i in i_remove:
        processed.pop(i - count)
        count += 1
    
    with open('../data/all_condensed.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(processed)
        
    return processed

test = condense_life_expect()
