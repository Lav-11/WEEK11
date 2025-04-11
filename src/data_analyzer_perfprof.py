import csv
import os
from collections import defaultdict

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# python3 ../src/perfprof.py -M 1.05 ../data/output.csv ../data/output.pdf

def process_csv(input_file, output_file):
    data = defaultdict(dict)
    algorithms = set()
    instances = set()

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile)
        for row in reader:
            algorithm, instance, score = row[0], row[1], float(row[2])
            algorithms.add(algorithm)
            instances.add(instance)
            data[algorithm][instance] = score

    algorithms = sorted(algorithms)
    instances = sorted(instances)

    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        header = [len(algorithms)] + algorithms
        writer.writerow(header)

        for instance in instances:
            row = [instance]
            for algorithm in algorithms:
                row.append(data[algorithm].get(instance, ''))
            writer.writerow(row)

input_file = '../data/results.csv'
output_file = '../data/output.csv'
process_csv(input_file, output_file)
