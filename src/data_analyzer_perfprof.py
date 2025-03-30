import csv
import os
from collections import defaultdict

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def process_csv(input_file, output_file):
    # Dizionario per memorizzare i dati
    data = defaultdict(dict)
    algorithms = set()
    instances = set()

    # Leggi il file CSV di input
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile)
        for row in reader:
            algorithm, instance, score = row[0], row[1], float(row[2])
            algorithms.add(algorithm)
            instances.add(instance)
            data[algorithm][instance] = score

    # Ordina algoritmi e istanze
    algorithms = sorted(algorithms)
    instances = sorted(instances)

    # Scrivi il file CSV di output
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)

        # Prima riga: numero di algoritmi e nomi degli algoritmi
        header = [len(algorithms)] + algorithms
        writer.writerow(header)

        # Scrivi i dati
        for instance in instances:
            row = [instance]
            for algorithm in algorithms:
                row.append(data[algorithm].get(instance, ''))
            writer.writerow(row)

    print("Output file path:", os.path.abspath(output_file))
    print(f"File di output scritto correttamente: {os.path.abspath(output_file)}")

# Esempio di utilizzo
input_file = '../data/results.csv'  # Sostituisci con il percorso del file di input
output_file = '../data/output.csv'  # Sostituisci con il percorso del file di output
process_csv(input_file, output_file)