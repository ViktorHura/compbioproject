import math
from random import randrange
import time
from typing import List
from colorama import init, Fore
init()

from util import genPSSM, calcScorePSSM

DNA_Letters = ["A", "C", "G", "T"]
motif_length = 10
DNA_seq_path = "data/motif_seq_small.txt"
DNA_seq_count = -1
PSSM_zero_correction = 0.1
ConvergedThreshold = 100

consolePrint = False
outputFile = "output/gibbs-small.csv"
runTest = 50

def load_sequences() -> List[str]:
    with open(DNA_seq_path) as file:
        lines = file.readlines()
        sequences = [line.rstrip() for line in lines]
        global DNA_seq_count
        DNA_seq_count = len(sequences)

    return sequences


def choose_new_instance(PSSM, sequence) -> int:
    max_pos = len(sequence) - motif_length
    instance_weights = [0] * (max_pos)

    for i in range(max_pos):
        substring = sequence[i:i+motif_length]
        instance_weights[i] = calcScorePSSM(PSSM, substring)

    best = max(instance_weights)
    return instance_weights.index(best)


def gibbs_iterate(s: List[int], sequences: List[str]):
    chosen_sequence = randrange(0, len(s))

    kmin1_instances = []

    for i, pos in enumerate(s):
        if i is chosen_sequence:
            continue

        seq = sequences[i]
        instance = seq[pos:pos+motif_length]
        kmin1_instances.append(instance)

    PSSM = genPSSM(kmin1_instances, motif_length, PSSM_zero_correction)

    s[chosen_sequence] = choose_new_instance(PSSM, sequences[chosen_sequence])


def printInstances(s: List[int], sequences: List[str]):
    output = []

    for i, pos in enumerate(s):
        seq = sequences[i]

        row = Fore.WHITE
        row += seq[:pos]
        row += Fore.GREEN
        row += seq[pos:pos+motif_length]
        row += Fore.WHITE
        row += seq[pos+motif_length:]
        output.append(row)

    print(" ")
    for l in output:
        print(l)
    print(" ")


def main():
    sequences = load_sequences()

    test_results = []

    for t in range(runTest):
        print("test " + str(t+1) + " / " + str(runTest))

        s: List[int] = []
        for i in sequences:
            s.append(randrange(0, len(i) - motif_length))

        if consolePrint:
            printInstances(s, sequences)

        gen = 0

        conv_counter = 0
        last = 0
        elapconv = 0

        starttime = time.perf_counter()

        while True:
            gibbs_iterate(s, sequences)

            k_instances = []

            for i, pos in enumerate(s):
                seq = sequences[i]
                instance = seq[pos:pos + motif_length]
                k_instances.append(instance)

            PSSM = genPSSM(k_instances, motif_length, PSSM_zero_correction)


            total = 0
            for inst in k_instances:
                total += -math.log(calcScorePSSM(PSSM, inst))
            elapsed = time.perf_counter() - starttime

            if total == last:
                if conv_counter == 0:
                    elapconv = elapsed
                conv_counter += 1
            else:
                conv_counter = 0

            if conv_counter >= ConvergedThreshold:
                if consolePrint:
                    print("{}: score {}, elapsed {:0.4f}".format(gen - ConvergedThreshold + 1, total, elapconv))
                    printInstances(s, sequences)
                test_results.append([total, elapconv])
                break

            last = total

            if consolePrint:
                print("{}: score {}, elapsed {:0.4f}".format(gen, total, elapsed))
            gen += 1

    import csv
    header = ["score", "time"]

    with open(outputFile, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for r in test_results:
            writer.writerow(r)
    file.close()


if __name__ == '__main__':
    main()