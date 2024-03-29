### Configuration variables ###
########################################################

motif_length = 20                           # length of motif to search for
DNA_seq_path = "data/UCSC_Cat_small.txt"    # path to file containin the sequences to search in
PSSM_zero_correction = 0.1                  # the pseudo-frequency to give unseen bases in the PSSM matrix
ConvergedThreshold = 20                     # how many iterations without score change must pass before concludes that the algorithm has converged

consolePrint = True                         # print results and motif instances
outputFile = "output/genetic-cat_small.csv" # path to output csv to store results of tests, empty string to not save
runTest = 1                                 # how many times to re-run the algorithm

population_size = 128                       # population size for genetic algorithm
tournament_size = 3                         # tournament size for tournament selection
elite_size = 1                              # number of elite organisms to always copy to the next generation
mutation_rate = 0.02                        # mutation rate for genetic algorithm

#######################################################


from src.organism import Organism
from src.evaluator import Evaluator
from src.geneticalgorithm import GeneticAlgorithm
from util import genPSSM, calcScorePSSM

from random import random, randrange, choice
import time
import math
from typing import List
from colorama import init, Fore
init()

DNA_seq_count = -1
DNA_Letters = ["G", "A", "T", "C"]

class motifOrg(Organism):
    def __init__(self, random = True):
        self.DNA: List[str] = []
        self.matches: List[List[int]] = []

        for i in range(DNA_seq_count):
            row = [0, 0]
            self.matches.append(row)

        super().__init__(random)

    def rndInit(self):
        self.DNA = [choice(DNA_Letters) for _ in range(motif_length)]

    def mutate(self, mutationrate: float):
        for i in range(len(self.DNA)):
            if random() < mutationrate:
                self.DNA[i] = choice([x for x in DNA_Letters if x != self.DNA[i]])


    def crossover(self, spouse: 'motifOrg') -> 'motifOrg':
        child = motifOrg(False)
        midpoint = randrange(0, len(self.DNA))
        child.DNA.extend(spouse.DNA[:midpoint])
        child.DNA.extend(self.DNA[midpoint:])
        child.setFitness(-1)
        return child

    def save(self, gen, sequences):
        output = []

        for i, match in enumerate(self.matches):
            output.append(Fore.WHITE + sequences[i])
            row = " " * match[0]
            for j in range(len(self.DNA)):
                l = self.DNA[j]
                color = Fore.RED
                if l == sequences[i][match[0]+j]:
                    color = Fore.GREEN
                row += color + l
            row += Fore.WHITE
            row += " - " + str(match[1])
            output.append(row)

        print(" ")
        for l in output:
            print(l)
        print(" ")

    def getScore(self, sequences):
        instances = []
        for i, m in enumerate(self.matches):
            seq = sequences[i]
            instance = seq[m[0]:m[0] + motif_length]
            instances.append(instance)

        PSSM = genPSSM(instances, motif_length, PSSM_zero_correction)

        total = 0
        for inst in instances:
            total += -math.log(calcScorePSSM(PSSM, inst))

        return total

    def __repr__(self):
        return "".join(self.DNA)

class motifEval(Evaluator):
    def __init__(self, ):

        with open(DNA_seq_path) as file:
            lines = file.readlines()
            self.sequences = [line.rstrip() for line in lines]
            global DNA_seq_count
            DNA_seq_count = len(self.sequences)

        super().__init__()

    def sequenceScore(self, seq: str, motif: motifOrg) -> tuple[int, int]:
        best_match: int = 0
        best_match_i = 0

        for i in range(len(seq) - motif_length):
            window = seq[i: i+motif_length]
            current_match = 0

            for j, letter in enumerate(window):
                current_match += int(letter == motif.DNA[j])

            if current_match > best_match:
                best_match = current_match
                best_match_i = i

        return best_match, best_match_i


    def evalMulti(self, pop: List[motifOrg]) -> float:
        avg_fit = 0

        for i in range(len(pop)):
            if pop[i].getFitness() == -1:
                fit = 0

                for j, seq in enumerate(self.sequences):
                    match, position = self.sequenceScore(seq, pop[i])
                    pop[i].matches[j][1] = match
                    pop[i].matches[j][0] = position

                    fit += float(match) / float(motif_length)

                fit = fit / len(self.sequences)
                pop[i].setFitness(fit)
            else:
                fit = pop[i].getFitness()
            avg_fit += fit

        avg_fit = avg_fit / len(pop)
        return avg_fit


def main():
    eval = motifEval()
    GA = GeneticAlgorithm(motifOrg, eval, populationSize=population_size, eliteSize=elite_size, tournamentSize=tournament_size, cutoffSize=0, mutationRate=mutation_rate)

    test_results = []

    for t in range(runTest):
        print("test " + str(t + 1) + " / " + str(runTest))
        GA.initGenerator()

        conv_counter = 0
        last = 0
        elapconv = 0

        starttime = time.perf_counter()

        while True:
            gen, avg, best_org = GA.nextGeneration()

            elapsed = time.perf_counter() - starttime

            bfit = best_org.getFitness()

            if bfit == last:
                if conv_counter == 0:
                    elapconv = elapsed
                conv_counter += 1
            else:
                conv_counter = 0

            if conv_counter >= ConvergedThreshold:
                total = best_org.getScore(eval.sequences)
                if consolePrint:
                    print("")
                    print("Converged at:")
                    print("{}: best_fit {}, PSSM score: {}, elapsed {:0.4f}".format(gen - ConvergedThreshold + 1, best_org.getFitness(), total/DNA_seq_count, elapconv))
                    best_org.save(gen, eval.sequences)
                test_results.append([total/DNA_seq_count, elapconv])
                break

            last = bfit
            if consolePrint:
                print("{}: best_fit {}, avg_fit {}, elapsed {:0.4f}".format(gen, bfit, avg, elapsed))

    if outputFile != "":
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