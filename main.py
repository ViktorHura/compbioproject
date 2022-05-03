import datetime
import os

from src.organism import Organism
from src.evaluator import Evaluator
from src.geneticalgorithm import GeneticAlgorithm

from random import random, randrange, choice
import string
import math
import time
from typing import List


DNA_Letters = ["G", "A", "T", "C"]
motif_length = 10
DNA_seq_path = "data/motif_seq_small.txt"
DNA_seq_count = -1

outputDir = "output"

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
            output.append(sequences[i])
            row = " " * match[0]
            for l in self.DNA:
                row += l

            row += " - " + str(match[1])
            output.append(row)

        textfile = open(outputDir + "/" + str(gen) + "_" + str(self.fitness) + ".txt", "w")
        for line in output:
            textfile.write(line + "\n")
        textfile.close()

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
    GA = GeneticAlgorithm(motifOrg, eval, populationSize=64, eliteSize=1, cutoffSize=12)
    GA.initGenerator()

    now = datetime.datetime.now().strftime("%d-%m-%H-%M_%S")
    global outputDir
    outputDir = os.path.join(os.getcwd(), outputDir, "out-" + now)
    os.mkdir(outputDir)

    starttime = time.perf_counter()
    bfit = 0
    best_org = None

    while bfit != 1:
        gen, avg, best_org = GA.nextGeneration()

        if gen % 1 == 0:
            elapsed = time.perf_counter() - starttime
            print("{}: best_fit {}, avg_fit {}, elapsed {:0.4f}".format(gen, best_org.getFitness(), avg, elapsed))
            print(best_org)
            if best_org.getFitness() > bfit:
                best_org.save(gen, eval.sequences)
                bfit = best_org.getFitness()

    print(best_org, GA.getGen())


if __name__ == '__main__':
    main()