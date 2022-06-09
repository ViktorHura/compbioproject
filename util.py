from typing import List
import texttable
import copy


def genPSSM(instances: List[str], motif_length, PSSM_zero_correction):
    matrix = [[0]*motif_length for x in range(4)]

    for i in range(motif_length):
        counts = [0, 0, 0, 0]

        for row in instances:
            counts[0] += row[i] == 'A'
            counts[1] += row[i] == 'C'
            counts[2] += row[i] == 'G'
            counts[3] += row[i] == 'T'

        zero_count = 0
        for j in range(4):
            zero_count += counts[j] == 0

        remainder = 1 - (zero_count * PSSM_zero_correction)

        for j in range(4):
            if counts[j] == 0:
                matrix[j][i] = PSSM_zero_correction
            else:
                matrix[j][i] = counts[j] / len(instances) * remainder

    return matrix


def calcScorePSSM(PSSM, instance):
    score = 1

    for i, letter in enumerate(instance):
        if letter == 'A':
            score *= PSSM[0][i]
        if letter == 'C':
            score *= PSSM[1][i]
        if letter == 'G':
            score *= PSSM[2][i]
        if letter == 'T':
            score *= PSSM[3][i]

    # copyPSSM = copy.deepcopy(PSSM)
    #
    # table = texttable.Texttable()
    # copyPSSM.insert(0, [x for x in range(motif_length)])
    # table.add_rows(copyPSSM)
    # print(table.draw())

    return score