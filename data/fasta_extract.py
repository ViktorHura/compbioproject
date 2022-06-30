from Bio import SeqIO

seq_count = 50
seq_length = 5000

def main():
    fasta_sequences = SeqIO.parse(open("felCat8.fa"), 'fasta')
    i = 0

    with open('UCSC_Cat.txt', 'w') as f:
        for fasta in fasta_sequences:
            if i == seq_count:
                break
            f.write(str(fasta.seq).replace("N", "")[:seq_length].upper())
            f.write('\n')
            i += 1

    f.close()

if __name__ == '__main__':
    main()