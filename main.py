from Bio import SeqIO
import random
import csv


# Function to randomly select k-mer motifs from sequences
def select_random_motifs(sequences, k):
    motif_list = []
    for seq in sequences:
        start_position = random.randint(0, len(seq) - k)
        motif = seq[start_position:start_position + k]
        motif_list.append(motif)
    return motif_list

# Function to calculate the profile of the motifs
def create_profile(motifs, k):
    profile = {base: [1] * k for base in "ACGT"}  # initializing with pseudocounts
    for motif in motifs:
        for i in range(k):
            profile[motif[i]][i] += 1
    return profile

# Function to calculate probabilities of k-mers based on the profile
def calculate_kmer_probabilities(seq, k, profile):
    probabilities = []
    for i in range(len(seq) - k + 1):
        prob = 1
        kmer = seq[i:i + k]
        for j in range(k):
            prob *= profile[kmer[j]][j]
        probabilities.append(prob)
    return probabilities

# Gibbs sampling function
def gibbs_sampler(sequences, k, num_iterations):
    motifs = select_random_motifs(sequences, k)
    best_motifs = motifs.copy()

    for _ in range(num_iterations):
        i = random.randint(0, len(sequences) - 1)
        current_motifs = motifs[:i] + motifs[i + 1:]
        profile = create_profile(current_motifs, k)
        probabilities = calculate_kmer_probabilities(sequences[i], k, profile)
        total = sum(probabilities)
        probabilities = [p / total for p in probabilities]
        new_position = random.choices(range(len(sequences[i]) - k + 1), weights=probabilities)[0]
        motifs[i] = sequences[i][new_position:new_position + k]

        # Optionally, check if the new set of motifs is better than the best found so far
        # and update best_motifs if so (based on a chosen scoring function)

    return best_motifs

# Read sequences from FASTA file
sequences = [str(record.seq) for record in SeqIO.parse("reads_motif.fa", "fasta")]

# Set parameters
k = 10  # Length of the motif
num_iterations = 1000  # Number of iterations for Gibbs sampling

# Run Gibbs sampler
motif_list = gibbs_sampler(sequences, k, num_iterations)

# Output motifs
# Save motifs to a CSV file
csv_filename = 'possible_motifs.csv'
with open(csv_filename, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Possible Motif'])  # Writing header
    for motif in motif_list:
        writer.writerow([motif])

print(f"Motifs found and saved in {csv_filename}")

