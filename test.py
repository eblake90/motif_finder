from Bio import SeqIO
from collections import defaultdict
import random
from tqdm import tqdm

# Function to calculate the Hamming distance (number of differing characters) between two strings
def hamming_distance(s1, s2):
    if len(s1) != len(s2):  # Check if the strings are of the same length
        raise ValueError("Strings must be of the same length")  # Raise an error if they are not
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))  # Count the number of differences

# Function to find the most frequently occurring motif of a given length in a list of sequences
def find_most_frequent_motif(sequences, motif_length=10):
    motifs = defaultdict(int)  # Create a dictionary to count occurrences of each motif

    # Loop over each sequence
    for seq in sequences:
        # Slide a window of size 'motif_length' over the sequence
        for i in range(len(seq) - motif_length + 1):
            motif = str(seq[i:i + motif_length])  # Extract the motif from the sequence
            motifs[motif] += 1  # Increment the count for this motif

    # Find the motif with the maximum count
    most_frequent_motif = max(motifs, key=motifs.get)
    return most_frequent_motif, motifs[most_frequent_motif]  # Return the most frequent motif and its count

# Function to find the subsequence in a sequence that is most similar to a given motif
def find_closest_subsequence(seq, motif):
    min_distance = len(motif)  # Initialize the minimum distance as the length of the motif
    closest_subseq = None  # Initialize the closest subsequence as None

    # Loop over the sequence to find the closest subsequence
    for i in range(len(seq) - len(motif) + 1):
        subseq = seq[i:i + len(motif)]  # Extract a subsequence
        distance = hamming_distance(subseq, motif)  # Calculate the Hamming distance to the motif
        if distance < min_distance:  # If this is the closest subsequence so far
            min_distance = distance  # Update the minimum distance
            closest_subseq = subseq  # Update the closest subsequence

    return closest_subseq  # Return the closest subsequence

# Function to create a probability matrix for a given motif based on a list of sequences
def create_probability_matrix(sequences, motif):
    matrix = {base: [0] * len(motif) for base in "ATCG"}  # Create a matrix to count nucleotide frequencies

    # Loop over each sequence
    for seq in sequences:
        closest_subseq = find_closest_subsequence(seq, motif)  # Find the subsequence most similar to the motif
        if closest_subseq:
            # Count the frequency of each nucleotide in each position of the closest subsequence
            for i, base in enumerate(closest_subseq):
                matrix[base][i] += 1

    # Convert the counts to probabilities
    for base in matrix:
        for i in range(len(matrix[base])):
            matrix[base][i] /= len(sequences)  # Divide the count by the total number of sequences

    return matrix  # Return the probability matrix

def mutate(kmer):
    position = random.randint(0, len(kmer) - 1)
    nucleotides = ['A', 'T', 'C', 'G']
    nucleotides.remove(kmer[position])
    kmer = list(kmer)
    kmer[position] = random.choice(nucleotides)
    return ''.join(kmer)

def evolutionary_motif_search(sequences, k=10, population_size=100, num_iterations=1000):
    # Initialization
    population = [random.choice(seq[i:i+k]) for seq in sequences for i in range(len(seq) - k + 1)]
    population = random.sample(population, population_size)

    best_kmer = None
    best_fitness = float('-inf')

    for iteration in tqdm(range(num_iterations), desc="Searching Motifs"):
        # Fitness Calculation
        fitness_scores = [calculate_fitness(kmer, sequences) for kmer in population]

        # Check for the best k-mer in this iteration
        for i, kmer in enumerate(population):
            if fitness_scores[i] > best_fitness:
                best_fitness = fitness_scores[i]
                best_kmer = kmer

        # Selection
        selected_kmers = select_top_kmers(population, fitness_scores, 10)

        # Mutation
        mutated_kmers = [mutate(kmer) for kmer in selected_kmers]

        # Update Population
        population = mutated_kmers + random.sample(population, population_size - len(mutated_kmers))

    # Return the best k-mer
    return best_kmer

# Implement the fitness function and selection function
def calculate_fitness(kmer, sequences):
    total_distance = 0
    for seq in sequences:
        min_distance = min(hamming_distance(kmer, seq[i:i+len(kmer)]) for i in range(len(seq) - len(kmer) + 1))
        total_distance += min_distance
    average_distance = total_distance / len(sequences)
    fitness = 1 / (1 + average_distance)  # Higher fitness for lower average distance
    return fitness

def select_top_kmers(population, fitness_scores, number):
    # Pair each k-mer with its fitness score
    paired_kmers = list(zip(population, fitness_scores))
    # Sort the k-mers based on fitness scores
    sorted_kmers = sorted(paired_kmers, key=lambda x: x[1], reverse=True)
    # Select the top 'number' k-mers
    top_kmers = [kmer for kmer, score in sorted_kmers[:number]]
    return top_kmers





# Main code
# Read sequences from a FASTA file
sequences = [str(record.seq) for record in SeqIO.parse("reads_motif.fa", "fasta")]

# Find the most frequent motif in these sequences
motif, count = find_most_frequent_motif(sequences)
print(f"Most Frequent Motif: {motif} (Count: {count})")  # Print the most frequent motif and its count

# Create a probability matrix for this motif
probability_matrix = create_probability_matrix(sequences, motif)
print("Probability Matrix:")
# Print the probability matrix
for base, probabilities in probability_matrix.items():
    print(f"{base}: {['{:.4f}'.format(p) for p in probabilities]}")


# Use the function on your sequences
best_motif = evolutionary_motif_search(sequences)
print(f"Identified Motif: {best_motif}")