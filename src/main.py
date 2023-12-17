from Bio import SeqIO
from collections import defaultdict

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



# Main code
# Read sequences from a FASTA file
sequences = [str(record.seq) for record in SeqIO.parse("data/reads_motif.fa", "fasta")]

# Find the most frequent motif in these sequences
motif, count = find_most_frequent_motif(sequences)
print(f"Most Frequent Motif: {motif} (Count: {count})")  # Print the most frequent motif and its count

# Create a probability matrix for this motif
probability_matrix = create_probability_matrix(sequences, motif)
print("Probability Matrix:")
# Print the probability matrix
for base, probabilities in probability_matrix.items():
    print(f"{base}: {['{:.4f}'.format(p) for p in probabilities]}")
