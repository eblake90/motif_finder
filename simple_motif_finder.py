# Import the SeqIO module from the Biopython library
from Bio import SeqIO

# Define a function to calculate the Hamming distance between two sequences
def hamming_distance(seq1, seq2):
    # Use a generator expression with sum to count the differing positions
    return sum(count1 != count2 for count1, count2 in zip(seq1, seq2))

# Define a function to find the best k-mer (10bp motif) in a set of sequences
def find_best_kmer(fasta_file):
    # Parse the FASTA file and convert it into a list of sequences
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    # Extract the sequence of the first record as a string
    first_sequence = str(sequences[0].seq)
    # Set the size of the k-mers
    kmer_size = 10

    # Generate all possible k-mers from the first sequence
    kmers = [first_sequence[i:i+kmer_size] for i in range(len(first_sequence) - kmer_size + 1)]

    # Initialize variables to track the best k-mer and its corresponding distance
    best_kmer = None
    lowest_total_distance = float('inf')

    # Iterate through each k-mer in the first sequence
    for kmer in kmers:
        # Initialize a variable to accumulate total distance for this k-mer
        total_distance = 0
        # Compare the k-mer with each subsequent sequence
        for seq in sequences[1:]:
            # Initialize a variable to find the minimum distance in this sequence
            min_distance = float('inf')
            # Slide the k-mer across the sequence to find the best alignment
            for i in range(len(seq) - kmer_size + 1):
                # Calculate the Hamming distance for the current alignment
                distance = hamming_distance(kmer, str(seq.seq[i:i+kmer_size]))
                # Update the minimum distance if a better alignment is found
                if distance < min_distance:
                    min_distance = distance
            # Add the best distance from this sequence to the total distance
            total_distance += min_distance

        # Update the best k-mer if this k-mer has a lower total distance
        if total_distance < lowest_total_distance:
            lowest_total_distance = total_distance
            best_kmer = kmer

    # Return the k-mer with the lowest total Hamming distance
    return best_kmer

# Output
motif = find_best_kmer("path_to_your_fasta_file.fasta")
print("Best motif:", motif)
