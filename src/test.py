from Bio import SeqIO
import random
import math
def parse_fasta(file_path):
    sequences = []  # Initialize an empty list to store sequences
    for record in SeqIO.parse(file_path, "fasta"):  # Parse the fasta file and iterate over each record
        sequences.append(str(record.seq))  # Convert each record's sequence to string and append to the list
    return sequences  # Return the list of sequences

def get_random_motif(sequence, length=10):
    start = random.randint(0, len(sequence) - length)  # Choose a random starting point for the motif
    return sequence[start:start + length]  # Return the substring (motif) of the specified length

def hamming_distance(s1, s2):
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))  # Calculate the Hamming distance between two strings

def score_motif(motif, sequences, max_distance=1):
    score = 0  # Initialize score to 0
    for seq in sequences:  # Iterate over all sequences
        for i in range(len(seq) - len(motif) + 1):  # Iterate over all possible starting points of the motif in the sequence
            if hamming_distance(motif, seq[i:i+len(motif)]) <= max_distance:  # Check if the Hamming distance is within the allowed max
                score += 1  # Increment score if the condition is met
    return score  # Return the total score

def mutate_motif(motif):
    bases = ['A', 'C', 'G', 'T']  # List of nucleotide bases
    mutation_site = random.randint(0, len(motif) - 1)  # Randomly select a mutation site in the motif
    current_base = motif[mutation_site]  # Get the current base at the mutation site
    possible_mutations = [b for b in bases if b != current_base]  # List all possible mutations (excluding the current base)
    new_base = random.choice(possible_mutations)  # Randomly select a new base for mutation
    return motif[:mutation_site] + new_base + motif[mutation_site + 1:]  # Return the motif with the mutated base

def simulated_annealing(sequences, temperature_range, cooling_rate):
    current_temp = temperature_range[0]  # Set the initial temperature
    final_temp = temperature_range[1]  # Set the final temperature

    current_motif = get_random_motif(random.choice(sequences))  # Get a random motif from a random sequence
    current_score = score_motif(current_motif, sequences)  # Score the current motif

    while current_temp > final_temp:  # Continue until the current temperature is reduced to the final temperature
        new_motif = mutate_motif(current_motif)  # Generate a mutated version of the current motif
        new_score = score_motif(new_motif, sequences)  # Score the new motif

        # Decide whether to accept the new motif
        if new_score > current_score or random.random() < math.exp((new_score - current_score) / current_temp):
            current_motif = new_motif  # Accept the new motif
            current_score = new_score  # Update the current score

        current_temp *= cooling_rate  # Reduce the temperature

    return current_motif  # Return the best motif found

def main():
    file_path = "data/reads_motif.fa"  # Set the file path for the fasta file
    sequences = parse_fasta(file_path)  # Parse the fasta file to get sequences

    temperature_range = (100.0, 0.1)  # Define the temperature range for simulated annealing
    cooling_rate = 0.95  # Define the cooling rate

    best_motif = simulated_annealing(sequences, temperature_range, cooling_rate)  # Find the best motif using simulated annealing

    # Write the results to a file
    with open("data/sa_output/motif_output.txt", "w") as file:
        file.write(f"Most likely motif: {best_motif}\n")  # Write the best motif to the file

if __name__ == "__main__":
    main()  # Run the main function if the script is executed as the main program
