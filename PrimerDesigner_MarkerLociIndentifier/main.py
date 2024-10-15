from abc import ABC, abstractmethod
import argparse
import primer3
import configparser
from Bio import SeqIO
from Bio.Seq import Seq
from typing import Dict, Tuple, List, Any
import heapq
import subprocess
import csv
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from torch.nn.utils.rnn import pad_sequence
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
import ray
from ray import tune
from ray.tune import Trainable
import os
from gensim.models import Word2Vec
from hyperopt import hp
from ray.tune.suggest.hyperopt import HyperOptSearch
import glob
import psycopg2
from psycopg2 import sql
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline, MafftCommandline, MuscleCommandline, PrankCommandline
from Bio import SeqIO
from Bio.Align import AlignInfo
from tempfile import NamedTemporaryFile
import re
from Bio import Entrez
import zipfile
import requests
import json
from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.api.assembly_metadata_api import AssemblyMetadataApi
from ncbi.datasets.openapi.api.genome_api import GenomeApi
from collections import defaultdict

# Parsing the input FASTA file to a dictionary
def parse_fasta_to_dict(fasta_file):
    sequence_dict = {}
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence_dict[record.id] = str(record.seq)
    return sequence_dict


def parse_fasta_by_amplicons(fasta_file):
    amplicon_dict = {}

    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            # Extract barcode from record.id

            identifier = record.id.split('_')[0]

            # Initialize a sub-dictionary for this identifier if it doesn't exist
            if identifier not in amplicon_dict:
                amplicon_dict[identifier] = {}

            # Add the sequence to the sub-dictionary
            amplicon_dict[identifier][record.id] = str(record.seq)

    return amplicon_dict


def extract_info_multiFASTA(header):
    """
    Extract the genome name and ID in the headers of FASTA files in the format of BV-BRC
    """

    # Definition of a regular expression to match the genome name and ID in FASTA files from BV-BRC
    pattern = r'\[([^\|]+) \| ([^\]]+)\]'

    # Match the header with the defined regular expression
    match = re.search(pattern, header)

    # Check if a match is found
    if match:
        # Extract and return Genome_Name and Genome_ID
        genome_name = match.group(1)
        genome_id = match.group(2)
        return genome_name, genome_id
    else:
        return None, None


def read_multifasta(file_path):
    """
    Read in Multi-FASTA files and write Genome Name, ID and sequence to a pandas dataframe
    """

    # Initialize an empty DataFrame with column names
    # columns = ['Genome Name', 'Genome ID', 'Sequence']
    gene_df = pd.DataFrame(columns=['genome_name', 'genome_id', 'sequence'])

    with open(file_path, 'r') as file:
        current_header = None
        current_sequence = ""

        for line in file:

            if line.startswith('>'):
                # The lines starting with > are the headers

                # If data has already been read in, we need to safe it to the data frame before defining the newly read header as "current_header"
                if current_header is not None:
                    genome_name, genome_id = extract_info_multiFASTA(current_header)
                    new_entry = {'genome_name': genome_name, 'genome_id': genome_id, 'sequence': current_sequence}
                    gene_df = gene_df.append(new_entry, ignore_index=True)
                    # reset the current_sequence
                    current_sequence = ""

                # Set the current_header to the newly read header line
                current_header = line[1:]
            else:
                # Add the sequence line to the sequence
                current_sequence += line.rstrip('\n')

        # Process the last sequence after the loop ends
        if current_header is not None:
            genome_name, genome_id = extract_info_multiFASTA(current_header)
            new_entry = {'genome_name': genome_name, 'genome_id': genome_id, 'sequence': current_sequence}
            gene_df = gene_df.append(new_entry, ignore_index=True)

    return gene_df


# Convert parameters to the correct type
def convert_param(value):
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


# Function to run BLAST for a single primer sequence
def blast_primer(primer_sequence, database):

    result = subprocess.run(
        ['blastn', '-query', '-', '-db', database, '-outfmt',
         '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'],
        input=primer_sequence, encoding='ascii', stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if result.returncode != 0:
        # If BLAST didn't execute successfully, print the error and return an empty result
        print(f"Error running BLAST: {result.stderr}")
        return []

    # Parse the CSV output from BLAST and return it
    return result.stdout.strip().split('\n')


def parse_blast_output(blast_output):
    reader = csv.DictReader(blast_output, fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    return [row for row in reader]


def is_specific(primer_blast_results, params):

    for result in primer_blast_results:
        if float(result['pident']) >= params['IDENTITY_THRESHOLD'] and float(result['evalue']) <= params['EVALUE_CUTOFF']:
            return True  # Primer is specific to at least one sequence
    return False  # No specific alignment found


def output_primers_to_csv(left_primers_info, right_primers_info, left_file_name, right_file_name):
    """
    Output left and right primer information to separate CSV files.

    Parameters:
    left_primers_info (dict): Dictionary containing left primer information.
    right_primers_info (dict): Dictionary containing right primer information.
    left_file_name (str): Name of the output CSV file for left primers.
    right_file_name (str): Name of the output CSV file for right primers.
    """

    # Function to write a single primer info to a CSV file
    # todo: possibly restructure as a standalone function because of the redefinition-related overhead, but most likely it not critical in this case
    def write_primers_to_csv(primers_info, file_name):
        with open(file_name, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Primer Sequence', 'Tm', 'GC Content', 'Dimer Partner', 'Dimer Potential',
                                '3-Prime End Complementarity'])

            for primer_seq, dimer_info in primers_info.items():
                if dimer_info:
                    first_partner = dimer_info[0]
                    csvwriter.writerow([primer_seq, first_partner['partner_tm'], first_partner['partner_gc'],
                                        first_partner['partner_sequence'], first_partner['dimer_potential'],
                                        first_partner['three_prime_complementarity']])

                for partner in dimer_info[1:]:
                    csvwriter.writerow(['', '', '', partner['partner_sequence'], partner['dimer_potential'],
                                        partner['three_prime_complementarity']])

    # Write left and right primers to separate CSV files
    write_primers_to_csv(left_primers_info, left_file_name)
    write_primers_to_csv(right_primers_info, right_file_name)


class PrimerSummarizerAlgorithm(ABC):
    @abstractmethod
    def summarize(self, primers: Dict[str, Dict[object]], params) -> Tuple[List[Tuple[Any, Any]], List[Tuple[Any, Any]]]:
        pass


# Implement the Frequency-based Strategy
class FrequencyBasedSummarizer(PrimerSummarizerAlgorithm):
    def summarize(self, primers: Dict[str, Dict[object]], params) -> Tuple[List[Tuple[Any, Any]], List[Tuple[Any, Any]]]:
        # Initialize dictionaries to count the frequency of each primer
        left_primer_freq = {}
        right_primer_freq = {}

        # Iterate over all sets of primer pairs for each sequence
        for primer_set in primers.values():
            for primer_pair in primer_set:
                # Count the left primer frequency
                left_primer = primer_pair['left']
                left_primer_tm = primer_pair['left_tm']
                left_primer_gc = primer_pair['left_gc']

                if left_primer in left_primer_freq:
                    left_primer_freq[left_primer]['freq'] += 1
                else:
                    left_primer_freq[left_primer]['freq'] = 1
                    left_primer_freq[left_primer]['tm'] = left_primer_tm
                    left_primer_freq[left_primer]['gc'] = left_primer_gc

                # Count the right primer frequency
                right_primer = primer_pair['right']
                right_primer_tm = primer_pair['right_tm']
                right_primer_gc = primer_pair['right_gc']

                if right_primer in right_primer_freq:
                    right_primer_freq[right_primer]['freq'] += 1
                else:
                    right_primer_freq[right_primer]['freq'] = 1
                    right_primer_freq[right_primer]['tm'] =  right_primer_tm
                    right_primer_freq[right_primer]['gc'] = right_primer_gc

        # Find the most frequent left and right primers
        most_frequent_left_primers = heapq.nlargest(params['n_most_frequent'], left_primer_freq, key=lambda k: left_primer_freq[k]['freq'])
        most_frequent_right_primers = heapq.nlargest(params['n_most_frequent'], right_primer_freq, key=lambda k: right_primer_freq[k]['freq'])

        # Return the most frequent primer pair with corresponding values of Tm and GC content
        most_frequent_left_primers_info = [(primer, left_primer_freq[primer]) for primer in most_frequent_left_primers]
        most_frequent_right_primers_info = [(primer, right_primer_freq[primer]) for primer in most_frequent_right_primers]

        return most_frequent_left_primers_info, most_frequent_right_primers_info


class ConsensusBasedSummarizer(PrimerSummarizerAlgorithm):
    def summarize(self, primers: Dict[str, Dict[object]], params) -> Dict[str, object]:
        pass


# Define the specificity check algorithm interface
class SpecificityCheckAlgorithm(ABC):
    @abstractmethod
    def specificity_check(self, primers, database, params):
        pass


class SpecificityCheckBLAST(SpecificityCheckAlgorithm):
    def specificity_check(self, primers, database, params):
        # Check specificity for left and right primers
        left_primers_specific = {}
        right_primers_specific = {}

        for left_primer, properties in primers[0]:
            if is_specific(parse_blast_output(blast_primer(left_primer, database)), params):
                left_primers_specific[left_primer] = properties

        for right_primer, properties in primers[0]:
            if is_specific(parse_blast_output(blast_primer(right_primer, database)), params):
                right_primers_specific[right_primer] = properties

        return left_primers_specific, right_primers_specific

class PrimerDimersCheck:
    def is_complementary(self, base1, base2):
        """Check if two bases are complementary."""
        base_pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return base_pairs.get(base1, '') == base2

    def check_3prime_complementarity(self, primer1, primer2, length=4):
        """
        Checks if the last 'length' bases of the 3'-end of two primers are complementary.

        Parameters:
        primer1 (str): Sequence of the first primer.
        primer2 (str): Sequence of the second primer (will be reversed).
        length (int): Number of bases at the 3'-end to check for complementarity.

        Returns:
        bool: True if complementary, False otherwise.
        """

        primer2 = primer2[::-1]

        # Create a translation table for complementarity
        complement = str.maketrans('ATCG', 'TAGC')

        # Check the last 'length' bases for complementarity
        return primer1[-length:].translate(complement) == primer2[:length]

    def sliding_window_dimer_check(self, primer1, primer2, window_size=3):
        """Check for dimer formation using a sliding window approach."""
        max_complementarity = 0
        for i in range(len(primer1) - window_size + 1):
            complementarity = sum(self.is_complementary(base1, base2)
                                  for base1, base2 in zip(primer1[i:i+window_size], primer2[:window_size]))
            max_complementarity = max(max_complementarity, complementarity)
        return max_complementarity

    def check_all_primer_combinations(self, primers_info):
        """
        Check all combinations of primers for potential dimer formation, exclude primers with high dimer potential,
        and return two dictionaries, one for left primers and one for right ones, with remaining primer sequences,
        their 'tm' and 'gc' info and their potential dimer partners with their 'tm' and 'gc' info.
        """
        # Initial primer dimer info
        left_primer_info = {}
        right_primer_info = {}

        # Process left primers against all other primers
        for i, (left_primer_seq, left_primer_data) in enumerate(primers_info[0]):
            for j, (other_primer_seq, other_primer_data) in enumerate(primers_info[0] + primers_info[1]):
                if i != j and self.check_for_dimers(left_primer_seq, other_primer_seq):
                    left_primer_info.setdefault(left_primer_seq, []).append({
                        'partner_sequence': other_primer_seq,
                        'partner_tm': other_primer_data['tm'],
                        'partner_gc': other_primer_data['gc']
                    })

        # Process right primers against all other primers
        for i, (right_primer_seq, right_primer_data) in enumerate(primers_info[1]):
            for j, (other_primer_seq, other_primer_data) in enumerate(primers_info[1] + primers_info[0]):
                if i != j and self.check_for_dimers(right_primer_seq, other_primer_seq):
                    right_primer_info.setdefault(right_primer_seq, []).append({
                        'partner_sequence': other_primer_seq,
                        'partner_tm': other_primer_data['tm'],
                        'partner_gc': other_primer_data['gc']
                    })

        # Identify primers to exclude
        left_to_exclude = self.identify_primers_to_exclude(left_primer_info)
        right_to_exclude = self.identify_primers_to_exclude(right_primer_info)

        # Filter out primers to be excluded
        remaining_left_primers = {seq: info for seq, info in left_primer_info.items() if seq not in left_to_exclude}
        remaining_right_primers = {seq: info for seq, info in right_primer_info.items() if seq not in right_to_exclude}

        return remaining_left_primers, remaining_right_primers

    def check_for_dimers(self, primer1, primer2):
        """
        Check if two primers have potential for dimer formation.
        """
        dimer_potential = self.sliding_window_dimer_check(primer1, primer2)
        three_prime_complementarity = self.check_3prime_complementarity(primer1, primer2)
        return dimer_potential > 0 or three_prime_complementarity

    def identify_primers_to_exclude(self, primer_dimer_info):
        """
        Identify primers with the highest potential for forming dimers to be excluded.
        """
        # Sort primers based on the number of potential dimer partners
        sorted_primers = sorted(primer_dimer_info.items(), key=lambda x: len(x[1]), reverse=True)

        exclusion_threshold = 3
        return [primer_id for primer_id, dimer_partners in sorted_primers if len(dimer_partners) > exclusion_threshold]


# Define the primer design algorithm interface
class PrimerDesignAlgorithm(ABC):
    @abstractmethod
    def design_primers(self, sequences, output_file, common_primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params, specificity_check_database, num_primers=1, n_most_frequent=5):
        pass


# Design primers using the primer3 library
class Primer3Algorithm(PrimerDesignAlgorithm):
    def design_primers(self, sequences, output_file, common_primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params, specificity_check_database, num_primers=1, n_most_frequent=5):
        # The code to design primers using the primer3 library
        """
        Select primer targets within a given DNA sequence.

        An example of the Primer3 config file:

        [Primer3Parameters]
        PRIMER_OPT_SIZE = 20
        PRIMER_PICK_INTERNAL_OLIGO = 1
        PRIMER_INTERNAL_MAX_SELF_END = 8
        PRIMER_MIN_SIZE = 18
        PRIMER_MAX_SIZE = 25
        PRIMER_OPT_TM = 60.0
        PRIMER_MIN_TM = 57.0
        PRIMER_MAX_TM = 63.0
        PRIMER_MIN_GC = 20.0
        PRIMER_MAX_GC = 80.0
        PRIMER_MAX_POLY_X = 100
        PRIMER_INTERNAL_MAX_POLY_X = 100
        PRIMER_SALT_MONOVALENT = 50.0
        PRIMER_DNA_CONC = 50.0
        PRIMER_MAX_NS_ACCEPTED = 0
        PRIMER_MAX_SELF_ANY = 12
        PRIMER_MAX_SELF_END = 8
        PRIMER_PAIR_MAX_COMPL_ANY = 12
        PRIMER_PAIR_MAX_COMPL_END = 8

        [SpecificityParameters]
        IDENTITY_THRESHOLD = 95.0
        EVALUE_CUTOFF = 0.01

        """

        primers_for_all_sequences = {}

        for seq_id, sequence in sequences.items():
            # Sequence-specific parameters
            primer_design_params = common_primer_design_params.copy()
            primer_design_params['SEQUENCE_ID'] = seq_id
            primer_design_params['SEQUENCE_TEMPLATE'] = sequence

            primer_results = primer3.bindings.designPrimers(
                primer_design_params,
                {'PRIMER_NUM_RETURN': num_primers}
            )

            # Process the primer results
            primers = set()
            for i in range(num_primers):
                primer_pair = {
                    'left': primer_results[f'PRIMER_LEFT_{i}_SEQUENCE'],
                    'right': primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                    'left_tm': primer_results[f'PRIMER_LEFT_{i}_TM'],
                    'right_tm': primer_results[f'PRIMER_RIGHT_{i}_TM'],
                    'left_gc': primer_results[f'PRIMER_LEFT_{i}_GC_PERCENT'],
                    'right_gc': primer_results[f'PRIMER_RIGHT_{i}_GC_PERCENT']
                }
                primers.add(primer_pair)

            primers_for_all_sequences[seq_id] = primers


        # Summarizing the primers based on a given algorithm
        summarized_primers = primer_summarizing_algorithm.summarize(primers_for_all_sequences, primer_summarizing_params)

        # Check the specificity of the selected primers using the specified database
        specific_primers = specificity_check_algorithm.specificity_check(summarized_primers, specificity_check_database, primer_specificity_params)

        primer_dimers_check = PrimerDimersCheck()
        remaining_primers = primer_dimers_check.check_all_primer_combinations(specific_primers)

        return remaining_primers


class CustomAlgorithm(PrimerDesignAlgorithm):
    def design_primers(self, sequences, output_file, common_primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params, specificity_check_database, num_primers=1, n_most_frequent=5):
        # The code for a custom primer design algorithm
        pass


# Implementation of the strategy interface
class Strategy(ABC):
    @abstractmethod
    def design_primers(self, input_file, output_file, primer_design_algorithm: PrimerDesignAlgorithm, primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params, specificity_check_database):
        pass


# Implementation of concrete strategies for primer design
class TargetedAmpliconSequencingStrategy(Strategy):
    def design_primers(self, input_file, output_file, primer_design_algorithm: PrimerDesignAlgorithm, primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params, specificity_check_database):
        # Logic for targeted amplicon sequencing
        print("Designing primers for targeted amplicon sequencing...")

        return primer_design_algorithm.design_primers(input_file, output_file, primer_design_params, primer_summarizing_algorithm, primer_summarizing_params, specificity_check_algorithm, specificity_check_database)


class QuasialignmentStrategy(Strategy):
    class Segment:
        """
        A class to represent an individual segment in a quasi-alignment. It includes information about the sequence ID, species, 
        and the position of the segment within the sequence.
        """
        def __init__(self, seq_id: str, species: str, start: int, end: int):
            """
            Initialize a Segment object.

            :param seq_id: The ID of the sequence from which the segment comes.
            :param species: The species of the sequence.
            :param start: The starting position of the segment within the sequence.
            :param end: The ending position of the segment within the sequence.
            """
            self.seq_id = seq_id
            self.species = species
            self.start = start
            self.end = end

        def __repr__(self):
            return f"Segment(seq_id={self.seq_id}, species={self.species}, start={self.start}, end={self.end})"
            
    class QuasiAlignment:
        """
        A class to represent a quasi-alignment that contains multiple segments.
        Each segment is associated with a sequence from a particular species.
        """
        def __init__(self, cluster_id: str):
            """
            Initialize a QuasiAlignment object.

            :param cluster_id: A unique identifier for the quasi-alignment cluster.
            """
            self.cluster_id = cluster_id
            self.segments = []

        def add_segment(self, seq_id: str, species: str, start: int, end: int):
            """
            Add a segment to the quasi-alignment.

            :param seq_id: The ID of the sequence from which the segment comes.
            :param species: The species of the sequence.
            :param start: The starting position of the segment within the sequence.
            :param end: The ending position of the segment within the sequence.
            """
            segment = Segment(seq_id, species, start, end)
            self.segments.append(segment)

        def get_segments_by_species(self, species: str):
            """
            Get all segments that belong to a specific species.

            :param species: The species to filter segments by.
            :return: A list of segments that belong to the specified species.
            """
            return [segment for segment in self.segments if segment.species == species]

        def __repr__(self):
            return f"QuasiAlignment(cluster_id={self.cluster_id}, segments={self.segments})"

            
    def __init__(self):
        pass

    

    
class MarkerLociIdentificationStrategy(Strategy):
    def __init__(self):
        # Initialize the candidate_species_markers attribute as an empty DataFrame
        self.candidate_species_markers = pd.DataFrame(columns=['species', 'consensusSequence', 'start', 'end'])
        self.species_markers = pd.DataFrame(columns=['species', 'consensusSequence', 'start', 'end'])
        self.coverage_matrix = pd.DataFrame()
        # Initialize an empty list to store the selected markers
        self.selected_markers = []
        # Initialize a dictionary to indicate species without markers
        self.species_without_markers = {}

        self.db_params = {
            "dbname": "genomics_db",
            "user": user,
            "password": password,
            "host": "localhost"
        }

    def run_cactus(self, config_path, seq_file_paths, output_dir):
        """
        Runs the Cactus program with the given configuration, sequence files, and output directory.

        :param config_path: Path to the Cactus configuration file.
        :param seq_file_paths: A list of paths to the input sequence files.
        :param output_dir: The directory where the aligned genome will be saved.
        """
        # Join the sequence file paths into a space-separated string
        seq_files_str = " ".join(seq_file_paths)

        # Construct the Cactus command
        cactus_cmd = f"cactus {output_dir} {config_path} {seq_files_str} --realTimeLogging"

        try:
            # Execute the Cactus command
            subprocess.run(cactus_cmd, shell=True, check=True)
            print("Cactus alignment completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running Cactus: {e}")

    def sliding_window_analysis(self, genomes_dict, reference_fasta, output_folder, window_size, step_size, threshold=0.95):
        """
        Performs a sliding window analysis over genomic alignments.

        Parameters:
        - genomes_dict: Dictionary with keys as species names and values as genomic sequences.
        - reference_fasta: Path to the reference database in fasta format.
        - output_folder: Base folder for Minimap2 output.
        - window_size: Size of the window to test.
        - step_size: Step size for moving the window.
        - threshold: Threshold for evaluating species identification reliability.

        Returns:
        - List of tuples, each representing a tested window with (start_position, window_length, proportion_reliable_species).
        """
        windows_results = []
        for start_pos in range(0, len(genomes_dict), step_size):
            window_genomes = {species: seq[start_pos:start_pos + window_size] for species, seq in genomes_dict.items()}
            minimap2_output_folder = os.path.join(output_folder, f"window_{start_pos}_{window_size}")
            os.makedirs(minimap2_output_folder, exist_ok=True)

            # Run Minimap2 for each window
            run_minimap2(window_genomes, reference_fasta, minimap2_output_folder)

            # Process Minimap2 output and calculate identity, filter, and evaluate
            species_results = {}
            for species, _ in window_genomes.items():
                paf_file = os.path.join(minimap2_output_folder, f"{species}_alignment.paf")
                filtered_hits = calculate_identity_and_filter(paf_file, species)
                species_results[species] = filtered_hits

            proportion_reliable_species = evaluate_species_identification_proportion(species_results, threshold)

            # Save the window results
            windows_results.append((start_pos, window_size, proportion_reliable_species))

        return windows_results

    def initialize_db(self, user, password):
        # Parameters
        db_params = {
            "dbname": "postgres",
            "user": user,
            "password": password,
            "host": "localhost"
        }

        # Connect to the existing 'postgres' database
        conn = psycopg2.connect(**db_params)
        conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)

        # Create a new database
        cur = conn.cursor()
        cur.execute(sql.SQL("CREATE DATABASE {};").format(sql.Identifier('genomics_db')))

        cur.close()
        conn.close()

        # Update db_params to connect to the new database
        db_params["dbname"] = "genomics_db"

        # Connect to the new database
        conn = psycopg2.connect(**db_params)
        cur = conn.cursor()

        # SQL command to create a table
        create_table_command = """
        CREATE TABLE annotated_genes (
            gene_id SERIAL PRIMARY KEY,
            species VARCHAR(255),
            gene_name VARCHAR(255),
            sequence TEXT,
            annotation TEXT,
            category VARCHAR(255)  -- e.g., 'single-copy', 'core', 'accessory'
        );
        """

        # Execute the command
        cur.execute(create_table_command)
        conn.commit()

        # Clean up
        cur.close()
        conn.close()

    def run_prokka_and_panaroo(self, genome_paths, output_dir):
        """
        Annotates genomes with Prokka and performs pan-genome analysis with Panaroo.

        Parameters:
        - genome_paths: List of paths to genome files in FASTA format.
        - output_dir: Base directory for output files. Subdirectories will be created for Prokka and Panaroo outputs.
        """
        # Ensure output directories exist
        prokka_output_dir = os.path.join(output_dir, "prokka")
        panaroo_output_dir = os.path.join(output_dir, "panaroo")
        os.makedirs(prokka_output_dir, exist_ok=True)
        os.makedirs(panaroo_output_dir, exist_ok=True)

        # Run Prokka for each genome
        prokka_outputs = []  # To collect paths to Prokka output directories
        for genome_path in genome_paths:
            genome_base_name = os.path.basename(genome_path).replace('.fasta', '')
            prokka_out_dir = os.path.join(prokka_output_dir, genome_base_name)
            prokka_cmd = [
                'prokka', '--outdir', prokka_out_dir, '--prefix', genome_base_name,
                genome_path
            ]
            subprocess.run(prokka_cmd, check=True)
            prokka_outputs.append(prokka_out_dir)

        # Prepare Panaroo input (list of annotated genome directories)
        gff_paths = [os.path.join(dir, f"{os.path.basename(dir)}.gff") for dir in prokka_outputs]
        panaroo_cmd = ['panaroo', '-i'] + gff_paths + ['-o', panaroo_output_dir, '--clean-mode', 'strict']

        # Run Panaroo
        subprocess.run(panaroo_cmd, check=True)

    def save_annotated_genes_from_panaroo(self, panaroo_output_dir):
        """
        Saves annotated genes from Panaroo's output into the PostgreSQL database.

        Parameters:
        - panaroo_output_dir: Directory containing Panaroo's GFF output files.
        """
        # Connect to the database
        conn = psycopg2.connect(**self.db_params)
        cur = conn.cursor()

        # Path to the GFF file
        gff_file_path = os.path.join(panaroo_output_dir,
                                     "gene_data.csv")  # Panaroo outputs a gene_data.csv that contains the consolidated annotations

        # Parse the GFF file
        with open(gff_file_path, 'r') as file:
            for line in file:
                if line.startswith("gene_id"):  # Skip header
                    continue
                parts = line.strip().split(',')
                gene_name = parts[0]
                species = parts[2]
                sequence = parts[3]
                annotation = parts[4]
                category = "not specified"

                # Insert into the database
                insert_command = sql.SQL("""INSERT INTO annotated_genes 
                                            (species, gene_name, sequence, annotation, category) 
                                            VALUES (%s, %s, %s, %s, %s);""")
                cur.execute(insert_command, (species, gene_name, sequence, annotation, category))

        # Commit changes and close the connection
        conn.commit()
        cur.close()
        conn.close()

    def run_pgap(self, genome_paths, output_dir, pgap_version='2021-07-01.build5508', container_runtime='docker'):
        """
        Annotates genomes using the NCBI Prokaryotic Genome Annotation Pipeline (PGAP).

        Parameters:
        - genome_paths: List of paths to genome files in FASTA format.
        - output_dir: Base directory for output files. Subdirectories will be created for PGAP outputs.
        - pgap_version: The version of PGAP to use, corresponding to the PGAP image version.
        - container_runtime: The container runtime to use ('docker' or 'singularity').
        """
        # Ensure output directory exists
        pgap_output_dir = os.path.join(output_dir, "pgap")
        os.makedirs(pgap_output_dir, exist_ok=True)

        # Run PGAP for each genome
        for genome_path in genome_paths:
            genome_base_name = os.path.basename(genome_path).replace('.fasta', '')
            output_subdir = os.path.join(pgap_output_dir, genome_base_name)
            os.makedirs(output_subdir, exist_ok=True)

            # Prepare YAML input file for PGAP
            yaml_input_path = os.path.join(output_subdir, f"{genome_base_name}.yaml")
            self.prepare_pgap_input_yaml(genome_path, yaml_input_path)

            if container_runtime == 'docker':
                pgap_cmd = [
                    'docker', 'run', '--rm', '-v', f"{output_subdir}:/pgap/output", '-v',
                    f"{genome_path}:/pgap/input:ro",
                    f"ncbi/pgap:{pgap_version}", 'cwltool', '--outdir', '/pgap/output', '/pgap/pgap.cwl',
                    '/pgap/input/pgap_input.yaml'
                ]
            elif container_runtime == 'singularity':
                singularity_image_path = f"ncbi-pgap-{pgap_version}.sif"
                pgap_cmd = [
                    'singularity', 'exec', '--bind', f"{output_subdir}:/pgap/output", '--bind',
                    f"{genome_path}:/pgap/input:ro",
                    singularity_image_path, 'cwltool', '--outdir', '/pgap/output', '/pgap/pgap.cwl',
                    '/pgap/input/pgap_input.yaml'
                ]
            else:
                raise ValueError("Invalid container_runtime specified. Use 'docker' or 'singularity'.")

            # Run PGAP
            subprocess.run(pgap_cmd, check=True)

    def prepare_pgap_input_yaml(self, genome_path, yaml_input_path):
        """
        Prepares the input YAML file required by PGAP.

        Parameters:
        - genome_path: Path to the genome file in FASTA format.
        - yaml_input_path: Path where the input YAML file for PGAP should be saved.
        """

        yaml_content = f"""
    fasta:
      - {genome_path}
    submol:
      - name: "Example"
        strain: "Generic Strain"
    """
        with open(yaml_input_path, 'w') as file:
            file.write(yaml_content)

    def save_annotated_genes_from_pgap(self, pgap_output_file):
        """
        Saves annotated genes from PGAP's GenBank output into the PostgreSQL database.

        Parameters:
        - pgap_output_file: Path to PGAP's GenBank output file containing annotated genes.
        """
        # Connect to the database
        conn = psycopg2.connect(**self.db_params)
        cur = conn.cursor()

        # Parse the GenBank file
        with open(pgap_output_file) as file:
            for record in SeqIO.parse(file, "genbank"):
                species = record.annotations.get("organism", "Unknown species")
                for feature in record.features:
                    if feature.type == "gene" or feature.type == "CDS":
                        gene_name = feature.qualifiers.get("locus_tag", ["Unknown"])[0]
                        sequence = str(feature.extract(record.seq))
                        annotation = feature.qualifiers.get("product", [""])[0]
                        category = "not specified"

                        # Insert into the database
                        insert_command = sql.SQL("""INSERT INTO annotated_genes 
                                                    (species, gene_name, sequence, annotation, category) 
                                                    VALUES (%s, %s, %s, %s, %s);""")
                        cur.execute(insert_command, (species, gene_name, sequence, annotation, category))

        # Commit changes and close the connection
        conn.commit()
        cur.close()
        conn.close()

    def install_gatk(install_dir):
        """
        Download and install GATK (Genome Analysis Toolkit).

        Parameters:
        - install_dir: Path to the directory where GATK will be installed.
        """
        # Define the URL of the GATK release package
        gatk_url = "https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip"

        # Create the installation directory if it doesn't exist
        os.makedirs(install_dir, exist_ok=True)

        # Download the GATK package
        response = requests.get(gatk_url)
        with open(os.path.join(install_dir, "gatk.zip"), "wb") as f:
            f.write(response.content)

        # Extract the downloaded package
        with zipfile.ZipFile(os.path.join(install_dir, "gatk.zip"), "r") as zip_ref:
            zip_ref.extractall(install_dir)

        # Add execute permissions to the pgap-download tool
        pgap_download_path = os.path.join(install_dir, "gatk-4.5.0.0", "pgap-download")
        os.chmod(pgap_download_path, 0o755)

        subprocess.run(["export", f"PATH=$PATH:{os.path.join(install_dir, 'gatk-4.5.0.0')}"], shell=True)

    def download_pgap_database(self, output_dir):
            """
            Downloads the PGAP database from NCBI.

            Parameters:
            - output_dir: Directory where the downloaded PGAP database will be saved.
            """
            # Run the pgap-download tool to download the PGAP database
            try:
                subprocess.run(['pgap-download', '-o', output_dir], check=True)
                print("PGAP database downloaded successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error downloading PGAP database: {e}")

    def restore_pgap_database(self, database_name, dump_file_path):
        """
        Restore a PostgreSQL database from a dump file using pg_restore.

        Parameters:
        - database_name: The name of the PostgreSQL database to restore.
        - dump_file_path: The path to the dump file containing the database dump.
        """
        # Command to restore the database using pg_restore
        restore_command = [
            'pg_restore', '--dbname', database_name, '--verbose', dump_file_path
        ]

        try:
            # Execute the command
            subprocess.run(restore_command, check=True)
            print("Database restore completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error restoring database: {e}")

    def check_existing_pgap_annotation(self, db_params, genome_name):
        try:
            # Connect to the PGAP database
            conn = psycopg2.connect(**db_params)
            cur = conn.cursor()

            # Execute SQL query to search for records
            cur.execute("SELECT * FROM annotated_genomes WHERE genome_name = %s", (genome_name,))
            records = cur.fetchall()

            # Check if any records were found
            if records:
                print(f"Annotation found for genome '{genome_name}'.")
            else:
                print(f"No annotation found for genome '{genome_name}'.")

            # Clean up
            cur.close()
            conn.close()

        except (Exception, psycopg2.DatabaseError) as error:
            print("Error while connecting to PostgreSQL:", error)

    def search_genome(self, email, species):
        Entrez.email = email

        # Define the database and construct the query
        database = "genome"
        query = f"{species}[Orgn] AND annotated[Title]"

        # Perform the search
        with Entrez.esearch(db=database, term=query, retmax=10) as handle:
            search_results = Entrez.read(handle)

        # Return the search results
        return search_results

    def fetch_genome_details(self, id_list):
        with Entrez.efetch(db="genome", id=id_list, rettype="gb", retmode="text") as handle:
            details = handle.read()
        return details

    def download_genome_data(tax_id, output_directory='genome_data'):
        """
        Download genome data for a given taxonomic ID using NCBI Datasets.

        Parameters:
        tax_id (str): Taxonomic ID of the species to download genome data for.
        output_directory (str): Directory where the genome data will be saved.

        Returns:
        str: Path to the downloaded genome data file.
        """
        # Setup the API client
        configuration = ApiClient()

        # Initialize the Genome API
        api_instance = GenomeApi(configuration)

        try:
            # Download genome data zip file
            print(f"Downloading genome data for tax ID {tax_id}...")
            api_response = api_instance.download_assembly_package(
                taxon=tax_id,
                include_sequence=True,
                _preload_content=False
            )

            # Define the output file path
            output_path = f"{output_directory}/ncbi_dataset_{tax_id}.zip"

            # Write the response content to a file
            with open(output_path, "wb") as out_file:
                out_file.write(api_response.data)

            print(f"Genome data downloaded successfully: {output_path}")
            return output_path

        except Exception as e:
            print(f"An error occurred: {e}")
            return None

    def load_db_config(config_path):
        """Load database configuration from a JSON file."""
        with open(config_path, 'r') as file:
            config = json.load(file)
        return config

    def query_local_database(species_name, config_path):
        # Load database configuration from JSON file
        config = load_db_config(config_path)

        # Connect to your PostgreSQL database using the loaded credentials
        conn = psycopg2.connect(
            dbname=config['dbname'],
            user=config['user'],
            password=config['password'],
            host=config['host']
        )
        cursor = conn.cursor()

        # SQL query to find annotated genomes
        query = """
        SELECT * FROM genomes
        WHERE species = %s AND annotation_available = True
        """
        cursor.execute(query, (species_name,))

        # Fetch results
        results = cursor.fetchall()

        # Close the connection
        cursor.close()
        conn.close()

        return results

    def extract_species_from_filename(self, genome_path):
        """
        Extracts the species name from a genome file path, assuming the file name contains the genus and species
        separated by an underscore, followed by an optional strain identifier and the .fasta extension.

        Example filename: "Escherichia_coli_strain_xyz.fasta" -> "Escherichia coli"

        Parameters:
        - genome_path: The file path of the genome in FASTA format.

        Returns:
        - The species name extracted from the file name.
        """
        filename = os.path.basename(genome_path)
        # Remove the .fasta extension and split the remaining name by underscores
        parts = filename.replace('.fasta', '').split('_')
        # The species name consists of the first two parts (genus and species)
        species_name = ' '.join(parts[:2])
        return species_name

    def process_genome(self, genome_path, output_dir, accession, email, annotation_method='prokka', pgap_version='2021-07-01.build5508',
                       container_runtime='docker'):
        """
        Decides whether to download the genome or annotate it using Prokka and Panaroo or PGAP based on the existence of NCBI annotations.

        Parameters:
        - genome_path: Path to the genome file in FASTA format (or a list of paths).
        - output_dir: Base directory for output files.
        - annotation_method: Specify 'prokka' for Prokka+Panaroo or 'pgap' for PGAP.
        - pgap_version: The version of PGAP to use (if PGAP is chosen).
        - container_runtime: The container runtime to use ('docker' or 'singularity') if PGAP is chosen.
        """
        species_name = self.extract_species_from_filename(genome_path)
        if check_ncbi_annotations(species_name):
            download_genome_annotation(accession, email, output_dir)
        else:
            print(f"No annotations found for {species_name}. Proceeding with local annotation.")
            if annotation_method == 'prokka':
                run_prokka_and_panaroo(genome_path, output_dir)
            elif annotation_method == 'pgap':
                run_pgap(genome_path, output_dir, pgap_version, container_runtime)
            else:
                raise ValueError(f"Invalid annotation method: {annotation_method}")

    def extract_single_copy_genes(self, genome_paths, output_dir):
        # Load the gene_presence_absence.csv file
        gpa_df = pd.read_csv('path/to/panaroo/output/gene_presence_absence.csv')

        # Filter for single-copy core genes
        presence_columns = gpa_df.columns[4:]
        single_copy_core_genes = gpa_df.loc[(gpa_df[presence_columns] != '').sum(axis=1) == len(presence_columns)]

        # Save the filtered genes to a new file
        single_copy_core_genes.to_csv('path/to/output/single_copy_core_genes.csv', index=False)

    def runMuscle(self,input_file,output_file):
        """
        Run muscle on an input file, creating an output file and return stdout and the stderr
        """

        # Command to run muscle
        cmd = ["muscle", "-align", input_file, "-output", output_file]

        try:
            # Run the command and capture stdout and stderr
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            stdout = stdout.decode('utf-8')
            stderr = stderr.decode('utf-8')

            return stdout, stderr

        except subprocess.CalledProcessError as e:
            return None, str(e)
        
    def runClustalw(self,input_file,output_file): 
        """
        Run clustalw on an input file, creating an output file and return stdout and the stderr
        """

        # Command to run clustalw
        cmd  = ['clustalw', '-INFILE=' + input_file, '-OUTFILE=' + output_file]

        try:
            # Run the command and capture stdout and stderr
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            stdout = stdout.decode('utf-8')
            stderr = stderr.decode('utf-8')

            return stdout, stderr

        except subprocess.CalledProcessError as e:
            return None, str(e)
        
    def runMafft(self,input_file,output_file):

        """
        Run mafft on an input file, creating an output file and return stdout and the stderr
        """

        # Command to run mafft
        cmd  = ['mafft', input_file, '>', output_file]

        try:
            # Run the command and capture stdout and stderr
            process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            stdout = stdout.decode('utf-8')
            stderr = stderr.decode('utf-8')

            return stdout, stderr

        except subprocess.CalledProcessError as e:
            return None, str(e)

    def align_sequences(self, input_file, alignment_tool='prank', output_format='fasta'):
        """
        Align sequences using a specified alignment tool.

        Parameters:
        - input_file: Path to the input FASTA file.
        - alignment_tool: The alignment tool to use ('clustalw', 'muscle', 'mafft', 'prank').
        - output_format: Format of the output alignment.

        Returns:
        - Path to the output alignment file.
        """

        output_file = f"{input_file.rsplit('.', 1)[0]}_{alignment_tool}.aln"

        if alignment_tool.lower() == 'clustalw':
            #clustalw_cline = ClustalwCommandline("clustalw2", infile=input_file, outfile=output_file)
            clustalw_cline = self.runClustalw(nfile=input_file, outfile=output_file)
            stdout, stderr = clustalw_cline()

        elif alignment_tool.lower() == 'muscle':
            #muscle_cline = MuscleCommandline(input=input_file, out=output_file)
            muscle_cline = self.runMuscle(nfile=input_file, outfile=output_file)
            stdout, stderr = muscle_cline()

        elif alignment_tool.lower() == 'mafft':
            # mafft_cline = MafftCommandline(input=input_file)
            # stdout, stderr, = subprocess.Popen(str(mafft_cline),
            #                                    stdout=subprocess.PIPE,
            #                                    stderr=subprocess.PIPE,
            #                                    shell=True,
            #                                    text=True).communicate()
            mafft_cline = self.runMafft(nfile=input_file, outfile=output_file)
            stdout, stderr = mafft_cline()
            with open(output_file, 'w') as f:
                f.write(stdout)

        elif alignment_tool.lower() == 'prank':
            prank_cline = PrankCommandline(d=input_file, o=output_file.rsplit('.', 1)[0])
            prank_cline()

        else:
            raise ValueError("Unsupported alignment tool. Choose 'clustalw', 'muscle', 'mafft', or 'prank'.")

        return output_file

    def save_single_copy_genes_to_fasta(self, output_file):
        conn = psycopg2.connect("dbname=genomics_db user=your_username")
        cur = conn.cursor()
        cur.execute("SELECT species, gene_name, sequence FROM annotated_genes WHERE category = 'single-copy'")

        with open(output_file, 'w') as fasta_file:
            for record in cur.fetchall():
                species, gene_name, sequence = record
                fasta_file.write(f'>{species}_{gene_name}\n{sequence}\n')

        cur.close()
        conn.close()

    def alignment_to_dataframe(self, alignment_file, format='fasta'):
        """
        Converts an alignment file into a simple Pandas DataFrame with two columns:
        'SeqName' for the sequence name and 'Seq' for the sequence itself.

        Parameters:
        - alignment_file: Alignment - Path to the alignment file.
        - format: Format of the alignment file (default 'fasta').

        Returns:
        - A Pandas DataFrame with 'SeqName' and 'Seq' columns.
        """
        # Read the alignment
        alignment = AlignIO.read(alignment_file, format)

        # Initialize an empty list to store sequence name and sequence data
        sequence_data = []

        # Loop through each record in the alignment
        for record in alignment:
            # Append a tuple with the sequence name (identifier) and the sequence
            sequence_data.append((record.id, str(record.seq)))

        # Create a DataFrame from the sequence data
        df = pd.DataFrame(sequence_data, columns=['SeqName', 'Seq'])

        return df

    def get_consensus_sequence_from_alignment(self, alignment_file, consensus_threshold=0.7):
        """
        Generates a consensus sequence from a given alignment file.

        Parameters:
        - alignment_file: Path to the alignment file (e.g., in FASTA format).
        - consensus_threshold: Threshold for determining the consensus (default is 0.7).

        Returns:
        - A consensus sequence as a string.
        """
        # Load the alignment from the file
        alignment = AlignIO.read(alignment_file, "fasta")

        # Create a summary object from the alignment
        summary_align = AlignInfo.SummaryInfo(alignment)

        # Generate the consensus sequence based on the given threshold
        consensus = summary_align.dumb_consensus(threshold=consensus_threshold, ambiguous="N")

        return str(consensus)

    def find_conserved_regions_consensus(self, alignment_file, threshold=0.7, min_length=5):
        """
        Identifies conserved regions in a set of sequences based on a consensus sequence from an existing alignment,
        including only those regions that meet a specified minimum length.

        Parameters:
        - alignment_file: Path to the alignment file.
        - threshold: Proportion of sequences that must contain a nucleotide at a position for it to be considered conserved.
        - min_length: Minimum length of conserved regions to be included in the result.

        Returns:
        - A tuple of two lists:
            1. The first list is a list of lists, where each inner list represents the positions of one conserved region.
            2. The second list is a list of lists of dominant letters for each conserved region.
        """
        # Read the alignment
        alignment = AlignIO.read(alignment_file, "fasta")

        # Generate a consensus sequence
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus_info = summary_align.summary()

        conserved_regions_positions = []
        conserved_regions_letters = []
        current_region_positions = []
        current_region_letters = []

        for i, base_dict in enumerate(consensus_info.items(), start=1):  # Positions are 1-based
            total_bases = sum(base_dict[1].values())
            dominant_bases = {base: count for base, count in base_dict[1].items() if count / total_bases >= threshold}

            if dominant_bases:
                current_region_positions.append(i)
                current_region_letters.append(list(dominant_bases.keys()))
            else:
                if current_region_positions and len(current_region_positions) >= min_length:
                    conserved_regions_positions.append(current_region_positions)
                    conserved_regions_letters.append(current_region_letters)
                    current_region_positions = []
                    current_region_letters = []

        # Check if there's an unclosed region at the end that meets the minimum length requirement
        if current_region_positions and len(current_region_positions) >= min_length:
            conserved_regions_positions.append(current_region_positions)
            conserved_regions_letters.append(current_region_letters)

        return conserved_regions_positions, conserved_regions_letters

    def shannon_entropy(self, letter_count_dict, base=2):
        """
        Calculate the Shannon Entropy based on a dictionary of the letters occurring at a given position and the amount of occurrences
        """

        total_sum = sum(letter_count_dict.values())
        for key in letter_count_dict:
            letter_count_dict[key] = (letter_count_dict[key] / total_sum) * math.log(letter_count_dict[key] / total_sum,
                                                                                     base)
        entropy = - sum(letter_count_dict.values())

        return entropy

    def find_conserved_regions_shannon_entropy(self, sequence_df, entropy_thr=0.2, length_thr=3):
        """
        Find conserved regions with set minimal length based on a entropy threshold
        """
        # find out how long the alignments are
        valuesList = list(sequence_df["Seq"])
        alignLength = len(valuesList[0])

        # initiate lists
        conserved_regions_positions = list()
        conserved_regions_dominant = list()
        current_conserved_region = list()
        current_conserved_region_letters = list()
        current_conserved_scores = list()
        current_dominant_gap_count = 0

        # calculate Shannon Entropy for all positions and track the conserved regions
        i = 0
        while i in range(0, alignLength):
            # count the different letters at each position
            letter_count_dict = {}
            for sequence in valuesList:
                letter_count_dict[sequence[i]] = letter_count_dict.get(sequence[i], 0) + 1
            # calculate Shannon Entropy
            entropy = shannon_entropy(letter_count_dict)

            # identify the dominant letter(s)
            max_value = max(letter_count_dict.values())
            dominant_letters = [key for key, value in letter_count_dict.items() if value == max_value]

            # if the calculated entropy is smaller than the threshold, initiate a new conserved region or append to the current one
            if entropy <= entropy_thr:
                current_conserved_region.append(i)
                current_conserved_region_letters.append(dominant_letters)
                current_conserved_scores.append(entropy)
                if "-" in dominant_letters:
                    current_dominant_gap_count = current_dominant_gap_count + 1
            # if the calculated entropy is bigger than the threshold, but we currently are in a conserved region, calculate the new average Shannon's Entropy
            elif len(current_conserved_region) > 0:
                current_conserved_scores.append(entropy)
                average = sum(current_conserved_scores) / len(current_conserved_scores)
                # continue the conserved region if it is smaller than the threshold
                if average <= entropy_thr:
                    current_conserved_region.append(i)
                    current_conserved_region_letters.append(dominant_letters)
                    if "-" in dominant_letters:
                        current_dominant_gap_count = current_dominant_gap_count + 1
                # otherwise end the conserved region and add it to the list of conserved regions if it meets the length threshold
                else:
                    if len(current_conserved_region) - current_dominant_gap_count >= length_thr:
                        conserved_regions_positions.append(current_conserved_region)
                        conserved_regions_dominant.append(current_conserved_region_letters)
                    current_conserved_region = list()
                    current_conserved_region_letters = list()
                    current_dominant_gap_count = 0

            i = i + 1

        # add the last current region to the list of conserved regions if it meets the length threshold
        if len(current_conserved_region) > 0:
            conserved_regions_positions.append(current_conserved_region)
            conserved_regions_dominant.append(current_conserved_region_letters)

        return conserved_regions_positions, conserved_regions_dominant

    def manhattan_distance(vector1,vector2):
        """
        Calculate the manhattan distance between two vectors
        """
        if len(vector1)==len(vector2):
            distance=0
            for i in range(0,len(vector1)):
                distance = distance + (abs(vector1[i]-vector2[i]))
                
        else:
            print(f"error: vectors do not have the same length")
        return distance
    
    def calculate_centroid_mean(vectors):
        """
        Calculate the centroid as the mean of all vectors in a cluster
        """
        mean_vector = [sum(values) / len(vectors) for values in zip(*vectors)]
        return mean_vector

    def find_clusters(segment_df,metric,threshhold):
        """
        Cluster the segments based on their metric distance
        """
        only_add_to_best_match = False
        # todo: maybe different way of determining which cluster to add this too if two clusters are equally good.

        # The cluster centers will hold an average vector for each cluster and be used to determine the distance of each segment to each cluster. 
        clusters = []
        cluster_centers = []
        
        for genome_id in segment_df['genome_id'].unique():

            # For each of the different genome IDs, look at each vector representing a sequence segment summary vector and either create a new cluster for it
            # or or put it inot one of the existing clusters
            for rowWithId in segment_df[segment_df['genome_id'] == genome_id].index:
                vector = segment_df.iloc[rowWithId]['vector']

                if len(cluster_centers) == 0:
                    # Create a new cluster for the first segment
                    clusters.append([vector])
                    # A cluster with only one sequence has that sequence as its center
                    cluster_centers.append(vector)
                    # The information which cluster we added the segment to is added to the df
                    segment_df.at[rowWithId, 'cluster'] = [0]
                    
                else:
                    if only_add_to_best_match:
                        # Calculate the metric distances of the segment to each center of a cluster
                        metric_distances = []
                        for center in cluster_centers:
                            if metric == "manhattan":
                                metric_distances.append(manhattan_distance(center,vector))
                            else:
                                raise Exception(f"unknown metric {metric} selected for conserved region identification based on quasi alignments")

                        if min(metric_distances) <= threshhold:
                            min_dist_cluster_index = metric_distances.index(min(metric_distances))
                            clusters[min_dist_cluster_index].append(vector)
                            cluster_centers[min_dist_cluster_index] = calculate_centroid_mean(clusters[min_dist_cluster_index])
                            segment_df.at[rowWithId, 'cluster'] = [min_dist_cluster_index]

                        else:
                            clusters.append([vector])
                            cluster_centers.append(vector)
                            segment_df.at[rowWithId, 'cluster'] = [len(clusters)]

                    else:
                        #todo: try to instead add to any cluster where the distance is smaller than the threshold and see if it is better
                        clusters_under_thr = []

                        for center in cluster_centers:
                            if metric == "manhattan":
                                dist = manhattan_distance(center,vector)
                            else:
                                raise Exception(f"unknown metric {metric} selected for conserved region identification based on quasi alignments")
                            if dist <= threshhold:
                                cluster_index = cluster_centers.index(center)
                                clusters[cluster_index].append(vector)
                                cluster_centers[cluster_index] = calculate_centroid_mean(clusters[cluster_index])
                                clusters_under_thr.append(cluster_index)

                        if len(clusters_under_thr) > 0:
                            segment_df.at[rowWithId, 'cluster'] = clusters_under_thr
                        else:
                            clusters.append([vector])
                            cluster_centers.append(vector)
                            segment_df.at[rowWithId, 'cluster'] = [len(clusters)]

        return segment_df

    def trim_sequence(source_df,seq_id,beginning,end,continuation_len_thr,segment_size):
        """
        based on beginning and end points, get a new sequence fragment and add a few bases at the beginning and end (determined by continuation_len_thr) if in bound of the sequence
        """
        seq_beginning = beginning - round(0.5*continuation_len_thr*segment_size) if (beginning - round(0.5*continuation_len_thr*segment_size) > 0) else 0
        seq_end = end + round(0.5*continuation_len_thr*segment_size) if (end + round(0.5*continuation_len_thr*segment_size) < len(source_df.loc[source_df["genome_id"]==seq_id]["sequence"].iloc[0])) else len(source_df.loc[source_df["genome_id"]==seq_id]["sequence"].iloc[0])
        sequence = source_df.loc[source_df["genome_id"]==seq_id]["sequence"].iloc[0][seq_beginning:seq_end]

        return sequence
                        
    def find_candidate_regions(segment_df, source_df, segment_size, conservation_thr,continuation_len_thr):
        """
        Identify the candidates for conserved regions based on the clusters.
        In this step, several segments of the same sequence within one cluster are merged if they overlap or are close together.
        """
        # Extract the numbers of the clusters as well as the total number of different sequences
        cluster_numbers = list(set([number for list in segment_df['cluster'] for number in list]))
        number_of_sequences = len(segment_df['genome_id'].unique())

        # Initiate an empty data frame to write the information to
        columns = ['cluster','genome_id', 'sequence', 'beginning', 'end']
        candidate_conserved_region_df = pd.DataFrame(columns=columns)

        for cluster_n in cluster_numbers:
            # For each cluster, we consider all the segments that are in the cluster
            current_cluster = segment_df[segment_df["cluster"].apply(lambda list: cluster_n in list)]
            # We calculate a consensus score that tells us how big a fraction of the different sequences has at least one segment in the cluster
            cluster_origins = current_cluster['genome_id'].unique()
            consensus_score = round(len(cluster_origins)/number_of_sequences,2)

            if consensus_score >= conservation_thr:
                # If the consensus sore is smaller than the threshold, we add the cluster as a candidate conserved region
                for seq_id in cluster_origins:
                    beginning = current_cluster[current_cluster["genome_id"]==seq_id]['beginning'].values
                    end = current_cluster[current_cluster["genome_id"]==seq_id]['end'].values

                    if(isinstance(beginning,np.ndarray)):
                        # If more than just one segment per sequence was added, we want to either consider them as one large segment 
                        # if they are overlapping or close we just adjust beginning and end, otherwise we ?
                        beginning_sorted = sorted(beginning)
                        end_sorted = sorted(end)

                        new_beginning = beginning[0] 
                        new_end = end[0]
                        for i in range(1,len(beginning_sorted)):
                            if (beginning_sorted[i]-beginning_sorted[i-1]-segment_size > continuation_len_thr*segment_size):
                                # todo: for now i just added both of them but i might want to come up with a better method... one possibility: in this case we want to align only one of the segments, we take the one that is closer to the center
                                
                                # the two segments are too far apart to be united so we add the first one to the df and continue looking at the next segment
                                sequence = trim_sequence(source_df,seq_id,new_beginning,new_end,continuation_len_thr,segment_size)
                                new_row = pd.DataFrame({'cluster' : [cluster_n],'genome_id' : [seq_id],'sequence' : [sequence], 'beginning': [new_beginning], 'end': [new_end]})
                                candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

                                # and set the new beginning and end to the new segment
                                new_beginning = beginning[i] 
                                new_end = end[i]
                            else:
                                # we want to unite the segments, so we overwrite the "new end" and continue looking at the next segment
                                new_end = end_sorted[i]
                        # after looking at each segment, we add the last one to our df
                        sequence = trim_sequence(source_df,seq_id,new_beginning,new_end,continuation_len_thr,segment_size)
                        new_row = pd.DataFrame({'cluster' : [cluster_n],'genome_id' : [seq_id],'sequence' : [sequence], 'beginning': [new_beginning], 'end': [new_end]})
                        candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
                    
                    else:
                        # In cases where we only have one segment for a sequence, we sinply add the segment to the df:
                        # We let the sequence be a bit longer than the candidate region so we can find the actual beginning of the conserved region in the alignment
                        sequence = trim_sequence(source_df,seq_id,beginning,end,continuation_len_thr,segment_size)
                        new_row = pd.DataFrame({'cluster' : [cluster_n],'genome_id' : [seq_id],'sequence' : [sequence], 'beginning': [beginning], 'end': [end]})
                        candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
            
        return candidate_conserved_region_df

    def find_overlapping_conserved_regions(candidate_conserved_region_df: pd.DataFrame, source_df, clusters_to_consider: list[int],clusters_in_final_df: list[int],continuation_len_thr,conservation_thr,number_of_sequences,segment_size):
        """
        Check if the candidate regions overlap or are very close to each other in any of the sequences and if so, combine them
        """
        for cluster_num in clusters_to_consider:

            cluster_of_interest = candidate_conserved_region_df.loc[candidate_conserved_region_df["cluster"]==cluster_num]

            edited = False

            # For now we only consider clusters that have only one segment per sequence, todo: is this wise?
            if len(set(cluster_of_interest["genome_id"]))==len(cluster_of_interest["genome_id"]):

                # We compare the cluster to each other cluster we still need to consider
                for cluster_to_compare_to in clusters_to_consider:
                    if (cluster_num != cluster_to_compare_to):
                        all_seq_distances:list = []
                        sequences_in_cluster_to_compare_to = candidate_conserved_region_df.loc[candidate_conserved_region_df["cluster"] == cluster_to_compare_to]["genome_id"]
                        # For now we only consider clusters that have only one segment per sequence, todo: is this wise?
                        if len(set(sequences_in_cluster_to_compare_to))==len(sequences_in_cluster_to_compare_to):
                            sequences_in_both_sets = set(cluster_of_interest["genome_id"]).intersection(set(sequences_in_cluster_to_compare_to))
                            # if a sufficient number of sequences in in both clusters, we check if the sequences overlap or are close to each other
                            if (len(sequences_in_both_sets) >= conservation_thr*number_of_sequences):
                                for seq_id in sequences_in_both_sets:
                                    beginning1 = candidate_conserved_region_df.loc[(candidate_conserved_region_df["cluster"] == cluster_to_compare_to)&(candidate_conserved_region_df["genome_id"] == seq_id)]["beginning"].values
                                    beginning0 = cluster_of_interest.loc[cluster_of_interest["genome_id"] == seq_id]["beginning"].values
                                    distance = min(abs(beginning1[0] - beginning0[0]),abs(beginning0[0] - beginning1[0]))
                                    all_seq_distances.append(distance)
                                averageDistance = sum(all_seq_distances)/len(all_seq_distances)
                                if abs(averageDistance) <= continuation_len_thr:
                                    # If this is the case, we merge the segments for each sequence and add them to the df
                                    for seq_id in sequences_in_both_sets:
                                        beginning0 = int(min(cluster_of_interest.loc[cluster_of_interest["genome_id"] == seq_id]["beginning"].values))
                                        beginning1 = int(min(candidate_conserved_region_df.loc[(candidate_conserved_region_df["cluster"] == cluster_to_compare_to)&(candidate_conserved_region_df["genome_id"] == seq_id)]["beginning"].values))
                                        new_beginning = int(min(beginning0,beginning1))
                                        end0 = int(max(cluster_of_interest.loc[cluster_of_interest["genome_id"] == seq_id]["end"].values))
                                        end1 = int(max(candidate_conserved_region_df.loc[(candidate_conserved_region_df["cluster"] == cluster_to_compare_to)&(candidate_conserved_region_df["genome_id"] == seq_id)]["end"].values))
                                        new_end = int(max(end0,end1))
                                        new_sequence = trim_sequence(source_df,seq_id,new_beginning,new_end,continuation_len_thr,segment_size)
                                        new_cluster_name = max(candidate_conserved_region_df['cluster'])+1
                                        new_row = pd.DataFrame({'cluster' : [new_cluster_name],'genome_id' : [seq_id],'sequence' : [new_sequence], 'beginning': [new_beginning], 'end': [new_end]})
                                        candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
                                    # remove the clusters from the list of clusters we still need to consider and add the new one instead 
                                    clusters_to_consider=np.append(clusters_to_consider[(clusters_to_consider != cluster_num)&(clusters_to_consider != cluster_to_compare_to)],new_cluster_name)
                                    edited = True
                                    (candidate_conserved_region_df,clusters_in_final_df) = find_overlapping_conserved_regions(candidate_conserved_region_df,source_df,clusters_to_consider,clusters_in_final_df,continuation_len_thr,conservation_thr,number_of_sequences,segment_size)
                                    break

            if not edited :
                clusters_to_consider=clusters_to_consider[clusters_to_consider != cluster_num]
                clusters_in_final_df = np.append(clusters_in_final_df,cluster_num)
            else: break

        return (candidate_conserved_region_df,clusters_in_final_df)

    def align_candidate_regions(candidate_conserved_region_df):
        """
        Align the sequences of the conserved region candidates
        """
        conserved_regions_positions = []
        conserved_regions_dominant = []
        cluster_nums = candidate_conserved_region_df['cluster'].unique()
        for i in cluster_nums:
            sequences_to_consider = candidate_conserved_region_df.loc[candidate_conserved_region_df["cluster"]==i]
            sequence_tuples:list = []
            for j in range(0,len(sequences_to_consider)):
                (genome_id,sequence) = (sequences_to_consider.iloc[j]["genome_id"],sequences_to_consider.iloc[j]["sequence"])
                sequence_tuples.append((genome_id,sequence))

            with NamedTemporaryFile(delete=False) as temp_file:
                for pair in sequence_tuples:
                    # Write the sequences to a temporary FASTA file
                    genome_in = pair[0]
                    sequence = pair[1]
                    temp_file.write(f">{genome_in} alignment_cluster_{i}_{genome_in}\n{sequence}\n".encode())
                temp_file.flush()
            alignment_path = align_sequences(temp_file.name)
            df = f2p.read_multifasta(alignment_path)

            conservedRegion = se.find_conserved_regions_shannon_entropy(df)
            conserved_regions_positions = conserved_regions_positions + conservedRegion[0]
            conserved_regions_dominant = conserved_regions_dominant + conservedRegion[1]
        
        return (conserved_regions_positions,conserved_regions_dominant)

    def quasi_align(source_df,threshold,overlap,conservation_thr,additional_sequence_symbols:list[str]=[],metric:str="manhattan",length_thr=0.8,segment_size:int=100,continuation_len_thr:int=2):
        """
        Find conserved regions without aligning sequences based on a clustering method.
        The sequences are divided into segments of a given length and the number of occurrences of the different possible base triplets is used to 
        determine the distance between two segments based on a given Metric. Similar segments are then clustered to find candidates for conserved regions.

        The parameters used here are:
        - segment_size: The length of segments the sequences are divided into (default: 100)
        - threshold: The maximum allowed metric distance that segments can have and still be clustered together
        - overlap: When dividing the sequence into segments to compare triplet occurrences in, an overlap of segments can be allowed and is given as 1/n to indicate the fraction of segments at which the overlap should begin.
        - conservation_thr: Sets the threshold for the fraction of all sequences that have to be in one cluster in order to consider the cluster as a potential conserved region
        - length_thr: Determines the minimal length that a segment can have relative to the segment size (default: 0.8)
        - continuation_len_thr: distance relative to segment size that two segments can have to be considered one continuous segment (default: 1)
        """
        # Todo: different distance for specific symbols? if e.g. n can be a or t then the distance of a n to these should be smaller than to g or c
        # Todo: add words beginning from the reverse end if the segment is to short? to make sure every part of the sequence is actually covered
        # Todo: find a better way to handle overlapping conserved region candidates. Right now we do not merge conserved region candidates if they have more than just one segment per sequence but in cases where we have very similar sequences this can lead to overlapping regions not being merged and we end up aligning more things than we would by just aligning the normal sequences
        # Todo: come up with good way to determine default parameters e.g. based on a few tests with the summary vectors. (maybe by calculating distances within and between a few sequences)

        # Generate all possible words of a given length (by default triplets) that we will later search for in the sequences. 
        alphabet = ["A","T","G","C"] + additional_sequence_symbols
        word_length = 3
        possible_words = []
        def generate(word:str):
            if (len(word)-word_length) == 0:
                possible_words.append(word)
                return
            else: 
                for base in alphabet:
                    generate(word + base)
        generate("")

        # Initiate a data frame to store the information on segments
        columns = ['genome_id', 'beginning', 'end', 'sequence', 'vector',"cluster"]
        segment_df = pd.DataFrame(columns=columns)

        # Get the number of rows in the data frame containing the source sequences
        num_samples = source_df.shape[0]

        # Divide the sequences into segments, count triplet occurrences and store information in the new dataframe
        for row in range(0,num_samples):

            # Extract nformation from source df
            sample = source_df.iloc[row]["sequence"]
            genome_id = source_df.iloc[row]["genome_id"]

            for i in range(0, len(sample), round(segment_size*overlap)):

                # Extract the actual sequence for each segment
                segment = sample[i:(i + segment_size)]

                if len(segment) >=  length_thr*segment_size:

                    # Construct a summary vector for triplet occurrences for each of the segments
                    vector = []
                    for word in possible_words:
                        vector.append(segment.count(word))

                    # Add the new information to the new df
                    new_row = pd.DataFrame({'genome_id': genome_id, 'beginning' : i, 'end' : i + segment_size, 'sequence' : segment, 'vector' : [vector], "cluster" : [np.nan]})
                    segment_df = pd.concat([segment_df, new_row], ignore_index=True)

        # Find clusters and update the segment_df 
        segment_df = find_clusters(segment_df,metric,threshold)

        # Find candidate clusters for conserved regions
        candidate_regions_df = find_candidate_regions(segment_df, source_df, segment_size, conservation_thr,continuation_len_thr)

        # Merge overlapping regions:
        (extended_candidate_regions_df,clusters_in_final_df) = find_overlapping_conserved_regions(candidate_regions_df, source_df, candidate_regions_df["cluster"].unique(), [],continuation_len_thr,conservation_thr,len(source_df),segment_size)

        candidate_regions_df = extended_candidate_regions_df[extended_candidate_regions_df["cluster"].isin(clusters_in_final_df)]

        #align the candidate regions and calculate their conservation score
        (conserved_regions_positions,conserved_regions_dominant) = align_candidate_regions(candidate_regions_df)

        return(conserved_regions_positions, conserved_regions_dominant)

    def extract_consensus_sequence_and_conserved_regions(self, folder_path):
        results = {}
        for file_name in os.listdir(folder_path):
            if file_name.endswith(".aln"):
                species_name = "_".join(file_name.split("_")[:-1])
                alignment_file = os.path.join(folder_path, file_name)
                consensus_sequence = self.get_consensus_sequence_from_alignment(alignment_file)

                alignment_df = alignment_to_dataframe(alignment_file)

                conserved_regions, dominant_letters = self.find_conserved_regions_shannon_entropy(alignment_df)
                results[species_name] = (consensus_sequence, conserved_regions, dominant_letters)

        return results

    def extract_conserved_region_sequences(self, species_data, species_name):
        """
        Creates a dictionary with individual conserved regions as keys (including the species name and a two-digit index)
        and the actual sequences of these regions as values for a selected species. The keys are formatted as
        "speciesName_XX_start-end" for each conserved region, where XX is the two-digit index of the region.

        Parameters:
        - species_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - species_name: The name of the species to process.

        Returns:
        - A dictionary with conserved regions (identified by species name, index, start and end positions)
          as keys and sequences for these regions as values.
        """
        conserved_regions_dict = {}
        if species_name in species_data:
            consensus_sequence, conserved_regions_positions, _ = species_data[species_name]

            for idx, region_positions in enumerate(conserved_regions_positions, start=1):
                start, end = region_positions[0], region_positions[-1]  # Extract start and end from the positions list
                # Extract the sequence part corresponding to the current conserved region using slicing
                region_sequence = consensus_sequence[start - 1:end]
                # Create a unique key for each region including the species name, a two-digit index, and its start and end positions
                region_key = f"{species_name}_{idx:02d}_{start}-{end}"
                conserved_regions_dict[region_key] = region_sequence
        else:
            print(f"Species '{species_name}' not found in the data.")

        return conserved_regions_dict

    def conserved_regions_to_dataframe(self, conserved_regions_info):
        """
        Processes information about conserved regions and stores it in the candidate_species_markers attribute as a DataFrame.

        Parameters:
        - conserved_regions_info: A dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        """
        # Temporary list to store data before converting to DataFrame
        temp_data = []

        for species, (consensus_sequence, regions_positions, _) in conserved_regions_info.items():
            for idx, region_positions in enumerate(regions_positions, start=1):
                start, end = region_positions[0], region_positions[-1]
                # Append the data to the temp_data list
                temp_data.append(
                    {'species': species, 'index': idx, 'consensusSequence': consensus_sequence, 'start': start, 'end': end})

        # Convert the temporary list of data into a DataFrame
        new_data = pd.DataFrame(temp_data)

        # Store this new DataFrame in the candidate_species_markers attribute
        self.candidate_species_markers = pd.concat([self.candidate_species_markers, new_data], ignore_index=True)

    def extract_sequences_and_run_minimap2(self, species_conservation_data, reference_db):
        """
        Loops through species and their conserved regions, extracts corresponding parts of the consensus sequence,
        and runs minimap2 against a reference database for these sequences.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - reference_db: Path to the reference database used by minimap2.
        """
        for species, (consensus_seq, conserved_regions, _) in species_conservation_data.items():
            for idx, region_positions in enumerate(conserved_regions, start=1):
                start_end_tuples = [(region[0], region[-1]) for region in region_positions if region]
                for start, end in start_end_tuples:
                    # Extract the sequence part corresponding to the current conserved region
                    region_seq = consensus_seq[start:end + 1]

                    # Create a temporary file to store the region sequence for minimap2 input
                    temp_seq_file = f"temp_{species}_{start}_{end}.fasta"
                    with open(temp_seq_file, "w") as f:
                        f.write(f">{species}_{start}_{end}\n{region_seq}\n")

                    # Define the output folder dynamically based on species and region
                    output_folder_specific = os.path.join("marker_loci/alignment/minimap2", species, f"_{idx:02d}_{start}_{end}_alignment.paf")
                    os.makedirs(output_folder_specific, exist_ok=True)

                    # Run minimap2 for the extracted sequence part against the reference database
                    self.run_minimap2(self.alignment_to_dataframe(temp_seq_file), reference_db, output_folder=output_folder_specific)

                    # Build the minimap2 command
                    minimap2_cmd = ["minimap2", reference_db, temp_fasta_filename, "-o", output_file]

                    # Run minimap2
                    subprocess.run(minimap2_cmd)

                    # Remove the temporary sequence file after minimap2 run
                    os.remove(temp_seq_file)

    def calculate_identity_and_filter(self, paf_folder_path, identity_threshold=98):
        """
        Processes all PAF files in a specified folder, calculating the identity of alignments,
        filtering based on a threshold, and calculating the proportion of alignments to the correct species.

        Parameters:
        - paf_folder_path: Path to the folder containing PAF files.
        - identity_threshold: Minimum identity percentage to filter alignments.

        Returns:
        - A dictionary where keys are sequence names (from filenames) and values are the proportion
          of alignments to the correct species.
        """
        results = {}

        # Loop through all PAF files in the folder
        for paf_file in glob.glob(os.path.join(paf_folder_path, '*.paf')):
            sequence_name = os.path.basename(paf_file).replace("_alignment.paf", "")
            species_name = sequence_name.replace("_", " ")

            total_alignments = 0
            correct_species_alignments = 0

            with open(paf_file, 'r') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    if len(parts) < 12:
                        continue  # Ensure there are enough fields for analysis

                    # Extract relevant information from the PAF line
                    query_length = int(parts[1])
                    num_matches = int(parts[9])
                    block_length = int(parts[10])

                    # Calculate identity as the percentage of matches over the alignment block length
                    identity = (num_matches / block_length) * 100

                    # Filter alignments based on identity threshold
                    if identity >= identity_threshold:
                        total_alignments += 1
                        ref_species_name = parts[5].split('_')[1]

                        if ref_species_name == species_name:
                            correct_species_alignments += 1

            # Calculate and store the proportion of correct species alignments
            if total_alignments > 0:
                results[sequence_name] = correct_species_alignments / total_alignments
            else:
                results[sequence_name] = 0  # No alignments met the identity threshold

        return results

    def filter_candidate_species_markers(self, target_folder="marker_loci/alignment/minimap2", threshold=0.95):
        """
        Modified to extract region index from sequence names, convert to number, and filter
        candidate_species_markers DataFrame based on species and region index.

        Parameters:
        - target_folder: Path to the target folder containing subdirectories for each species.
        - threshold: Proportion threshold for selecting sequences.
        """
        filtered_sequences = pd.DataFrame()

        # Loop through each subdirectory in the target folder
        for species in os.listdir(target_folder):
            species_dir = os.path.join(target_folder, species)
            if os.path.isdir(species_dir):
                species_results = self.calculate_identity_and_filter(species_dir)

                # Filter sequences based on the threshold and extract region index
                for seq_name, prop in species_results.items():
                    if prop >= threshold:
                        # Extracting the region index from the sequence name
                        _, region_index_str, _ = seq_name.split('_')[-3:]
                        region_index = int(region_index_str)  # Convert index to integer

                        # Select matching rows from candidate_species_markers DataFrame
                        matching_rows = self.candidate_species_markers[
                            (self.candidate_species_markers['species'] == species) &
                            (self.candidate_species_markers['index'] == region_index)
                        ]

                        # Append these rows to the filtered_sequences DataFrame
                        filtered_sequences = pd.concat([filtered_sequences, matching_rows], ignore_index=True)

        return filtered_sequences

    def evaluate_species_identification(self, results, threshold=0.95):
        """
        Evaluates the identification reliability of species based on their alignment proportions.

        Parameters:
        - results: A dictionary where keys are species names and values are proportions of correctly identified alignments.
        - threshold: The proportion threshold above which a species is considered reliably identified.

        Returns:
        - The number of species reliably identified based on the given threshold.
        """
        reliable_identifications = {species: proportion >= threshold for species, proportion in results.items()}
        # Calculate the proportion of species identified reliably
        proportion_reliable_species = sum(reliable_identifications.values()) / len(results) if results else 0
        return proportion_reliable_species

    def cross_species_markers_test_minimap2(self, species_markers, reference_db_base_path, output_base_path):
        """
        Executes minimap2 for marker regions of each species against all other species.

        Parameters:
        - species_markers: DataFrame with columns 'species', 'index', 'consensusSequence', 'start', 'end'.
        - reference_db_base_path: Base path where reference databases for each species are stored.
        - output_base_path: Base path where minimap2 output files will be stored.
        """
        species_list = species_markers['species'].unique()

        for target_species in species_list:
            markers = species_markers[species_markers['species'] == target_species]

            for _, marker in markers.iterrows():
                # Prepare the query sequence
                query_sequence = marker['consensusSequence'][marker['start'] - 1:marker['end']]
                region_index = marker['index']

                with NamedTemporaryFile(mode='w', delete=False) as temp_file:
                    temp_file.write(f">{target_species}_{region_index}\n{query_sequence}\n")
                    temp_file_path = temp_file.name

                for reference_species in species_list:
                    if reference_species != target_species:
                        # Define the reference database path for the current reference_species
                        reference_db = os.path.join(reference_db_base_path, f"{reference_species}.db")

                        # Define the output folder path for the results
                        output_folder = os.path.join(output_base_path, target_species, f"against_{reference_species}")
                        os.makedirs(output_folder, exist_ok=True)

                        # Define the output file path
                        output_file = os.path.join(output_folder, f"{target_species}_{str(region_index).zfill(2)}_vs_{reference_species}.paf")

                        # Construct the minimap2 command
                        minimap2_cmd = [
                            'minimap2',
                            # '-a',
                            reference_db,
                            temp_file_path,
                            '-o',
                            output_file
                        ]

                        # Execute minimap2
                        try:
                            subprocess.run(minimap2_cmd, check=True)
                            print(f"minimap2 alignment completed: {output_file}")
                        except subprocess.CalledProcessError as e:
                            print(f"Error running minimap2: {e}")

                # Cleanup the temporary file
                os.remove(temp_file_path)

    def create_coverage_matrix(self, alignment_threshold=0.5):
        """
        Creates a coverage matrix indicating which marker regions can identify which species,
        based on the analysis of minimap2 results.

        Parameters:
        - alignment_threshold: Threshold for the proportion of correct alignments above which
          a marker is considered reliable for identifying a species.
        """
        # Prepare the coverage matrix structure
        species_list = self.species_markers['species'].unique()
        self.coverage_matrix = pd.DataFrame(columns=['Marker', *species_list])

        # Loop through each species and its markers
        for _, row in self.species_markers.iterrows():
            species = row['species']
            index = row['index']
            marker_id = f"{species}_{str(index).zfill(2)}"  # Construct marker ID as used in PAF filenames

            # Initialize a dictionary to hold the coverage data for this marker
            coverage_data = {'Marker': marker_id}

            # Loop through all species to compare against this marker
            for target_species in species_list:
                if target_species == species:
                    coverage_data[target_species] = None  # Skip self-comparison
                    continue

                # Construct the path to the PAF files for this marker against the target_species
                paf_folder_path = os.path.join(self.output_base_path, species, f"against_{target_species}")

                # Use calculate_identity_and_filter to process the PAF files and get results
                results = self.calculate_identity_and_filter(paf_folder_path)

                # Determine if the marker reliably identifies the target_species based on alignment_threshold
                is_reliable = 1 if results.get(marker_id, 0) >= alignment_threshold else 0
                coverage_data[target_species] = is_reliable

            # Append the coverage data for this marker to the coverage matrix
            self.coverage_matrix = self.coverage_matrix.append(coverage_data, ignore_index=True)

    def save_coverage_matrix_to_csv(self, filename):
        """
        Saves the coverage matrix DataFrame to a CSV file.

        Parameters:
        - filename: The path and name of the file where the coverage matrix will be saved.
        """
        if not self.coverage_matrix.empty:
            self.coverage_matrix.to_csv(filename, index=False)
            print(f"Coverage matrix saved to {filename}.")
        else:
            print("Coverage matrix is empty. No file was created.")

    def select_optimal_markers(self):
        """
        Selects markers that reliably identify the most species, then markers that cover most of the
        remaining species, and so forth, until markers for all species are selected. Species that are determined
        to be uncovered by any available markers are recorded.
        """
        # Convert coverage entries to numeric values for easier processing
        coverage_df = self.coverage_matrix.set_index('Marker').apply(pd.to_numeric, errors='coerce')
        species_list = coverage_df.columns

        # Initialize a set to keep track of covered species
        covered_species = set()

        while covered_species != set(species_list):
            # Sum the coverage for each marker across all species
            marker_coverage = coverage_df.sum(axis=1)

            # Select the marker with the maximum coverage of uncovered species
            best_marker = marker_coverage.idxmax()
            best_coverage = marker_coverage.max()

            # If no marker adds coverage, break the loop
            if best_coverage == 0:
                break

            # Update the list of selected markers
            self.selected_markers.append(best_marker)

            # Update the set of covered species
            newly_covered_species = coverage_df.loc[best_marker][coverage_df.loc[best_marker] == 1].index
            covered_species.update(newly_covered_species)

            # Remove the selected marker to not consider it in the next iteration
            coverage_df.drop(index=best_marker, inplace=True)

        # Identify species without markers
        uncovered_species = set(species_list) - covered_species
        for species in uncovered_species:
            self.species_without_markers[species] = "No available markers"

    def save_optimal_markers_to_csv(self, markers_output_file, uncovered_species_output_file):
        """
        Saves the list of selected markers and the list of species not covered by any markers to separate CSV files.

        Parameters:
        - markers_output_file: Path to the output CSV file for selected markers.
        - uncovered_species_output_file: Path to the output CSV file for species not covered by any markers.
        """
        # Save selected markers
        pd.DataFrame(self.selected_markers, columns=['Selected Markers']).to_csv(markers_output_file, index=False)

        # Save species not covered by any markers
        if self.species_without_markers:
            uncovered_species_df = pd.DataFrame(list(self.species_without_markers.keys()),
                                                columns=['Uncovered Species'])
            uncovered_species_df.to_csv(uncovered_species_output_file, index=False)
        else:
            print("All species are covered by the selected markers.")

    def identify_snps_and_indels(self, alignment_file, format='fasta'):
        """
        Identifies SNPs and indels in a multiple sequence alignment.

        Parameters:
        - alignment_file: Path to the alignment file.
        - format: Format of the alignment file (default 'fasta').

        Returns:
        - A dictionary with 'snps' and 'indels' as keys and lists of positions as values.
        """
        alignment = AlignIO.read(alignment_file, format)
        snps = []
        indels = []

        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            if '-' in column:
                indels.append(i)
            elif len(set(column)) > 1:
                snps.append(i)

        return {'snps': snps, 'indels': indels}
        
    # Function to split multi-FASTA file and select random entries
    def select_random_entries(self, fasta_file, num_entries=1000):
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if len(records) <= num_entries:
            return records
        return random.sample(records, num_entries)
        
    def run_prokka_fast(self, directory="."):
        # Loop over each file in the directory
        for filename in os.listdir(directory):
            if filename.endswith(".fasta"):
                # Extract the species name from the filename
                species_name = filename.replace('chromosome_sequences_', '').replace('.fasta', '')
                
                # Set the output directory for Prokka
                output_dir = f"prokka_{species_name}"
                os.makedirs(output_dir, exist_ok=True)
                
                # Select n random entries or all if the number of entries < n (default n = 1000) from the multi-FASTA file
                fasta_file = os.path.join(directory, filename)
                selected_entries = select_random_entries(fasta_file)
                
                # Write each selected entry to a separate FASTA file and run Prokka
                for i, entry in enumerate(selected_entries, 1):
                    entry_filename = f"{species_name}_{i}.fasta"
                    entry_filepath = os.path.join(output_dir, entry_filename)
                    
                    with open(entry_filepath, "w") as entry_file:
                        SeqIO.write(entry, entry_file, "fasta")
                    
                    # Build the Prokka command
                    prokka_cmd = [
                        'prokka',
                        '--kingdom', 'Bacteria',
                        '--outdir', os.path.join(output_dir, f"entry_{i}"),
                        '--prefix', f"{species_name}_{i}",
                        '--cpus', "0",
                        '--fast',
                        entry_filepath
                    ]
                    
                    # Execute the Prokka command
                    try:
                        print(f"Running Prokka for {species_name}_{i}...")
                        subprocess.run(prokka_cmd, check=True)
                        print(f"Prokka completed successfully for {species_name}_{i}")
                    except subprocess.CalledProcessError as e:
                        print(f"Prokka failed for {species_name}_{i}: {str(e)}")

        # Print completion message
        print("Processing of all files with Prokka is complete.")
        
    def run_prokka(self, directory="."):
        # Loop over each file in the directory
        for filename in os.listdir(directory):
            if filename.endswith(".fasta"):
                # Extract the species name from the filename
                species_name = filename.replace('chromosome_sequences_', '').replace('.fasta', '')
                
                # Set the output directory for Prokka
                output_dir = f"prokka_{species_name}"
                os.makedirs(output_dir, exist_ok=True)
                
                # Select n random entries or all if the number of entries < n (default n = 1000) from the multi-FASTA file
                fasta_file = os.path.join(directory, filename)
                selected_entries = select_random_entries(fasta_file)
                
                # Write each selected entry to a separate FASTA file and run Prokka
                for i, entry in enumerate(selected_entries, 1):
                    entry_filename = f"{species_name}_{i}.fasta"
                    entry_filepath = os.path.join(output_dir, entry_filename)
                    
                    with open(entry_filepath, "w") as entry_file:
                        SeqIO.write(entry, entry_file, "fasta")
                    
                    # Build the Prokka command
                    prokka_cmd = [
                        'prokka',
                        '--kingdom', 'Bacteria',
                        '--outdir', os.path.join(output_dir, f"entry_{i}"),
                        '--prefix', f"{species_name}_{i}",
                        '--cpus', "0",
                        entry_filepath
                    ]
                    
                    # Execute the Prokka command
                    try:
                        print(f"Running Prokka for {species_name}_{i}...")
                        subprocess.run(prokka_cmd, check=True)
                        print(f"Prokka completed successfully for {species_name}_{i}")
                    except subprocess.CalledProcessError as e:
                        print(f"Prokka failed for {species_name}_{i}: {str(e)}")
    
    def gene_data_to_fasta(self, csv_file):
        # Dictionary to store sequences by (species, gene) combination
        sequences = defaultdict(list)
        
        with open(csv_file, 'r') as file:
            reader = csv.reader(file)
            header = next(reader)  # Skip the header row

            for row in reader:
                species_full = row[0]
                dna_sequence = row[5]
                gene_name = row[6]
                
                # Extract the species name (strip the second underscore and everything after it)
                species_name = species_full.split('_')[0] + '_' + species_full.split('_')[1]

                gene_name = gene_name.replace('/', '___')
                
                # Construct the file key
                file_key = f"{species_name}_{gene_name}"
                
                # Prepare the FASTA header and sequence
                fasta_header = f">{species_full}"
                fasta_sequence = dna_sequence
                
                # Add to the dictionary
                sequences[file_key].append((fasta_header, fasta_sequence))
        
        # Write each (species, gene) combination to a separate FASTA file
        for file_key, seq_list in sequences.items():
            with open(f"{file_key}.fasta", 'w') as fasta_file:
                for header, sequence in seq_list:
                    fasta_file.write(f"{header}\n")
                    fasta_file.write(f"{sequence}\n")
                    
    def extract_sequences_and_run_blat(self, species_conservation_data, reference_db):
        """
        Loops through species and their conserved regions, extracts corresponding parts of the consensus sequence,
        and runs BLAT against a reference database for these sequences.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - reference_db: Path to the reference database used by BLAT.
        """
        # Iterate over each species in the conservation data
        for species, (consensus_seq, conserved_regions, _) in species_conservation_data.items():
            print(f"Processing species: {species}")
        
            # Iterate over each conserved region for the current species
            for idx, region_positions in enumerate(conserved_regions, start=1):
                print(f"  Processing conserved region {idx} with positions: {region_positions}")
        
                # Extract the start and end positions of the conserved region
                start_end_tuple = (region_positions[0], region_positions[-1])
                start = start_end_tuple[0]
                end = start_end_tuple[1]
                print(f"    Extracting sequence from positions {start} to {end}")

                # Extract the sequence part corresponding to the current conserved region
                region_seq = consensus_seq[start:end + 1]
                print(f"    Extracted sequence: {region_seq[:50]}...")  # Print the first 50 characters for brevity

                # Create a temporary file to store the region sequence for BLAT input
                temp_seq_file = f"temp_{species}_{start}_{end}.fasta"
                print(f"    Writing sequence to temporary file: {temp_seq_file}")
                with open(temp_seq_file, "w") as f:
                    f.write(f">{species}_{start}_{end}\n{region_seq}\n")

                # Define the output folder dynamically based on species and region
                output_folder_specific = os.path.join("marker_loci/classification/blat")
                print(f"    Ensuring output folder exists: {output_folder_specific}")
                os.makedirs(output_folder_specific, exist_ok=True)

                # Build the BLAT command
                output_file = os.path.join(output_folder_specific, f"{species}_{idx:02d}_{start}_{end}_alignment.psl")
                blat_cmd = ["blat", reference_db, temp_seq_file, output_file]
                print(f"    Running BLAT with command: {' '.join(blat_cmd)}")

                # Run BLAT
                subprocess.run(blat_cmd)

                # Remove the temporary sequence file after BLAT run
                print(f"    Removing temporary file: {temp_seq_file}")
                os.remove(temp_seq_file)

    def extract_sequences_and_run_kraken2(self, species_conservation_data, kraken2_db, num_threads):
        """
        Loops through species and their conserved regions, finds the corresponding alignment files in the current directory,
        extracts the part of the alignment specified by the conserved region, and runs Kraken2 against a reference database
        for these sequences using multiple threads.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - kraken2_db: Path to the Kraken2 database.
        - num_threads: Number of threads to use for Kraken2 classification.
        """
        # Ensure the output folder exists
        output_folder_specific = os.path.join("marker_loci/classification/kraken2")
        print(f"Ensuring output folder exists: {output_folder_specific}")
        os.makedirs(output_folder_specific, exist_ok=True)

        # Iterate over each species in the conservation data
        for species, (consensus_seq, conserved_regions, _) in species_conservation_data.items():
            print(f"Processing species: {species}")

            # Find alignment files for the current species
            alignment_files = glob.glob(f"./{species}*")
            
            for alignment_file in alignment_files:
                print(f"  Found alignment file: {alignment_file}")

                # Process each conserved region
                for idx, region_positions in enumerate(conserved_regions, start=1):
                    print(f"  Processing conserved region {idx} with positions: {region_positions}")
                    
                    # Extract the start and end positions of the conserved region
                    start_end_tuple = (region_positions[0], region_positions[-1])
                    start = start_end_tuple[0]
                    end = start_end_tuple[1]
                    print(f"    Extracting sequence from positions {start} to {end}")

                    # Read the alignment file and extract the relevant part
                    with open(alignment_file, "r") as infile:
                        lines = infile.readlines()
                    
                    extracted_seq_file = f"temp_{species}_{start}_{end}.fasta"
                    with open(extracted_seq_file, "w") as outfile:
                        write_sequence = False
                        for line in lines:
                            if line.startswith(">"):
                                header = line.strip()
                                outfile.write(f"{header}_{start}_{end}\n")
                                write_sequence = True
                            elif write_sequence:
                                seq_part = line.strip()[start:end+1]
                                outfile.write(seq_part + "\n")

                    # Define the output file for Kraken2
                    base_filename = os.path.basename(alignment_file)
                    output_file = os.path.join(output_folder_specific, f"{base_filename}_{idx:02d}_{start}_{end}_classification.txt")
                    
                    # Build the Kraken2 command
                    kraken2_cmd = [
                        "kraken2", 
                        "--db", kraken2_db, 
                        "--output", output_file, 
                        "--report", output_file.replace('.txt', '_report.txt'), 
                        "--threads", str(num_threads), 
                        extracted_seq_file
                    ]
                    print(f"  Running Kraken2 with command: {' '.join(kraken2_cmd)}")

                    # Run Kraken2
                    subprocess.run(kraken2_cmd)

                    # Remove the temporary sequence file after Kraken2 run
                    print(f"  Removing temporary file: {extracted_seq_file}")
                    os.remove(extracted_seq_file)

    def extract_sequences_and_run_pblat(self, species_conservation_data, reference_db, num_cores):
        """
        Loops through species and their conserved regions, extracts corresponding parts of the consensus sequence,
        and runs pBLAT against a reference database for these sequences using the specified number of cores.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - reference_db: Path to the reference database used by pBLAT.
        - num_cores: Number of cores to use for pBLAT.
        """
        # Ensure the output folder exists
        output_folder_specific = os.path.join("marker_loci/classification/pblat")
        print(f"Ensuring output folder exists: {output_folder_specific}")
        os.makedirs(output_folder_specific, exist_ok=True)

        # Iterate over each species in the conservation data
        for species, (consensus_seq, conserved_regions, _) in species_conservation_data.items():
            print(f"Processing species: {species}")

            # Iterate over each conserved region for the current species
            for idx, region_positions in enumerate(conserved_regions, start=1):
                print(f"  Processing conserved region {idx} with positions: {region_positions}")

                # Extract the start and end positions of the conserved region
                start_end_tuple = (region_positions[0], region_positions[-1])
                start = start_end_tuple[0]
                end = start_end_tuple[1]
                print(f"    Extracting sequence from positions {start} to {end}")

                # Extract the sequence part corresponding to the current conserved region
                region_seq = consensus_seq[start:end + 1]
                print(f"    Extracted sequence: {region_seq[:50]}...")  # Print the first 50 characters for brevity

                # Create a temporary file to store the region sequence for pBLAT input
                temp_seq_file = f"temp_{species}_{start}_{end}.fasta"
                print(f"    Writing sequence to temporary file: {temp_seq_file}")
                with open(temp_seq_file, "w") as f:
                    f.write(f">{species}_{start}_{end}\n{region_seq}\n")

                # Define the output file dynamically based on species and region
                output_file = os.path.join(output_folder_specific, f"{species}_{idx:02d}_{start}_{end}_alignment.psl")

                # Build the pBLAT command
                pblat_cmd = ["pblat", "-threads=" + str(num_cores), reference_db, temp_seq_file, output_file]
                print(f"    Running pBLAT with command: {' '.join(pblat_cmd)}")

                # Run pBLAT
                subprocess.run(pblat_cmd)

                # Remove the temporary sequence file after pBLAT run
                print(f"    Removing temporary file: {temp_seq_file}")
                os.remove(temp_seq_file)

    def extract_sequences_and_run_pblat_batch(self, species_conservation_data, reference_db, num_cores, min_identity=98):
        """
        Collects all query sequences and runs pBLAT against a reference database for these sequences using the specified number of cores.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - reference_db: Path to the reference database used by pBLAT.
        - num_cores: Number of cores to use for pBLAT.
        - min_identity:  Minimum sequence identity for pBLAT alignments.
        """
        # Ensure the output folder exists
        output_folder_specific = os.path.join("marker_loci/classification/pblat")
        print(f"Ensuring output folder exists: {output_folder_specific}")
        os.makedirs(output_folder_specific, exist_ok=True)

        # Create a temporary file to store all region sequences for pBLAT input
        combined_seq_file = "combined_temp_sequences.fasta"
        with open(combined_seq_file, "w") as combined_f:
            # Iterate over each species in the conservation data
            for species, (consensus_seq, conserved_regions, _) in species_conservation_data.items():
                print(f"Processing species: {species}")

                # Iterate over each conserved region for the current species
                for idx, region_positions in enumerate(conserved_regions, start=1):
                    print(f"  Processing conserved region {idx} with positions: {region_positions}")

                    # Extract the start and end positions of the conserved region
                    start_end_tuple = (region_positions[0], region_positions[-1])
                    start = start_end_tuple[0]
                    end = start_end_tuple[1]
                    print(f"    Extracting sequence from positions {start} to {end}")

                    # Extract the sequence part corresponding to the current conserved region
                    region_seq = consensus_seq[start:end + 1]
                    print(f"    Extracted sequence: {region_seq[:50]}...")  # Print the first 50 characters for brevity

                    # Write the region sequence to the combined temporary file for pBLAT input
                    combined_f.write(f">{species}_{start}_{end}\n{region_seq}\n")

        # Define the output file for pBLAT
        output_file = os.path.join(output_folder_specific, "combined_alignment.psl")

        # Build the pBLAT command
        pblat_cmd = ["pblat", "-threads=" + str(num_cores), "-minIdentity=" + str(min_identity), reference_db, combined_seq_file, output_file]
        print(f"Running pBLAT with command: {' '.join(pblat_cmd)}")

        # Run pBLAT
        subprocess.run(pblat_cmd)

        # Remove the combined temporary sequence file after pBLAT run
        print(f"Removing temporary file: {combined_seq_file}")
        os.remove(combined_seq_file)

    def extract_sequences_and_run_pblat_batch_partitioned_db(self, species_conservation_data, reference_db_dir, num_cores, min_identity):
        """
        Collects all query sequences and runs pBLAT against each partition of a reference database using the specified number of cores.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - reference_db_dir: Directory containing the partitioned reference database files in .2bit format.
        - num_cores: Number of cores to use for pBLAT.
        - min_identity: Minimum sequence identity for pBLAT alignments.
        """
        # Ensure the output folder exists
        output_folder_specific = os.path.join("marker_loci/classification/pblat")
        print(f"Ensuring output folder exists: {output_folder_specific}")
        os.makedirs(output_folder_specific, exist_ok=True)

        # Create a temporary file to store all region sequences for pBLAT input
        combined_seq_file = "combined_temp_sequences.fasta"
        with open(combined_seq_file, "w") as combined_f:
            # Iterate over each species in the conservation data
            for species, (consensus_seq, conserved_regions, _) in species_conservation_data.items():
                print(f"Processing species: {species}")

                # Iterate over each conserved region for the current species
                for idx, region_positions in enumerate(conserved_regions, start=1):
                    print(f"  Processing conserved region {idx} with positions: {region_positions}")

                    # Extract the start and end positions of the conserved region
                    start_end_tuple = (region_positions[0], region_positions[-1])
                    start = start_end_tuple[0]
                    end = start_end_tuple[1]
                    print(f"    Extracting sequence from positions {start} to {end}")

                    # Extract the sequence part corresponding to the current conserved region
                    region_seq = consensus_seq[start:end + 1]
                    print(f"    Extracted sequence: {region_seq[:50]}...")  # Print the first 50 characters for brevity

                    # Write the region sequence to the combined temporary file for pBLAT input
                    combined_f.write(f">{species}_{start}_{end}\n{region_seq}\n")

        # Iterate over each partitioned reference database file
        for ref_db_part in os.listdir(reference_db_dir):
            if ref_db_part.endswith(".2bit"):
                ref_db_path = os.path.join(reference_db_dir, ref_db_part)

                # Define the output file for pBLAT
                output_file = os.path.join(output_folder_specific, f"combined_alignment_{ref_db_part}.psl")

                # Build the pBLAT command
                pblat_cmd = ["pblat", "-threads=" + str(num_cores), "-minIdentity=" + str(min_identity), ref_db_path, combined_seq_file, output_file]
                print(f"Running pBLat for {ref_db_part} with command: {' '.join(pblat_cmd)}")

                # Run pBLAT
                subprocess.run(pblat_cmd)

        # Remove the combined temporary sequence file after all pBLAT runs
        print(f"Removing temporary file: {combined_seq_file}")
        os.remove(combined_seq_file)
        
    def extract_sequences_and_run_pblat_batch_partitioned(self, species_conservation_data, reference_db, num_cores):
        """
        Collects all query sequences and runs pBLAT against a reference database for these sequences using the specified number of cores.
        The sequences are split into four parts and processed in four separate pBLAT runs.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - reference_db: Path to the reference database used by pBLAT.
        - num_cores: Number of cores to use for pBLAT.
        """
        # Ensure the output folder exists
        output_folder_specific = os.path.join("marker_loci/classification/pblat")
        print(f"Ensuring output folder exists: {output_folder_specific}")
        os.makedirs(output_folder_specific, exist_ok=True)

        # Collect all sequences into a list
        sequences = []
        for species, (consensus_seq, conserved_regions, _) in species_conservation_data.items():
            print(f"Processing species: {species}")

            for idx, region_positions in enumerate(conserved_regions, start=1):
                print(f"  Processing conserved region {idx} with positions: {region_positions}")

                start_end_tuple = (region_positions[0], region_positions[-1])
                start = start_end_tuple[0]
                end = start_end_tuple[1]
                print(f"    Extracting sequence from positions {start} to {end}")

                region_seq = consensus_seq[start:end + 1]
                print(f"    Extracted sequence: {region_seq[:50]}...")  # Print the first 50 characters for brevity

                sequences.append(f">{species}_{start}_{end}\n{region_seq}\n")

        # Split sequences into four parts
        num_sequences = len(sequences)
        part_size = ceil(num_sequences / 4)

        for i in range(4):
            part_sequences = sequences[i * part_size:(i + 1) * part_size]
            part_file = f"combined_temp_sequences_part_{i+1}.fasta"
            
            with open(part_file, "w") as part_f:
                part_f.writelines(part_sequences)

            # Define the output file for pBLAT
            output_file = os.path.join(output_folder_specific, f"combined_alignment_part_{i+1}.psl")

            # Build the pBLAT command
            pblat_cmd = ["pblat", "-threads=" + str(num_cores), reference_db, part_file, output_file]
            print(f"Running pBLAT for part {i+1} with command: {' '.join(pblat_cmd)}")

            # Run pBLAT
            subprocess.run(pblat_cmd)

            # Remove the part temporary sequence file after pBLAT run
            print(f"Removing temporary file: {part_file}")
            # os.remove(part_file)

    def process_conservation_data(self, species_conservation_data, folder_path="."):
        """
        Processes the species conservation data by finding consensus sequences for related regions and appending new entries
        to the output dictionary.

        Parameters:
        - species_conservation_data: Dictionary with species names as keys and 3-tuples (consensus sequence,
          conserved region positions, and dominant letters) as values.
        - folder_path: Path to the directory containing .aln files. Default is the current directory.

        Returns:
        - Dictionary with updated entries including consensus sequences for related regions.
        """
        output_data = {}

        for key, (consensus_sequence, conserved_regions, dominant_letters) in species_conservation_data.items():
            # Extract region_name and current_species_name from the key
            parts = key.split("_")
            region_name = "_".join(parts[2:])
            current_species_name = "_".join(parts[:2])

            for file_name in os.listdir(folder_path):
                if file_name.endswith(".aln") and region_name in file_name and current_species_name not in file_name:
                    alignment_file = os.path.join(folder_path, file_name)

                    # Find the consensus sequence for the alignment file
                    new_consensus_sequence = self.get_consensus_sequence_from_alignment(alignment_file)

                    # Build the new key
                    new_species_name = "_".join(file_name.split("_")[:2])
                    new_key = f"{current_species_name}_{region_name}_{new_species_name}"

                    # Append the new entry to the output dictionary
                    output_data[new_key] = (new_consensus_sequence, conserved_regions, dominant_letters)

        return output_data

    def calculate_frequencies(self, alignment):
        num_sequences = len(alignment)
        num_positions = alignment.get_alignment_length()
        frequencies = []

        for pos in range(num_positions):
            nucleotide_counts = defaultdict(float)
            total_count = 0  # To count only valid nucleotides

            for record in alignment:
                nucleotide = record.seq[pos]
                if nucleotide not in ['-', 'N']:
                    nucleotide_counts[nucleotide] += 1
                    total_count += 1

            # Normalize frequencies
            for nucleotide in nucleotide_counts:
                nucleotide_counts[nucleotide] /= total_count

            frequencies.append(nucleotide_counts)

        return frequencies

    def renyi_entropy(self, frequencies, alpha=2):
        if frequencies is None:
            return None  # Skip positions with only gaps and Ns
        if alpha == 1:
            return -sum(p * np.log2(p) for p in frequencies.values() if p > 0)
        else:
            sum_p_alpha = sum(p ** alpha for p in frequencies.values())
            return (1 / (1 - alpha)) * np.log2(sum_p_alpha)

    def find_species_conserved_regions(self, directory="separate_species"):
        data = []
        
        for filename in os.listdir(directory):
            if filename.endswith(".aln"):
                filepath = os.path.join(directory, filename)
                alignment = AlignIO.read(filepath, "fasta")
                
                # Parse the filename to extract species_name and region
                parts = filename.split("_")
                species_name = "_".join(parts[:2])
                region = "_".join(parts[2:]).rsplit("_", 1)[0]
                
                # Find conserved regions
                conserved_regions = self.find_low_entropy_regions(alignment, threshold=0.3, alpha=1.5)
                
                # Append the result to the data list
                data.append({
                    "region": region,
                    "species": species_name,                    
                    "conserved_regions": conserved_regions
                })
        
        # Create a DataFrame from the data list
        df = pd.DataFrame(data, columns=["region", "species", "conserved_regions"])
        return df

    def find_high_entropy_regions(self, alignment, threshold=1.4, alpha=2, min_length=3):
        frequencies = self.calculate_frequencies(alignment)
        num_positions = alignment.get_alignment_length()
        
        entropy_values = [self.renyi_entropy(freq, alpha) for freq in frequencies]
        
        regions = []
        in_region = False
        start = None
        region_entropies = []

        for i in range(num_positions):
            if entropy_values[i] is not None:
                if entropy_values[i] > threshold:
                    if not in_region:
                        start = i
                        in_region = True
                        region_entropies = [entropy_values[i]]
                    else:
                        region_entropies.append(entropy_values[i])
                    avg_entropy = sum(region_entropies) / len(region_entropies)
                    if avg_entropy <= threshold:
                        in_region = False
                        if i - start >= min_length:
                            regions.append((start, i - 1))
                else:
                    if in_region:
                        avg_entropy = sum(region_entropies) / len(region_entropies)
                        if avg_entropy > threshold:
                            in_region = False
                            if i - start >= min_length:
                                regions.append((start, i - 1))
            else:
                if in_region:
                    avg_entropy = sum(region_entropies) / len(region_entropies)
                    if avg_entropy > threshold:
                        in_region = False
                        if i - start >= min_length:
                            regions.append((start, i - 1))

        # Check if the last region extends to the end
        if in_region and num_positions - start >= min_length:
            regions.append((start, num_positions - 1))

        return regions

    def find_low_entropy_regions(self, alignment, threshold=0.2, alpha=2, min_length=3):
        frequencies = self.calculate_frequencies(alignment)
        num_positions = alignment.get_alignment_length()
        
        entropy_values = [self.renyi_entropy(freq, alpha) for freq in frequencies]
        
        regions = []
        in_region = False
        start = None
        region_entropies = []

        for i in range(num_positions):
            if entropy_values[i] is not None:
                if entropy_values[i] <= threshold:
                    if not in_region:
                        start = i
                        in_region = True
                        region_entropies = [entropy_values[i]]
                    else:
                        region_entropies.append(entropy_values[i])
                    avg_entropy = sum(region_entropies) / len(region_entropies)
                    if avg_entropy > threshold:
                        in_region = False
                        if i - start >= min_length:
                            regions.append((start, i - 1))
                else:
                    if in_region:
                        avg_entropy = sum(region_entropies) / len(region_entropies)
                        if avg_entropy <= threshold:
                            in_region = False
                            if i - start >= min_length:
                                regions.append((start, i - 1))
            else:
                if in_region:
                    avg_entropy = sum(region_entropies) / len(region_entropies)
                    if avg_entropy > threshold:
                        in_region = False
                        if i - start >= min_length:
                            regions.append((start, i - 1))

        # Check if the last region extends to the end
        if in_region and num_positions - start >= min_length:
            regions.append((start, num_positions - 1))

        return regions

    def find_cross_species_conserved_and_variable_regions(self, filepath, min_length=3):
        data = []
        
        # Read the alignment file
        alignment = AlignIO.read(filepath, "fasta")
        
        # Parse the filename to extract region
        filename = os.path.basename(filepath)
        
        # Step 1: Remove the file extension
        base_filename = filename.rsplit(".", 1)[0]
        # base_filename = 'genus_species_some_region_with_underscores'

        # Step 2: Split the base filename by underscores
        parts = base_filename.split("_")
        # parts = ['genus', 'species', 'some', 'region', 'with', 'underscores']

        # Step 3: The region starts after the genus and species names
        region = "_".join(parts[2:])
        # region = 'some_region_with_underscores'
        
        # Extract species names
        species_names = set()
        for record in alignment:
            species_name = "_".join(record.id.split("_")[:2])
            species_names.add(species_name)
        
        # Find variable regions with high entropy
        high_entropy_regions = self.find_high_entropy_regions(alignment, threshold=1.2, alpha=2, min_length=min_length)
        
        # Find conserved regions with low entropy
        low_entropy_regions = self.find_low_entropy_regions(alignment, threshold=0.4, alpha=2, min_length=min_length)
        
        # Append the result to the data list
        data.append({
            "region": region,
            "species": species_names,
            "variable_regions": high_entropy_regions,
            "conserved_regions": low_entropy_regions
        })
        
        # Create a DataFrame from the data list
        df = pd.DataFrame(data, columns=["region", "species", "variable_regions", "conserved_regions"])
        return df
        
    def process_cross_species_alignments(self, directory="."):
        all_data = []
        
        for filename in os.listdir(directory):
            if filename.endswith(".fas"):
                filepath = os.path.join(directory, filename)
                alignment = AlignIO.read(filepath, "fasta")

                # Parse the filename to extract region
                filename = os.path.basename(filepath)                
                
                print(f"processing {filepath}")
                
                # Extract species names
                species_names = set("_".join(record.id.split("_")[:2]) for record in alignment)

                # Step 1: Remove the file extension
                region = filename.rsplit(".", 1)[0]
                # base_filename = 'genus_species_some_region_with_underscores'
                
                # Step 2: Split the base filename by underscores
                # parts = base_filename.split("_")
                # parts = ['genus', 'species', 'some', 'region', 'with', 'underscores']

                # Step 3: The region starts after the genus and species names
                # region = "_".join(parts[2:])
                # region = 'some_region_with_underscores'
                
                # Generate all combinations of species names of different lengths
                for r in range(2, len(species_names) + 1):
                    for combination in itertools.combinations(species_names, r):
                        # print(combination)

                        # Create a temporary file for this combination
                        temp_alignment = [record for record in alignment if "_".join(record.id.split("_")[:2]) in combination]
                        temp_filename = f"tmp/species_combination_{region}.aln"
                        with open(temp_filename, "w") as temp_file:
                            write(temp_alignment, temp_file, "fasta")
                        
                        # Run the analysis function
                        df = self.find_cross_species_conserved_and_variable_regions(temp_filename)
                        
                        # Append the result to the all_data list
                        all_data.append(df)
        
        # Combine all the data into a single DataFrame
        result_df = pd.concat(all_data, ignore_index=True)
        return result_df
        
    def merge_species_and_cross_species_conservation_data(self, df1, df2):
        # Create a new column in df1 to store the list of lists of conserved regions for each species
        df1["species_conserved_regions"] = None

        # Iterate over the rows of df1
        for i, row in df1.iterrows():
            region = row["region"]
            species_set = row["species"]
            
            species_conserved_regions = []
            
            # For each species in the set, find the corresponding row in df2 and get the conserved regions
            for species in species_set:
                matching_rows = df2[(df2["species"] == species) & (df2["region"] == region)]
                
                # There should be only one matching row for each species and region
                for _, match in matching_rows.iterrows():
                    species_conserved_regions.append(match["conserved_regions"])
            
            # Store the list of lists in the new column
            df1.at[i, "species_conserved_regions"] = species_conserved_regions
        
        return df1

    # Function to convert string to list of tuples
    def convert_string_to_object(self, string):
        try:
            return ast.literal_eval(string)
        except (ValueError, SyntaxError) as e:
            print(f"Error converting string to object: {string}")
            return []

    # Function to merge adjacent regions
    def merge_adjacent_regions(self, regions, gap=5):
        if not regions:
            return []
        
        # Sort regions by start position
        regions = sorted(regions, key=lambda x: x[0])
        combined_regions = []
        current_region = regions[0]

        for start, end in regions[1:]:
            if start - current_region[1] <= gap:
                # Merge regions
                current_region = (current_region[0], max(current_region[1], end))
            else:
                # Add current region to combined regions and start a new region
                combined_regions.append(current_region)
                current_region = (start, end)
        
        # Add the last region
        combined_regions.append(current_region)
        return combined_regions

    # Function to process the DataFrame
    def combine_variable_regions(self, df, gap=5):
        df['variable_regions_combined'] = df['variable_regions'].apply(lambda x: self.merge_adjacent_regions(x, gap))
        return df

    # Function to process the DataFrame
    def combine_conserved_regions(self, df, gap=5):
        df['conserved_regions_combined'] = df['conserved_regions'].apply(lambda x: self.merge_adjacent_regions(x, gap))
        return df

    # Function to combine species conserved regions
    def combine_species_conserved_regions(self, df, gap=5):
        def combine_species_regions(species_regions):
            return [self.merge_adjacent_regions(regions, gap) for regions in species_regions]

        df['species_conserved_regions_combined'] = df['species_conserved_regions'].apply(combine_species_regions)
        return df

    # Function to filter regions based on minimum length
    def filter_regions_by_length(self, regions, min_length):
        return [region for region in regions if region[1] - region[0] + 1 >= min_length]

    # Function to process the DataFrame
    def filter_candidate_marker_regions(self, df, min_length):
        df = df.copy()  # Ensure we are working with a copy of the DataFrame
        df.loc[:, 'candidate_marker_regions_filtered'] = df['candidate_marker_regions'].apply(lambda x: self.filter_regions_by_length(x, min_length))
        return df
        
    def intersect_two_regions(self, region1, region2):
        """Find the intersection of two regions."""
        try:
            start1, end1 = region1
            start2, end2 = region2
        except ValueError as e:
            print(f"Error unpacking regions: {region1}, {region2}")
            raise e
        if end1 < start2 or end2 < start1:
            return None  # No overlap
        return (max(start1, start2), min(end1, end2))

    def find_intersection(self, regions):
        """Find the intersection of multiple lists of regions."""
        if not regions:
            return []

        # Initialize with the first list of regions
        intersection = regions[0]
        # print(f"regions = {regions}")
        # print(f"{intersection}")
        for other_regions in regions[1:]:
            new_intersection = []
            for region1 in intersection:
                for region2 in other_regions:
                    # print(f"region1 = {region1}")
                    # print(f"region2 = {region2}")
                    intersected_region = self.intersect_two_regions(region1, region2)
                    if intersected_region:
                        new_intersection.append(intersected_region)
            intersection = new_intersection

            # If at any point the intersection is empty, we can return immediately
            if not intersection:
                return []

        return sorted(intersection)

    def calculate_intersections(self, df):
        df['species_conserved_regions_intersection'] = df['species_conserved_regions'].apply(self.find_intersection)
        return df

    # Function to calculate intersections between species_conserved_regions_intersection and variable_regions
    def intersect_region_lists(self, list1, list2):
        """Find the intersection of two lists of regions."""
        intersections = []
        for region1 in list1:
            for region2 in list2:
                intersected_region = self.intersect_two_regions(region1, region2)
                if intersected_region:
                    intersections.append(intersected_region)
        return sorted(intersections)

    def calculate_candidate_marker_regions(self, df):
        df['candidate_marker_regions'] = df.apply(
            lambda row: self.intersect_region_lists(row['species_conserved_regions_intersection'], row['variable_regions']), axis=1)
        return df
        
    def is_subregion(self, region, conserved_regions, window_size):
        """Check if a subregion of a given window size is within any of the conserved regions."""
        start, end = region
        for conserved_start, conserved_end in conserved_regions:
            if conserved_start <= start <= conserved_end and conserved_start <= start + window_size - 1 <= conserved_end:
                return True
            if conserved_start <= end - window_size + 1 <= conserved_end and conserved_start <= end <= conserved_end:
                return True
        return False
        
    def find_adjusted_marker_regions(self, row, max_bases=30, initial_window_size=20, min_window_size=18):
        candidate_regions = row['candidate_marker_regions']
        conserved_regions = row['conserved_regions']
        region = row['region']
        species = row['species']

        output_data = []

        for start, end in candidate_regions:
            shift_left = 0
            shift_right = 0
            left_window_size = initial_window_size
            right_window_size = initial_window_size

            for window_size in range(initial_window_size, min_window_size - 1, -1):
                for shift in range(max_bases + 1):
                    # Test left end
                    new_start = start - shift
                    if self.is_subregion((new_start, new_start + window_size - 1), conserved_regions, window_size):
                        shift_left = shift
                        left_window_size = window_size
                        break
                else:
                    continue
                break

            for window_size in range(initial_window_size, min_window_size - 1, -1):
                for shift in range(max_bases + 1):
                    # Test right end
                    new_end = end + shift
                    if self.is_subregion((new_end - window_size + 1, new_end), conserved_regions, window_size):
                        shift_right = shift
                        right_window_size = window_size
                        break
                else:
                    continue
                break

            if shift_left or shift_right:
                sequence_for_primer_design = (start - shift_left, end + shift_right)
                output_data.append({
                    'region': region,
                    'species': species,
                    'candidate_marker_region': (start, end),
                    'sequence_for_primer_design': sequence_for_primer_design,
                    'shift_left': shift_left,
                    'shift_right': shift_right,
                    'primer_binding_site_length_left': left_window_size,
                    'primer_binding_site_length_right': right_window_size
                })

        return output_data
        
    def process_marker_regions(self, df, max_bases=30, initial_window_size=20, min_window_size=18):
        output_rows = []

        for index, row in df.iterrows():
            adjusted_regions = self.find_adjusted_marker_regions(row, max_bases, initial_window_size, min_window_size)
            output_rows.extend(adjusted_regions)

        return pd.DataFrame(output_rows, columns=[
            'region', 'species', 'candidate_marker_region', 'sequence_for_primer_design', 
            'shift_left', 'shift_right', 'primer_binding_site_length_left', 
            'primer_binding_site_length_right'
        ])
        
    def slice_alignment(self, alignment_file, start, end, format="fasta"):
        # Read the alignment file
        alignment = AlignIO.read(alignment_file, format)
        
        # Slice the alignment
        sliced_records = []
        for record in alignment:
            sliced_seq = record.seq[start:end]
            sliced_record = SeqRecord(sliced_seq, id=record.id, description=record.description)
            sliced_records.append(sliced_record)
        
        # Create a new alignment from the sliced records
        sliced_alignment = MultipleSeqAlignment(sliced_records)
        
        return sliced_alignment

    def identify_markers(self, ):
        self.species_markers = filter_candidate_species_markers()
        save_optimal_markers_to_csv()

    def design_primers(self, input_file, output_file, primer_design_algorithm: PrimerDesignAlgorithm, primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params, specificity_check_database):
        pass

class MarkerLociGroupDifferentiationStrategy(Strategy):
    @staticmethod
    def sequence_to_kmers(sequence, k):
        return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

    @staticmethod
    def integer_encode_kmer(kmer):
        encoding_map = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
        return [encoding_map[nucleotide] for nucleotide in kmer]

    @staticmethod
    def positional_encode_sequence(sequence, k):
        encoded_kmers = []
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            encoded_kmer = MarkerLociIdentificationStrategy.integer_encode_kmer(kmer)
            # Append positional information
            positional_encoded_kmer = [i] + encoded_kmer  # Prepending the position index
            encoded_kmers.append(positional_encoded_kmer)
        return encoded_kmers

    class _FeatureExtractor(ABC):
        @abstractmethod
        def extract_features(self, sequence):
            pass

    class _KmerFeatureExtractor(_FeatureExtractor):
        def extract_features(self, sequence, k=3):
            return MarkerLociIdentificationStrategy.sequence_to_kmers(sequence, k)

    class _Word2VecFeatureExtractor(_FeatureExtractor):

        def extract_features(self, sequence):
            kmer_extractor = MarkerLociIdentificationStrategy._KmerFeatureExtractor()

            kmers = kmer_extractor.extract_features(sequence)

            model = Word2Vec(sentences=kmers, vector_size=100, window=5, min_count=1, workers=4)

            processed_embeddings = self.process_embeddings(model, kmers)

            return processed_embeddings

        def process_embeddings(self, model, kmers):
            # todo: possibly implement a more advanced logic to process the embeddings

            return [model.wv[kmer].mean(axis=0) for kmer in kmers]

    class _GenomicTransformer(nn.Module):
        def __init__(self, embedding_size, nhead, num_layers, num_classes):
            super().__init__()
            self.transformer = nn.Transformer(d_model=embedding_size, nhead=nhead, num_layers=num_layers)
            self.fc = nn.Linear(embedding_size, num_classes)  # Final classification layer

        def forward(self, x):
            # x: [seq_len, batch_size, embedding_size]
            x = self.transformer(x, x)

            x = x.mean(dim=0)
            x = self.fc(x)
            return x

    class _GenomicTransformer(nn.Module):
        def __init__(self, embedding_size, nhead, num_layers, num_classes):
            super().__init__()
            config = BertConfig(
                hidden_size=embedding_size,
                num_attention_heads=nhead,
                num_hidden_layers=num_layers,
            )
            self.transformer = BertModel(config)
            self.fc = nn.Linear(embedding_size, num_classes)  # Final classification layer

        def forward(self, x):
            # x: [batch_size, seq_len, embedding_size]
    
            # Adjust x to match BERT input dimensions
            x = x.permute(1, 0, 2)  # BERT expects [seq_len, batch_size, embedding_size]

            # Get transformer outputs (last hidden state and attentions)
            outputs = self.transformer(x, output_attentions=True)
            hidden_state = outputs.last_hidden_state
            attentions = outputs.attentions

            # Aggregate over sequence and pass through the classification layer
            aggregated_output = hidden_state.mean(dim=1)
            logits = self.fc(aggregated_output)

            return logits, attentions

    class _TrainGenomicTransformer(Trainable):
        def _setup(self, config):
            # Model and optimizer initialization
            self.model = MarkerLociIdentificationStrategy._GenomicTransformer(config["embedding_size"], config["nhead"], config["num_layers"],
                                            config["num_classes"])
            self.optimizer = torch.optim.Adam(self.model.parameters(), lr=config["lr"])
            self.criterion = torch.nn.CrossEntropyLoss()

            self.device = "cuda" if torch.cuda.is_available() else "cpu"
            self.model.to(self.device)

            # Load and preprocess data
            file_path = config["file_path"]
            train_df, val_df = self._load_and_label_data(file_path)

            # Initialize the label encoder
            label_encoder = LabelEncoder()

            # Create separate DataLoaders for training and validation
            self.train_loader = self.create_data_loader(train_df, config["batch_size"])
            self.val_loader = self.create_data_loader(val_df, config["batch_size"])

        def _create_data_loader(self, df, batch_size):
            # Convert string labels to integers using the instance's label_encoder
            labels = self.label_encoder.fit_transform(df['label'])

            # Convert embeddings list of lists to tensors and pad sequences
            features_list = [torch.tensor(embeddings) for embeddings in df['extracted_features']]
            features_padded = pad_sequence(features_list, batch_first=True, padding_value=0)
            labels = torch.tensor(labels)

            # Create and return a DataLoader
            dataset = TensorDataset(features_padded, labels)
            return DataLoader(dataset, batch_size=batch_size, shuffle=True)

        def _train(self):
            num_epochs = self.config.get("num_epochs", 10)

            for epoch in range(num_epochs):
                # Training phase
                self.model.train()
                total_train_loss, total_train_accuracy = 0, 0
                for batch in self.train_loader:
                    inputs = batch[0].to(self.device)
                    labels = batch[1].to(self.device)
                    self.optimizer.zero_grad()
                    outputs = self.model(inputs)
                    loss = self.criterion(outputs, labels)
                    loss.backward()
                    self.optimizer.step()
                    total_train_loss += loss.item()

                    # Calculate accuracy
                    predicted = torch.argmax(outputs, dim=1)
                    total_train_accuracy += accuracy_score(labels.cpu(), predicted.cpu())

                avg_train_loss = total_train_loss / len(self.train_loader)
                avg_train_accuracy = total_train_accuracy / len(self.train_loader)

                # Validation phase with attention weights collection
                self.model.eval()
                total_val_loss, total_val_accuracy = 0, 0
                with torch.no_grad():
                    for batch in self.val_loader:
                        inputs = batch[0].to(self.device)
                        labels = batch[1].to(self.device)
                        outputs, attentions = self.model(inputs)  # Get outputs and attention weights
                        loss = self.criterion(outputs, labels)
                        total_val_loss += loss.item()

                        # Calculate accuracy
                        predicted = torch.argmax(outputs, dim=1)
                        total_val_accuracy += accuracy_score(labels.cpu(), predicted.cpu())

                        # Collect attention weights for analysis
                        collected_attention_weights.append(attentions)

                avg_val_loss = total_val_loss / len(self.val_loader)
                avg_val_accuracy = total_val_accuracy / len(self.val_loader)

                # Report metrics for each epoch
                tune.report(train_loss=avg_train_loss, train_accuracy=avg_train_accuracy,
                            val_loss=avg_val_loss, val_accuracy=avg_val_accuracy)

            return {"train_loss": avg_train_loss, "train_accuracy": avg_train_accuracy,
                    "val_loss": avg_val_loss, "val_accuracy": avg_val_accuracy}
        
        def analyze_attention_weights(self, attention_weights):
            pass

    def __init__(self, input_file, feature_extractor):
        self._input_file = input_file
        self._feature_extractor = feature_extractor

        # todo: elaborate, connect to CLI parameters
        self._ray_config = {
            "file_path": self._input_file,
            "batch_size": 128,
            "num_epochs": 10,
            "embedding_size": 100,
            "nhead": hp.choice("nhead", [2, 4, 6, 8]),
            "num_layers": 2,
            "num_classes": 2,
            "lr": hp.uniform("lr", 0.0001, 0.01),
        }

    def _load_and_label_data(self, file_path):
        """Loads and labels data from a single FASTA file based on sequence IDs."""
        data = []
        sequence_dict = parse_fasta_to_dict(file_path)
        for seq_id, seq in sequence_dict.items():

            if seq_id.startswith("pathogenic_"):
                pathogenicity = "pathogenic"
            elif seq_id.startswith("non-pathogenic_"):
                pathogenicity = "non-pathogenic"
            else:
                continue

            extracted_features = self._feature_extractor.extract_features(seq)
            data.append({'sequence_id': seq_id, 'extracted_features': extracted_features, 'label': pathogenicity})

            data_df = pd.DataFrame(data)

            # Split data into training and validation sets
            train_df, val_df = train_test_split(data_df, test_size=0.2)

        return train_df, val_df

    def identify_markers(self):
        print("Identifying marker loci...")

        hyperopt_search = HyperOptSearch(space=self._ray_config, metric="accuracy", mode="max")

        # todo: elaborate, connect to CLI parameters
        analysis = tune.run(MarkerLociIdentificationStrategy._TrainGenomicTransformer,
                            config=self._ray_config,
                            num_samples=10,
                            search_alg=hyperopt_search,
                            resources_per_trial={"cpu": 1, "gpu": 1},
                            metric="accuracy",
                            mode="max")

        # Find the best hyperparameters
        best_config = analysis.get_best_config(metric="accuracy", mode="max")

    def design_primers(self, input_file, output_file, primer_design_algorithm: PrimerDesignAlgorithm, primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params, specificity_check_database):
        print("Designing primers for detected marker loci...")


# Context class to use the strategies
class StrategyContext:
    def __init__(self, strategy: Strategy, input_file, output_file, database):
        self._strategy = strategy
        self._input_file = input_file
        self._output_file = output_file
        self._database = database

        if isinstance(self._strategy, TargetedAmpliconSequencingStrategy):
        # self._sequences = parse_fasta_to_dict(input_file)
            self._amplicon_sequences = parse_fasta_by_amplicons(self._input_file)

    def design_primers(self, primer_design_algorithm: PrimerDesignAlgorithm, primer_design_params, primer_summarizing_algorithm: PrimerSummarizerAlgorithm, primer_summarizing_params, specificity_check_algorithm: SpecificityCheckAlgorithm, primer_specificity_params):
        if isinstance(self._strategy, MarkerLociIdentificationStrategy):
            self._strategy.identify_markers()

            # todo:  loop through marker loci, design primers for each of them and output them to a file

        elif isinstance(self._strategy, TargetedAmpliconSequencingStrategy):
            # designed_primers = {}

            for barcode, sequences in self._amplicon_sequences.items():
                primers = self._strategy.design_primers(sequences, self._output_file, primer_design_algorithm, primer_design_params, primer_summarizing_algorithm, primer_summarizing_params, specificity_check_algorithm, primer_specificity_params, self._database)

                # designed_primers[barcode] = primers

                output_primers_to_csv(primers[0], primers[1], self._output_file+barcode+'left_primers.csv', self._output_file+barcode+'right_primers.csv')


if __name__ == "__main__":
    ray.init()

    parser = argparse.ArgumentParser(description="")

    # Add arguments related to marker loci identification to the parser
    parser.add_argument("-a", "--alignment-tool", choices=['clustalw', 'muscle', 'mafft', 'prank'],
                        help="Alignment tool to use. Possible values: clustalw, muscle, mafft, prank.")

    parser.add_argument("-r", "--reference-db", required=True,
                        help="Path to the reference database for minimap2.")

    parser.add_argument("-I", "--identity-threshold", type=int, choices=range(0, 101),
                        help="Identity threshold for minimap2 (0-100).")

    parser.add_argument("-P", "--proportion-threshold", type=int, choices=range(0, 101),
                        help="Threshold for proportion of correctly identified sequences (0-100).")

    parser.add_argument("-l", "--min-length", type=int,
                        help="Minimum length of conserved regions.")

    parser.add_argument("-g", "--algorithm", choices=['consensus_sequence', 'shannon_entropy', 'quasi_alignment'],
                        help="Algorithm for finding conserved regions. Possible values: consensus_sequence, shannon_entropy, quasi_alignment.")

    # Add arguments related to primer design to the parser
    parser.add_argument('-i', '--input_file', required=True, help='The input file to process.')
    parser.add_argument('-o', '--output_file', required=True, default='output.txt', help='The output file to write to.')
    parser.add_argument('-s', '--strategy', required=True, default='amplicon', choices=['amplicon', 'marker'], help='Choose a general strategy - primer design only for the WGS dataset (amplicon) or identification of marker loci in the WGS dataset followed by the primer design for the identified marker loci (marker).')
    parser.add_argument('-a', '--primer_design_algorithm', required=True, default='primer3', choices=['primer3', 'custom'], help='Choose a primer design algorithm.')
    parser.add_argument('-p', '--primer3_parameters', type=str, help='Path to the Primer3 parameters config file.')
    parser.add_argument('-S', '--primer_summarizing_algorithm', required=True, choices=['frequency', 'consensus'], help='An algorithm for summarizing primers.')
    parser.add_argument('-c', '--specificity_check_algorithm', required=True, choices=['blast'], help='An algorithm for primer specificity checking.')
    parser.add_argument('-d', '--database_for_specificity_check', required=True, help = 'A database for checking the specificity of primers.')
    parser.add_argument('-n', '--n_most_frequent', help='A number of the most frequently occurring primers to further work with.')

    args = parser.parse_args()

    if args.strategy == 'amplicon':
        strategy = TargetedAmpliconSequencingStrategy()
    elif args.strategy == 'marker':
        feature_extractor = MarkerLociIdentificationStrategy._Word2VecFeatureExtractor()  # todo: select feature extraction algorithm based on the CLI parameters
        strategy = MarkerLociIdentificationStrategy(args.input_file, feature_extractor)
    else:
        raise ValueError("Unknown strategy specified.")

    context = StrategyContext(strategy, args.input_file, args.output_file, args.database_for_specificity_check)

    if args.strategy == 'marker':
        strategy.identify_markers()

    # Initialize configparser
    config = configparser.ConfigParser()

    if args.primer_design_algorithm == 'primer3':
        if args.primer3_parameters is None:
            raise ValueError("The primer3 config file was not specified.")

        print(f"Primer3 parameters will be loaded from: {args.primer3_parameters}")

        # Read the configuration file
        config.read(args.primer3_parameters)

        primer_design_algorithm = Primer3Algorithm()

        primer_design_params = {param.upper(): config.get('Primer3Parameters', param) for param in config.options('Primer3Parameters')}
        primer_design_params_converted = {key: convert_param(value) for key, value in primer_design_params.items()}
    elif args.primer_design_algorithm == 'custom':
        primer_design_algorithm = CustomAlgorithm()
    else:
        raise ValueError("Unknown primer design algorithm specified.")

    if args.primer_summarizing_algorithm == 'frequency':
        primer_summarizing_algorithm = FrequencyBasedSummarizer()
        primer_summarizing_params = {'n_most_frequent': args.n_most_frequent}
    elif args.primer_summarizing_algorithm == 'consensus':
        primer_summarizing_algorithm = ConsensusBasedSummarizer()
    else:
        raise ValueError("Unknown primer summarizing algorithm specified.")

    if args.specificity_check_algorithm == 'blast':
        specificity_check_algorithm = SpecificityCheckBLAST()
        primer_specificity_params = {param.upper(): config.get('SpecificityParameters', param) for param in config.options('SpecificityParameters')}
        primer_specificity_params_converted = {key: convert_param(value) for key, value in primer_specificity_params.items()}
    else:
        raise ValueError("Unknown primer summarizing algorithm specified.")

    context.design_primers(primer_design_algorithm, primer_design_params_converted, primer_summarizing_algorithm, primer_summarizing_params, specificity_check_algorithm, primer_specificity_params_converted)

    ray.shutdown()
