#!/usr/bin/env python3

import os
import sys
import csv
import math
import shutil
import argparse
import subprocess
from collections import Counter, defaultdict
import itertools
import pandas as pd
import numpy as np
import primer3
from datasketch import WeightedMinHashGenerator
import time


@dataclass
class InputPaths:
    raw_dir: Optional[str] = None
    prokka_dir: Optional[str] = None
    panaroo_dir: Optional[str] = None


@dataclass
class Primer3Params:
    global_params: Dict[str, object] = field(default_factory=dict)
    design_params: Dict[str, object] = field(default_factory=dict)


@dataclass
class Config:
    input_type: str
    input_paths: InputPaths
    output_dir: str
    max_cores: int
    database_path: str
    primer3_config_file: Optional[str] = None
    primer3: Primer3Params = field(default_factory=Primer3Params)

    def __post_init__(self):
        allowed_inputs = {"raw", "prokka", "panaroo"}
        if self.input_type not in allowed_inputs:
            raise ValueError(f"`input_type` must be one of {allowed_inputs}, got: {self.input_type}")
        
        active_dir = {
            "raw": self.input_paths.raw_dir,
            "prokka": self.input_paths.prokka_dir,
            "panaroo": self.input_paths.panaroo_dir,
        }[self.input_type]
        if not active_dir:
            raise ValueError(f"{self.input_type}_dir must be set in `input_paths` for input_type = {self.input_type}")
        
        # Primer3 configuration logic
        if self.primer3_config_file and self.primer3.global_params:
            print("Warning: primer3_config_file is set. Inline primer3 parameters will be ignored.")


class ConfigLoader:
    @staticmethod
    def load(path: str) -> Config:
        with open(path, "r") as f:
            raw = yaml.safe_load(f)

        # Parse nested dataclasses manually
        raw["input_paths"] = InputPaths(**raw.get("input_paths", {}))
        raw["primer3"] = Primer3Params(**raw.get("primer3", {}))

        return Config(**raw)


class MSAStrategy:
    # --- Internal Methods ---
    
    def run_prokka(self, input_dir=".", output_dir="output/prokka"):
        """
        Runs Prokka on all FASTA files in the specified input directory, saving results in a single output directory.

        :param input_dir: Path to the directory containing FASTA files.
        :param output_dir: Path to the directory where Prokka output will be saved.
        """
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Iterate over all files in the input directory
        for i, filename in enumerate(os.listdir(input_dir)):
            if filename.endswith(".fasta") or filename.endswith(".fa"):
                entry_filepath = os.path.join(input_dir, filename)
                species_name = os.path.splitext(filename)[0]  # Remove file extension for the prefix

                # Construct Prokka command with the same output directory for all files
                prokka_cmd = [
                    'prokka',
                    '--kingdom', 'Bacteria',
                    '--outdir', output_dir,  # Use the same output directory
                    '--prefix', f"{species_name}_{i}",  # Unique prefix for each file
                    '--cpus', self.max_cores,
                    entry_filepath
                ]

                # Run Prokka
                print(f"Running Prokka for {filename} with command: {' '.join(prokka_cmd)}")
                subprocess.run(prokka_cmd, check=True)
                print(f"Prokka completed for {filename}")

        print("All Prokka runs completed.")

    

    # --- Internal Classes ---
    
    

    # --- Public Methods ---

    def __init__(self, config: Config):
        self.config = config

        self.input_type = config.input_type
        self.output_dir = config.output_dir
        self.max_cores = config.max_cores
        self.database_path = config.database_path

        self.input_dir = {
            "raw": config.input_paths.raw_dir,
            "prokka": config.input_paths.prokka_dir,
            "panaroo": config.input_paths.panaroo_dir
        }[self.input_type]

        if config.primer3_config_file:
            self.use_external_primer3_config = True
            self.primer3_config_path = config.primer3_config_file
        else:
            self.use_external_primer3_config = False
            self.primer3_global = config.primer3.global_params
            self.primer3_design = config.primer3.design_params        

    def run(self):
        print(f"Running with input type: {self.input_type}")
        print(f"Input directory: {self.input_dir}")
        print(f"Using {self.max_cores} cores")
        if self.use_external_primer3_config:
            print(f"Using external Primer3 config at {self.primer3_config_path}")
        else:
            print(f"Designing primers with Tm range: {self.primer3_design.get('PRIMER_MIN_TM')} - {self.primer3_design.get('PRIMER_MAX_TM')}")
            print(f"Product size range: {self.primer3_design.get('PRIMER_PRODUCT_SIZE_RANGE')}")
            print(f"Primer size range: {self.primer3_design.get('PRIMER_MIN_SIZE')} - {self.primer3_design.get('PRIMER_MAX_SIZE')}")
            print(f"Number of primers per locus: {self.primer3_global.get('PRIMER_NUM_RETURN')}")      


class QuasiAlignmentStrategy():
    # --- Internal Methods ---
    
    def parse_species_csv(self, file_path):
        """
        Parses a CSV file containing target species and their associated non-target species.

        :param file_path: Path to the input CSV file.
                          The file should have two columns: target species and non-target species (semicolon-separated).
        :return: A dictionary with target species as keys and lists of non-target species as values.
        """
        species_dict = {}

        with open(file_path, mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                # Skip empty rows
                if not row:
                    continue
                                
                target_species = row[0].strip()
                non_target_species = []
                
                if len(row) > 1:
                    non_target_species = [species.strip() for species in row[1].split(';') if species.strip()]
                
                species_dict[target_species] = non_target_species

        return species_dict
        
    @staticmethod
    def manhattan_distance(point1, point2):
        """
        Calculate the Manhattan distance between two points.

        :param point1: A list or tuple representing the coordinates of the first point.
        :param point2: A list or tuple representing the coordinates of the second point.
        :return: The Manhattan distance between the two points.
        """
        if len(point1) != len(point2):
            raise ValueError("Both points must have the same number of dimensions")

        distance = sum(abs(p1 - p2) for p1, p2 in zip(point1, point2))
        return distance

    @staticmethod
    def cosine_similarity(vector_a, vector_b):
        """
        Calculates the cosine similarity between two vectors.

        :param vector_a: A list or tuple representing the first vector.
        :param vector_b: A list or tuple representing the second vector.
        :return: The cosine similarity between vector_a and vector_b.
        """
        # Ensure inputs are lists or tuples
        if not isinstance(vector_a, (list, tuple)) or not isinstance(vector_b, (list, tuple)):
            raise TypeError("Input vectors must be lists or tuples.")
        
        # Convert to NumPy arrays
        vector_a = np.array(vector_a)
        vector_b = np.array(vector_b)
        
        # Calculate dot product and magnitudes
        dot_product = np.dot(vector_a, vector_b)
        magnitude_a = np.linalg.norm(vector_a)
        magnitude_b = np.linalg.norm(vector_b)
        
        # Avoid division by zero
        if magnitude_a == 0 or magnitude_b == 0:
            return 0  # Return 0 for cosine similarity if either vector is zero
        
        # Compute cosine similarity
        return dot_product / (magnitude_a * magnitude_b)
    
    @staticmethod 
    def jaccard_similarity_with_frequencies(vector_a, vector_b):
        """
        Calculates the Jaccard similarity between two lists or tuples containing frequencies of individual items (e.g., p-mers).

        :param list_a: A list or tuple representing the first vector.
        :param list_b: A list or tuple representing the second vector.
        :return: The Jaccard similarity between the two collections, accounting for frequencies.
        """
        # Convert lists to Counters to handle frequencies
        counter_a = Counter(vector_a)
        counter_b = Counter(vector_b)
        
        # Calculate intersection and union based on minimum and maximum frequencies
        intersection = sum((min(counter_a[item], counter_b[item]) for item in counter_a if item in counter_b))
        union = sum((max(counter_a[item], counter_b[item]) for item in set(counter_a) | set(counter_b)))
        
        # Avoid division by zero in case both counters are empty
        if union == 0:
            return 0
        
        # Compute Jaccard similarity with frequencies
        return intersection / union
    
    @staticmethod
    def calculate_weighted_minhash_similarity(frequency_vector_a, frequency_vector_b, num_perm=128):
        """
        Calculates the Weighted MinHash Jaccard similarity between two frequency vectors.

        :param frequency_vector_a: A list or tuple representing the frequency of p-mers for the first vector.
        :param frequency_vector_b: A list or tuple representing the frequency of p-mers for the second vector.
        :param num_perm: The number of permutations (hash functions) for the MinHash.
        :return: The Weighted MinHash Jaccard similarity between the two vectors.
        """
        if not isinstance(frequency_vector_a, (list, tuple)) or not isinstance(frequency_vector_b, (list, tuple)):
            raise TypeError("Input vectors must be lists or tuples.")
        
        if len(frequency_vector_a) != len(frequency_vector_b):
            raise ValueError("The two frequency vectors must be of the same length.")
        
        wmg = WeightedMinHashGenerator(len(frequency_vector_a), sample_size=num_perm)
        
        wm_a = wmg.minhash(frequency_vector_a)
        wm_b = wmg.minhash(frequency_vector_b)
        
        return wm_a.jaccard(wm_b)
    
    def run_prokka(self, input_dir=".", output_dir="output/prokka"):
        """
        Runs Prokka on all FASTA files in the specified input directory, saving results in a single output directory.

        :param input_dir: Path to the directory containing FASTA files.
        :param output_dir: Path to the directory where Prokka output will be saved.
        """
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Iterate over all files in the input directory
        for i, filename in enumerate(os.listdir(input_dir)):
            if filename.endswith(".fasta") or filename.endswith(".fa"):
                entry_filepath = os.path.join(input_dir, filename)
                species_name = os.path.splitext(filename)[0]  # Remove file extension for the prefix

                # Construct Prokka command with the same output directory for all files
                prokka_cmd = [
                    'prokka',
                    '--kingdom', 'Bacteria',
                    '--outdir', output_dir,  # Use the same output directory
                    '--prefix', f"{species_name}_{i}",  # Unique prefix for each file
                    '--cpus', "0",
                    entry_filepath
                ]

                # Run Prokka
                print(f"Running Prokka for {filename} with command: {' '.join(prokka_cmd)}")
                subprocess.run(prokka_cmd, check=True)
                print(f"Prokka completed for {filename}")

        print("All Prokka runs completed.")

    def run_panaroo_no_alignment(self, input_dir="output/prokka", output_dir="output/panaroo"):
        """
        Run Panaroo with no alignments and strict cleaning mode.

        Parameters:
            input_dir (str): Directory containing input files (default: "output/prokka").
            output_dir (str): Directory for Panaroo output (default: "output/panaroo").
        """
        # Construct the Panaroo command
        command = [
            "panaroo",
            "-i", f"{input_dir}/*.gff",
            "-o", output_dir,
            "--clean-mode", "strict"
        ]

        # Run the command
        try:
            subprocess.run(command, check=True)
            print(f"Panaroo completed successfully. Output is in '{output_dir}'.")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running Panaroo: {e}")
    
    def process_panaroo_output(self, presence_absence_file="output/panaroo/gene_presence_absence.csv", gene_data_file="output/panaroo/gene_data.csv"):
        """
        Process Panaroo output files and generate individual CSV files for each gene.

        This function processes two input CSV files generated by Panaroo:
        - `gene_presence_absence.csv`: Contains gene presence/absence information across samples.
        - `gene_data.csv`: Contains detailed information for each gene, including species, DNA sequence, and gene name.

        For each row in `gene_presence_absence.csv`, a separate CSV file is generated in the directory
        `output/panaroo/processed/`. The file is named after the gene name and contains the annotation ID, species name,
        DNA sequence, and gene name in the order: annotation ID, gene name, species name, DNA sequence.

        Parameters:
        ----------
        presence_absence_file : str
            The path to the `gene_presence_absence.csv` file.
        gene_data_file : str
            The path to the `gene_data.csv` file.

        Output:
        ------
        Multiple CSV files, one for each row in `gene_presence_absence.csv`, are created in the
        `output/panaroo/processed/` directory. Each file contains the following columns:
        - Sequence ID (== Annotation ID)
        - Gene Name
        - Species Name
        - DNA Sequence
        """
        
        # Ensure the output directory exists
        output_dir = 'output/panaroo/processed'
        os.makedirs(output_dir, exist_ok=True)

        # Step 1: Read the gene_data.csv into a dictionary for quick look-up
        gene_data = {}
        with open(gene_data_file, mode='r', newline='') as gene_data_csv:
            gene_data_reader = csv.reader(gene_data_csv)
            next(gene_data_reader)  # Skip the header

            for row in gene_data_reader:
                # Split the 4th column (annotation IDs) into individual IDs
                annotation_ids = row[3].split(";")
                for annotation_id in annotation_ids:
                    annotation_id = annotation_id.strip()
                    if annotation_id:  # Avoid empty annotation IDs
                        # print(annotation_id, row[0], row[6])
                        gene_data[annotation_id] = {
                            'species_name': row[0],  # 1st column - species name
                            'dna_sequence': row[5]  #,  # 6th column - DNA sequence
                            #'gene_name': row[6],     # 7th column - gene name
                        }

        # Step 2: Process the gene_presence_absence.csv file
        with open(presence_absence_file, mode='r', newline='') as presence_absence_csv:
            presence_absence_reader = csv.reader(presence_absence_csv)
            
            # Skip the header of the gene_presence_absence.csv file
            header = next(presence_absence_reader)

            for row in presence_absence_reader:
                gene_name = row[0]  # The gene name from the first column
                
                # Track if this gene has at least one valid annotation
                has_valid_annotation = False
                
                # Temporary storage for rows to write
                output_rows = []

                # Start looping from the 4th column onwards in gene_presence_absence.csv
                for cell in row[3:]:
                    if cell:  # Non-empty cell, meaning there's an annotation ID
                        annotation_ids = cell.split(";")
                        for annotation_id in annotation_ids:
                            annotation_id = annotation_id.strip()
                            # Step 3: Look up the annotation ID in the gene_data dictionary
                            if annotation_id in gene_data:
                                species_name = gene_data[annotation_id]['species_name']
                                dna_sequence = gene_data[annotation_id]['dna_sequence']
                                # gene_name = gene_data[annotation_id]['gene_name']
                                
                                # Add the row to output rows
                                output_rows.append([annotation_id, gene_name, species_name, dna_sequence])
                                has_valid_annotation = True  # At least one valid annotation found
                            else:
                                print(f"Warning: Annotation ID {annotation_id} not found in gene_data.csv")

                # Only create the file if there is at least one valid annotation
                if has_valid_annotation:
                    output_file_path = os.path.join(output_dir, f"{gene_name}.csv")
                    with open(output_file_path, mode='w', newline='') as output_csv:
                        output_writer = csv.writer(output_csv)
                        print("gene name: ", gene_name)
                        # Write the header for the output CSV file
                        output_writer.writerow(['SeqID', 'Gene Name', 'Species Name', 'DNA Sequence'])
                        # Write the rows with valid annotations
                        output_writer.writerows(output_rows)
    
    def create_quasialignments(self, species_dict, segment_length, step, distance_func, threshold, input_dir='output/panaroo/processed'):
        """
        Processes all CSV files in the specified directory by removing the '_\\d+' suffix from species names in the third column.
        For each target species in the dictionary, selects rows where the species name in the third column matches
        the target species or any species in the associated list of non-target species. Then extracts segments
        using the extract_segments_from_list function. Finally, clusters segments based on p-mer profiles.

        :param species_dict: A dictionary where keys are target species and values are lists of non-target species.
        :param segment_length: The length of each segment to cut from the sequences.
        :param step: The number of positions to move to the right after each segment is cut.
        :param distance_func: A function to calculate the distance between two pmer_profiles.
        :param threshold: The maximum allowable distance for adding a Segment to an existing QuasiAlignment.
        :param input_dir: The directory containing the CSV files to be processed (default is 'output/panaroo/processed').
        :return: A dictionary where keys are file names and values are dictionaries with positions as keys
                 and lists of QuasiAlignment objects as values.
        """
        result_dict = {}

        # Loop through all CSV files in the specified directory
        for file_name in os.listdir(input_dir):
            if file_name.endswith(".csv"):
                print(file_name)
                file_path = os.path.join(input_dir, file_name)
                
                # Load the CSV file into a DataFrame
                df = pd.read_csv(file_path)
                
                # Remove the suffix '_\d+' from the third column (Species name)
                df.iloc[:, 2] = df.iloc[:, 2].str.replace(r'_\d+', '', regex=True)
                
                # Initialize an empty DataFrame to store filtered rows for this file
                filtered_df = pd.DataFrame()
                
                print(species_dict.keys())
                
                # Loop through the keys (target species) and filter rows based on species name
                for target_species, non_target_species_list in species_dict.items():
                    print("test")
                    species_to_keep = [target_species] + non_target_species_list
                    temp_df = df[df.iloc[:, 2].isin(species_to_keep)]
                    filtered_df = pd.concat([filtered_df, temp_df], ignore_index=True)
                
                print(species_to_keep)
                
                # Prepare the data for extract_segments_from_list
                entries = []
                for _, row in filtered_df.iterrows():
                    species_name = row[2]  # Third column is the species name
                    gene_name = row[1]     # Second column is the gene name
                    seq_id = row[0]        # First column is the sequence ID
                    sequence = row[3]      # Fourth column is the DNA sequence
                    entries.append((species_name, gene_name, seq_id, sequence))
                
                # Extract segments from the filtered data
                segments_dict = self.extract_segments_from_list(entries, segment_length, step)
                
                # Initialize a dictionary to store the QuasiAlignment objects for each position
                file_quasialignments = {}

                # Loop through each position in segments_dict
                for position, segment_list in segments_dict.items():
                    quasi_alignments = []  # List to hold quasi-alignments for the current position

                    # Process each Segment in the segment_list at this position
                    for segment in segment_list:
                        # Use add_segment_to_quasialignment to assign the segment to an existing or new cluster
                        quasi_alignments = self.assign_segment_to_closest_quasialignment(quasi_alignments, segment, distance_func, threshold)
                    
                    # Store the list of QuasiAlignments for this position
                    file_quasialignments[position] = quasi_alignments
                
                # Store the clusters for this file in the result dictionary
                result_dict[file_name] = file_quasialignments
        
        return result_dict
    
    @staticmethod
    def extract_segments(file_path, segment_length, step):
        """
        Processes a CSV file to cut segments of a specified length from each sequence in the file.
        It cuts the first segment starting at the beginning, then moves to the right by a specified 
        number of positions (step) and cuts the second segment and so forth until it reaches the end of the sequence.
        Each segment is stored using the Segment class, containing gene name, species name, and sequence ID, 
        and the dictionary returned will have keys as positions within the original sequence, and values as Segment objects.

        :param file_path: The path to the CSV file with columns 'SeqID', 'Gene Name', 'Species Name', and 'DNA Sequence'.
        :param segment_length: The length of each segment to cut from the sequences.
        :param step: The number of positions to move to the right after each segment is cut.
        :return: A dictionary where the keys are the starting positions within the sequence, and the values are Segment objects.
        """
        all_segments = {}
        
        # Open and read the CSV file
        with open(file_path, mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Loop through each row in the CSV
            for row in reader:
                sequence = row['DNA Sequence']
                gene_name = row['Gene Name']
                species_name = row['Species Name']
                seq_id = row['SeqID']
                
                # Loop through the sequence and extract segments of the specified length
                for i in range(0, len(sequence), step):
                    # Ensure we don't go beyond the end of the sequence
                    if i + segment_length <= len(sequence):
                        segment = sequence[i:i + segment_length]
                        
                        # Create a Segment object for each extracted segment
                        segment_obj = Segment(seq_id=seq_id, species=species_name, gene_name=gene_name, segment=segment)
                        
                        # Store the Segment object in the dictionary with the starting position as the key
                        all_segments[i] = segment_obj
    
        return all_segments
    
    def extract_segments_from_list(self, entries, segment_length, step):
        """
        Processes a list of sequences and extracts segments of a specified length.
        It cuts the first segment starting at the beginning, then moves to the right by a specified 
        number of positions (step) and cuts the second segment, and so forth until it reaches the end of the sequence.
        Each segment is stored in a dictionary where the keys are starting positions, and the values
        are lists of Segment objects (to handle multiple segments starting at the same position).

        :param entries: A list of tuples, where each tuple contains species name, gene name, sequence ID, and the DNA sequence.
        :param segment_length: The length of each segment to cut from the sequences.
        :param step: The number of positions to move to the right after each segment is cut.
        :return: A dictionary where the keys are starting positions within the sequence, and the values are lists of Segment objects.
        """
        all_segments = {}

        # Loop through the provided list of entries
        for species_name, gene_name, seq_id, sequence in entries:
            # Ensure the sequence is in uppercase
            sequence = sequence.upper()

            # Loop through the sequence and extract segments of the specified length
            for i in range(0, len(sequence), step):
                # Ensure we don't go beyond the end of the sequence
                if i + segment_length <= len(sequence):
                    segment = sequence[i:i + segment_length]
                    
                    # Create a Segment object for each extracted segment
                    segment_obj = self.Segment(seq_id=seq_id, species_name=species_name, gene_name=gene_name, start=i, length=segment_length, sequence=segment)
                    
                    # If the position is already a key, append the segment to the list
                    if i in all_segments:
                        all_segments[i].append(segment_obj)
                    else:
                        # If it's a new position, create a list and add the segment
                        all_segments[i] = [segment_obj]
        
        return all_segments
    
    def assign_segment_to_closest_quasialignment(self, quasi_alignments, segment, distance_func, threshold):
        """
        Adds a Segment object to the closest quasi-alignment based on its distance to each quasi-alignment's medoid.
        If no quasi-alignment is found within the threshold distance, creates a new QuasiAlignment object and adds the Segment object to it.

        :param quasi_alignments: A list of QuasiAlignment objects.
        :param segment: A Segment object to be added.
        :param distance_func: A function to calculate the distance between two pmer_profiles.
        :param threshold: The maximum allowable distance for adding the Segment to an existing QuasiAlignment.
        :return: The updated list of QuasiAlignment objects.
        """
        best_fit_quasi_alignment = None
        min_distance = float('inf')
        
        # Iterate through each QuasiAlignment to find the best fit within the threshold
        for quasi_alignment in quasi_alignments:
            if quasi_alignment.medoid is not None:
                # Calculate the distance between the segment and the medoid of the current QuasiAlignment
                distance = segment.get_distance(quasi_alignment.medoid, self.manhattan_distance)
                
                # Check if this QuasiAlignment is a better fit (within threshold and lower distance)
                if distance < threshold and distance < min_distance:
                    best_fit_quasi_alignment = quasi_alignment
                    min_distance = distance
        
        # Add the segment to the best-fit QuasiAlignment or create a new one if none are within threshold
        if best_fit_quasi_alignment:
            best_fit_quasi_alignment.add_segment(segment)
        else:
            # Create a new QuasiAlignment and add the segment to it
            new_quasi_alignment = QuasiAlignmentStrategy.QuasiAlignment(int(time.time_ns()))
            new_quasi_alignment.add_segment(segment)
            quasi_alignments.append(new_quasi_alignment)
        
        return quasi_alignments

    def evaluate_informative_positions(self, quasialignment_data, species_dict, threshold, proportion_threshold,
                                   output_dir="output/panaroo/processed/selected"):
        """
        Evaluates each file based on the proportion of informative positions per target species.
        A file is selected only if all target species meet the informative threshold.
        Outputs a CSV summary with proportions per species per file.

        :param quasialignment_data: Dict[filename -> Dict[position -> List[QuasiAlignment]]]
        :param species_dict: Dict[target_species -> List[non_target_species]]
        :param threshold: Proportion of target segments required at a position to be informative.
        :param proportion_threshold: Proportion of informative positions required per species.
        :param output_dir: Path to store selected files and the summary CSV.
        :return: Dict[filename -> Dict[target_species -> proportion of informative positions]]
        """
        file_species_proportions = {}

        os.makedirs(output_dir, exist_ok=True)

        for file_name, positions_dict in quasialignment_data.items():
            print(f"Evaluating file: {file_name}")
            total_positions = len(positions_dict)
            species_informative_counts = {species: 0 for species in species_dict}

            # Check informativeness per position per target species
            for position, quasi_alignments in positions_dict.items():
                for target_species, non_target_species_list in species_dict.items():
                    total_target_segments = 0
                    target_only_count = 0

                    for quasi_alignment in quasi_alignments:
                        contains_non_target = any(
                            segment.species_name in non_target_species_list for segment in quasi_alignment.segments
                        )

                        target_segments = [s for s in quasi_alignment.segments if s.species_name == target_species]
                        total_target_segments += len(target_segments)

                        if not contains_non_target:
                            target_only_count += len(target_segments)

                    if total_target_segments > 0:
                        target_proportion = target_only_count / total_target_segments
                        if target_proportion >= threshold:
                            species_informative_counts[target_species] += 1

            # Calculate proportions and store them
            species_proportions = {
                species: species_informative_counts[species] / total_positions if total_positions > 0 else 0
                for species in species_dict
            }
            file_species_proportions[file_name] = species_proportions

            # Check if all species meet the proportion threshold
            if all(p >= proportion_threshold for p in species_proportions.values()):
                source_file_path = os.path.join("output/panaroo/processed", os.path.basename(file_name))
                destination_file_path = os.path.join(output_dir, os.path.basename(file_name))
                shutil.copy(source_file_path, destination_file_path)
                print(f"{file_name} was evaluated as informative for all species.")

        # Write summary CSV
        summary_csv_path = os.path.join(output_dir, "informative_positions.csv")
        with open(summary_csv_path, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = ['File'] + list(species_dict.keys())
            writer.writerow(header)

            for file_name, proportions in file_species_proportions.items():
                row = [file_name] + [f"{proportions[species]:.4f}" for species in species_dict]
                writer.writerow(row)

        return file_species_proportions

    def calculate_consensus(self, primer_list):
        """Calculates the consensus sequence from a list of primer sequences."""
        if not primer_list:
            return ""
        consensus = ""
        for position_bases in zip(*primer_list):
            base_counts = Counter(position_bases)
            most_common_base, freq = base_counts.most_common(1)[0]
            consensus += self.get_iupac_code(base_counts) if len(base_counts) > 1 else most_common_base
        return consensus

    def get_iupac_code(self, base_counts):
        """Gets the IUPAC code for degenerate positions based on base counts."""
        iupac_dict = {
            frozenset(['A']): 'A', frozenset(['C']): 'C', frozenset(['G']): 'G', frozenset(['T']): 'T',
            frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y', frozenset(['G', 'C']): 'S',
            frozenset(['A', 'T']): 'W', frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
            frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D', frozenset(['A', 'C', 'T']): 'H',
            frozenset(['A', 'C', 'G']): 'V', frozenset(['A', 'T', 'C', 'G']): 'N'
        }
        return iupac_dict[frozenset(base_counts.keys())]
    
    @staticmethod
    def add_missing_keys_with_zero(dict1, dict2):
        """
        Adds keys that are present in dict1 but not in dict2 (and vice versa) to each dictionary with values set to zero.
        
        :param dict1: The first dictionary.
        :param dict2: The second dictionary.
        :return: Two dictionaries with added keys and values set to zero where keys were missing.
        """
        # Find keys unique to each dictionary
        keys_in_dict1_not_in_dict2 = dict1.keys() - dict2.keys()
        keys_in_dict2_not_in_dict1 = dict2.keys() - dict1.keys()
        
        # Add missing keys to each dictionary with values set to zero
        for key in keys_in_dict1_not_in_dict2:
            dict2[key] = 0
        for key in keys_in_dict2_not_in_dict1:
            dict1[key] = 0
        
        return dict1, dict2
        
    def find_consensus_binding_regions(self, all_primers, window_size=10, threshold_proportion=0.5):
        """
        Identifies consensus binding regions for forward and reverse primers separately
        where a threshold proportion of primers bind across sequences.

        :param all_primers: Dictionary with sequence IDs as keys and lists of primer dictionaries as values.
        :param window_size: The range of positions to consider in the histogram for grouping similar binding sites.
        :param threshold_proportion: The minimum proportion of sequences that need to bind near a position for consensus.
        :return: Two lists of position ranges representing consensus binding regions for forward and reverse primers.
        """
        forward_position_histogram = defaultdict(int)
        reverse_position_histogram = defaultdict(int)
        total_sequences = len(all_primers)

        for seq_id, primers in all_primers.items():
            for primer in primers:
                left_pos = primer['left_pos']
                right_pos = primer['right_pos']
                
                for pos in range(left_pos - window_size // 2, left_pos + window_size // 2 + 1):
                    forward_position_histogram[pos] += 1
                
                for pos in range(right_pos - window_size // 2, right_pos + window_size // 2 + 1):
                    reverse_position_histogram[pos] += 1

        forward_consensus_regions = []
        current_region = []
        for pos in sorted(forward_position_histogram):
            count = forward_position_histogram[pos]
            if count / total_sequences >= threshold_proportion:
                if not current_region:
                    current_region = [pos, pos]
                else:
                    current_region[1] = pos
            else:
                if current_region:
                    forward_consensus_regions.append(tuple(current_region))
                    current_region = []
        if current_region:
            forward_consensus_regions.append(tuple(current_region))

        reverse_consensus_regions = []
        current_region = []
        for pos in sorted(reverse_position_histogram):
            count = reverse_position_histogram[pos]
            if count / total_sequences >= threshold_proportion:
                if not current_region:
                    current_region = [pos, pos]
                else:
                    current_region[1] = pos
            else:
                if current_region:
                    reverse_consensus_regions.append(tuple(current_region))
                    current_region = []
        if current_region:
            reverse_consensus_regions.append(tuple(current_region))

        return forward_consensus_regions, reverse_consensus_regions
        
    def select_consensus_primers_old(self, primers, forward_regions, reverse_regions):
        # Pre-compute the best reverse primer for each reverse consensus region
        reverse_best_primers = {}
        for rev_region_start, rev_region_end in reverse_regions:
            reverse_candidates = [
                p for p in primers if rev_region_start <= p['right_pos'] <= rev_region_end
            ]
            if reverse_candidates:
                reverse_best_primers[(rev_region_start, rev_region_end)] = min(
                    reverse_candidates,
                    key=lambda p: abs(p['right_pos'] - (rev_region_start + rev_region_end) // 2)
                )
            else:
                reverse_best_primers[(rev_region_start, rev_region_end)] = None

        # Process forward regions and pair with pre-computed reverse primers
        for region_start, region_end in forward_regions:
            # Select the best forward primer
            forward_candidates = [
                p for p in primers if region_start <= p['left_pos'] <= region_end
            ]
            if forward_candidates:
                selected_forward = min(
                    forward_candidates,
                    key=lambda p: abs(p['left_pos'] - (region_start + region_end) // 2)
                )
            else:
                selected_forward = None

            # Loop through pre-computed reverse primers
            for (rev_region_start, rev_region_end), selected_reverse in reverse_best_primers.items():
                # Store the pair of forward and reverse primers
                consensus_primers[seq_id].append({
                    'forward_region': (region_start, region_end),
                    'reverse_region': (rev_region_start, rev_region_end),
                    'forward': selected_forward['left_seq'] if selected_forward else None,
                    'reverse': selected_reverse['right_seq'] if selected_reverse else None
                })
        
        return consensus_primers
        
    def select_consensus_primers(self, primers, forward_regions, reverse_regions):
        """
        Selects consensus primers for forward and reverse regions for each sequence.

        :param primers: A dictionary where keys are sequence IDs and values are lists of dictionaries, 
                        each describing individual primer pairs.
        :param forward_regions: A list of tuples representing forward consensus regions (start, end).
        :param reverse_regions: A list of tuples representing reverse consensus regions (start, end).
        :return: A dictionary where keys are sequence IDs and values are lists of selected primer pairs for consensus regions.
        """
        consensus_primers = {}

        for seq_id, primer_list in primers.items():
            consensus_primers[seq_id] = []

            # Pre-compute the best reverse primer for each reverse consensus region
            reverse_best_primers = {}
            for rev_region_start, rev_region_end in reverse_regions:
                reverse_candidates = [
                    p for p in primer_list if rev_region_start <= p['right_pos'] <= rev_region_end
                ]
                if reverse_candidates:
                    reverse_best_primers[(rev_region_start, rev_region_end)] = min(
                        reverse_candidates,
                        key=lambda p: abs(p['right_pos'] - (rev_region_start + rev_region_end) // 2)
                    )
                else:
                    reverse_best_primers[(rev_region_start, rev_region_end)] = None

            # Process forward regions and pair with pre-computed reverse primers
            for region_start, region_end in forward_regions:
                # Select the best forward primer
                forward_candidates = [
                    p for p in primer_list if region_start <= p['left_pos'] <= region_end
                ]
                if forward_candidates:
                    selected_forward = min(
                        forward_candidates,
                        key=lambda p: abs(p['left_pos'] - (region_start + region_end) // 2)
                    )
                else:
                    selected_forward = None

                # Loop through pre-computed reverse primers
                for (rev_region_start, rev_region_end), selected_reverse in reverse_best_primers.items():
                    # Store the pair of forward and reverse primers
                    consensus_primers[seq_id].append({
                        'forward_region': (region_start, region_end),
                        'reverse_region': (rev_region_start, rev_region_end),
                        'forward': selected_forward['left_seq'] if selected_forward else None,
                        'reverse': selected_reverse['right_seq'] if selected_reverse else None
                    })

        return consensus_primers
        
    def generate_degenerate_primers(self, consensus_primers):
        """
        Generates degenerate primers for each consensus binding region across all sequences.

        :param consensus_primers: The output dictionary from select_consensus_primers.
                                  Keys are sequence IDs, values are lists of primer dictionaries.
        :return: A list of dictionaries, each containing degenerate forward and reverse primers for a consensus region.
        """
        degenerate_primers = []

        # Collect all primers grouped by consensus region
        region_primers = {}

        for seq_id, regions in consensus_primers.items():
            for region in regions:
                forward_region = region['forward_region']
                reverse_region = region['reverse_region']

                # Use the region as a key to group primers
                region_key = (forward_region, reverse_region)
                if region_key not in region_primers:
                    region_primers[region_key] = {'forward': [], 'reverse': []}

                # Add primers from this sequence to the region's primer pool
                if region['forward']:
                    region_primers[region_key]['forward'].append(region['forward'])
                if region['reverse']:
                    region_primers[region_key]['reverse'].append(region['reverse'])

        # Generate degenerate primers for each region
        for (forward_region, reverse_region), primers in region_primers.items():
            forward_consensus = self.calculate_consensus(primers['forward']) if primers['forward'] else None
            reverse_consensus = self.calculate_consensus(primers['reverse']) if primers['reverse'] else None

            degenerate_primers.append({
                'forward_region': forward_region,
                'reverse_region': reverse_region,
                'forward_degenerate': forward_consensus,
                'reverse_degenerate': reverse_consensus
            })

        return degenerate_primers
        
    def create_species_dict(self, csv_file):
        """
        Reads a CSV file with target and non-target species and creates a dictionary where
        keys are target species and values are lists of non-target species.

        :param csv_file: Path to the CSV file with two columns: target species and semicolon-separated non-target species.
        :return: A dictionary with target species as keys and lists of non-target species as values.
        """
        species_dict = {}

        with open(csv_file, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.reader(file)

            for row in reader:
                # Skip empty rows
                if not row:
                    continue

                target_species = row[0].strip()
                non_target_species = []

                # Check if the row has a second column
                if len(row) > 1:
                    non_target_species = [species.strip() for species in row[1].split(';') if species.strip()]

                species_dict[target_species] = non_target_species

        return species_dict
        
    def parse_primer3_config(self, file_path):
        """
        Reads a Primer3 configuration file and stores the parameters in a dictionary.
        
        :param file_path: Path to the Primer3 configuration file.
        :return: A dictionary containing the Primer3 parameters.
        """
        primer3_params = {}

        # Open and read the file line by line
        with open(file_path, 'r') as file:
            for line in file:
                # Ignore comments and empty lines
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # Split the line into key and value
                if "=" in line:
                    key, value = line.split("=", 1)
                    primer3_params[key.strip()] = value.strip()
        
        return primer3_params
        
    def design_primers_from_csv_api(self, input_csv, output_csv):
        """
        Designs primers for each sequence in a CSV file using Primer3, identifies consensus primer regions,
        and generates degenerate primers for each target species.

        :param input_csv: Path to the input CSV file containing columns 'Sequence ID', 'Gene Name', 'Species Name', 'DNA Sequence'.
        :param output_csv: Path to the output CSV file to save degenerate primer sequences.
        """
        df = pd.read_csv(input_csv)
        species_primers = {}
        
        df.iloc[:, 2] = df.iloc[:, 2].str.replace(r'_\d+', '', regex=True)

        # Step 1: Collect primers for each sequence ID
        for target_species in self.species_dict.keys():
            species_primers[target_species] = {}

            # Filter the DataFrame to get sequences for the target species
            target_df = df[df['Species Name'] == target_species]
            print("design of primers1", target_species, df['Species Name'])
            for _, row in target_df.iterrows():
                sequence_id = row['SeqID']
                sequence = row['DNA Sequence']
                print("design of primers")
                # Set sequence in primer parameters and design primers
                self.primer_parameters['SEQUENCE_TEMPLATE'] = sequence
                print(self.primer_parameters)
                primers = primer3.bindings.design_primers({}, self.primer_parameters)

                # Store primers if available
                if 'PRIMER_LEFT_0_SEQUENCE' in primers and 'PRIMER_RIGHT_0_SEQUENCE' in primers:
                    primer_entry = {
                        'left_seq': primers['PRIMER_LEFT_0_SEQUENCE'],
                        'right_seq': primers['PRIMER_RIGHT_0_SEQUENCE'],
                        'left_pos': primers['PRIMER_LEFT_0'][0],
                        'right_pos': primers['PRIMER_RIGHT_0'][0]
                    }
                    if sequence_id not in species_primers[target_species]:
                        species_primers[target_species][sequence_id] = []
                    species_primers[target_species][sequence_id].append(primer_entry)
                    print(primer_entry)

        # Step 2: Identify consensus binding regions for each species
        degenerate_primers = {}
        for target_species, primer_dict in species_primers.items():
            # Identify forward and reverse consensus binding regions
            forward_regions, reverse_regions = self.find_consensus_binding_regions(primer_dict)
            
            # Generate degenerate primers based on identified regions
            degenerate_primers[target_species] = self.generate_degenerate_primers(
                self.select_consensus_primers(primer_dict, forward_regions, reverse_regions)
            )

        # Step 3: Write degenerate primers to output CSV
        with open(output_csv, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['Species', 'Primer Type', 'Consensus Sequence'])
            for species, primer_data in degenerate_primers.items():
                for region in primer_data:
                    writer.writerow([species, 'Forward', region['forward_degenerate']])
                    writer.writerow([species, 'Reverse', region['reverse_degenerate']])
                    
    def design_primers_from_csv(self, input_csv, output_csv):
        """
        Designs primers for each sequence in a CSV file using Primer3_core, identifies consensus primer regions,
        and generates degenerate primers for each target species.

        :param input_csv: Path to the input CSV file containing columns 'Sequence ID', 'Gene Name', 'Species Name', 'DNA Sequence'.
        :param output_csv: Path to the output CSV file to save degenerate primer sequences.
        """
        df = pd.read_csv(input_csv)
        species_primers = {}
        df.iloc[:, 2] = df.iloc[:, 2].str.replace(r'_\d+', '', regex=True)

        # Step 1: Collect primers for each sequence ID
        for target_species in self.species_dict.keys():
            species_primers[target_species] = {}

            # Filter the DataFrame to get sequences for the target species
            target_df = df[df['Species Name'] == target_species]
            print("Designing primers for:", target_species)

            for _, row in target_df.iterrows():
                sequence_id = row['SeqID']
                sequence = row['DNA Sequence']

                # Prepare Primer3 input file
                primer3_input_file = "primer3_input.txt"
                primer3_output_file = "primer3_output.txt"

                with open(primer3_input_file, "w") as f:
                    f.write(f"SEQUENCE_ID={sequence_id}\n")
                    f.write(f"SEQUENCE_TEMPLATE={sequence}\n")
                    for param, value in self.primer_parameters.items():
                        f.write(f"{param}={value}\n")
                    f.write("=\n")  # End of the Primer3 input file

                # Run Primer3
                try:
                    subprocess.run(
                        ["primer3_core", primer3_input_file, "--output=" + primer3_output_file],
                        check=True,
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error running primer3_core for {sequence_id}: {e}")
                    continue

                # Parse Primer3 output
                primers = {}
                with open(primer3_output_file, "r") as f:
                    for line in f:
                        if line.startswith("PRIMER_LEFT_0_SEQUENCE"):
                            primers["left_seq"] = line.strip().split("=")[1]
                        elif line.startswith("PRIMER_RIGHT_0_SEQUENCE"):
                            primers["right_seq"] = line.strip().split("=")[1]
                        elif line.startswith("PRIMER_LEFT_0="):
                            primers["left_pos"] = int(line.strip().split("=")[1].split(",")[0])
                        elif line.startswith("PRIMER_RIGHT_0="):
                            primers["right_pos"] = int(line.strip().split("=")[1].split(",")[0])

                # Store primers if available
                if "left_seq" in primers and "right_seq" in primers:
                    primer_entry = {
                        "left_seq": primers["left_seq"],
                        "right_seq": primers["right_seq"],
                        "left_pos": primers["left_pos"],
                        "right_pos": primers["right_pos"],
                    }
                    if sequence_id not in species_primers[target_species]:
                        species_primers[target_species][sequence_id] = []
                    species_primers[target_species][sequence_id].append(primer_entry)
                
        print(species_primers)

        # Step 2: Identify consensus binding regions for each species
        degenerate_primers = {}
        for target_species, primer_dict in species_primers.items():
            # Identify forward and reverse consensus binding regions
            forward_regions, reverse_regions = self.find_consensus_binding_regions(primer_dict)
            # print(forward_regions)

            # Generate degenerate primers based on identified regions
            degenerate_primers[target_species] = self.generate_degenerate_primers(
                self.select_consensus_primers(primer_dict, forward_regions, reverse_regions)
            )

        # Step 3: Write degenerate primers to output CSV
        with open(output_csv, mode="w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Species", "Primer Type", "Consensus Sequence"])
            for species, primer_data in degenerate_primers.items():
                for region in primer_data:
                    writer.writerow([species, "Forward", region["forward_degenerate"]])
                    writer.writerow([species, "Reverse", region["reverse_degenerate"]])

        # Cleanup temporary files
        os.remove(primer3_input_file)
        os.remove(primer3_output_file)
        
    # --- Internal Classes ---
    
    class Segment:
        """
        A class to represent an individual segment in a quasi-alignment. It includes information about the sequence ID, species, 
        and the position of the segment within the sequence.
        """
        def __init__(self, seq_id: str, species_name: str, gene_name: str, start: int, length: int, sequence: str, p: int = 3):
            """
            Initialize a Segment object.

            :param seq_id: The ID of the sequence from which the segment comes.
            :param species_name: The species of the sequence.
            :param gene_name: The gene region of the sequence.
            :param start: The starting position of the segment within the sequence.
            :param length: The length of the segment within the sequence.
            """
            self.seq_id = seq_id
            self.species_name = species_name
            self.gene_name = gene_name
            self.start = start      # 0-indexed
            self.length = length
            self.sequence = sequence
            self.pmer_profile = self.get_pmer_composition(self, p)
            
        def get_segment(self, sequence):
            """Returns the segment of the sequence based on start and length."""
            end = self.start + self.length
            return sequence[self.start:end]

        def get_pmer_composition(self, segment, p):
            """
            Takes a Segment object and computes the p-mer composition of the segment.
            The output is a dictionary where keys are only the p-mers present in the segment's sequence, and values are their counts.
            P-mers are stored in uppercase.

            :param segment: A Segment object containing the start and length information.
            :param p: The length of the p-mer (substring).
            :return: A dictionary with present p-mers as keys and their counts as values.
            """
            sequence = segment.sequence.upper()  # Convert the entire sequence to uppercase
            pmer_dict = {}

            # Loop through the sequence and extract p-mers
            for i in range(len(sequence) - p + 1):  # Ensure we don't go beyond the sequence
                pmer = sequence[i:i + p]
                
                # Add or update the count of the p-mer in the dictionary
                if pmer in pmer_dict:
                    pmer_dict[pmer] += 1
                else:
                    pmer_dict[pmer] = 1
            
            return pmer_dict
        
        def get_distance(self, other_segment, distance_func):
            """
            Calculates the distance between this segment and another segment using a specified distance function.
            
            :param other_segment: The other Segment object to calculate the distance to.
            :param distance_func: A function that takes two lists (or tuples) as input and returns the distance between them.
            :return: The distance between the two p-mer profiles.
            """
            if not callable(distance_func):
                raise ValueError("distance_func must be a callable function")
            
            aligned_profile_self, aligned_profile_other = QuasiAlignmentStrategy.add_missing_keys_with_zero(self.pmer_profile, other_segment.pmer_profile)
            
            profile_values_self = tuple(aligned_profile_self.values())
            profile_values_other = tuple(aligned_profile_other.values())
            
            # Calculate and return the distance using the specified distance function
            return distance_func(profile_values_self, profile_values_other)

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
            self.medoid = None

        def add_segment(self, new_segment: "QuasiAlignmentStrategy.Segment"):
            """
            Adds a new Segment to the quasi-alignment and recalculates the medoid.
            The medoid is the segment with the smallest average distance to all other segments.
            
            :param new_segment: Segment object to be added to the quasi-alignment.
            """
            # Add the new segment to the list
            self.segments.append(new_segment)
            
            # Recalculate the medoid if there is more than one segment
            if len(self.segments) > 1:
                # Calculate the average distance of each segment to all others
                min_avg_distance = float('inf')
                new_medoid = None

                for segment in self.segments:
                    distances = [segment.get_distance(other, QuasiAlignmentStrategy.manhattan_distance) for other in self.segments if other != segment]
                    avg_distance = np.mean(distances)
                    
                    # Update the medoid if a smaller average distance is found
                    if avg_distance < min_avg_distance:
                        min_avg_distance = avg_distance
                        new_medoid = segment

                # Set the new medoid
                self.medoid = new_medoid
            else:
                # If only one segment, it's the medoid by default
                self.medoid = new_segment

        def get_segments_by_species(self, species: str):
            """
            Get all segments that belong to a specific species.

            :param species: The species to filter segments by.
            :return: A list of segments that belong to the specified species.
            """
            return [segment for segment in self.segments if segment.species_name == species]
            
    # --- Public Methods ---
            
    def __repr__(self):
        return f"QuasiAlignment(cluster_id={self.cluster_id}, segments={self.segments}, medoid={self.medoid})"

    def __init__(self, species_dict_file, primer_parameters_file):
        self.species_dict = self.parse_species_csv(species_dict_file)
        self.primer_parameters = self.parse_primer3_config(primer_parameters_file)
        
    def design_primers(self, input_folder=".", output_folder="output/primers", starting_point="genome_fasta_files"):
        """
        Runs the primer design process for all CSV files in the input folder.

        :param input_folder: Directory containing input CSV files (default: 'output/panaroo/processed/selected').
        :param output_folder: Directory where output CSV files will be saved (default: 'output/primers').
        """
        if starting_point == "genome_fasta_files":
            self.run_prokka(input_folder)
            self.run_panaroo_no_alignment()

        if starting_point == "genome_fasta_files" or starting_point == "panaroo_output":
            print("Processing Panaroo output...")
            self.process_panaroo_output()
        
        if starting_point == "genome_fasta_files" or starting_point == "panaroo_output" or starting_point == "processed_panaroo_output":
            print("Creating quasialignments...")
            quasialignment_data = self.create_quasialignments(self.species_dict, segment_length=100, step=50,
                                                                  distance_func=self.manhattan_distance, threshold=30)
            print("Evaluating informative positions...")
            informative_positions = self.evaluate_informative_positions(quasialignment_data, self.species_dict,
                                                                        threshold=0.9, proportion_threshold=0.8)

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        for filename in os.listdir(input_folder + "/output/panaroo/processed/selected"):
            if filename.endswith(".csv"):
                input_csv = os.path.join("output/panaroo/processed/selected", filename)
                output_csv = os.path.join(output_folder, f"primers_{filename}")
                print(f"Designing primers for {filename}")
                self.design_primers_from_csv(input_csv, output_csv)
        print("All primer designs completed.")
  
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    
    parser.add_argument('-i', '--input_directory', default=".", help='An input directory.')
    parser.add_argument('-s', '--species_csv', required=True, help='A CSV file with target and non-target species.')
    parser.add_argument('-p', '--primer3_params', required=True, help='The Primer3 config file.')
    parser.add_argument('-S', '--starting_point', default="panaroo_output", help='A starting point of the analysis. [genome_fasta_files | panaroo_output | processed_panaroo_output | selected_marker_regions]')

    args = parser.parse_args()

    strategy = QuasiAlignmentStrategy(args.species_csv, args.primer3_params)
    strategy.design_primers(args.input_directory, starting_point=args.starting_point)

