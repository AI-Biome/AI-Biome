import numpy as np
import pandas as pd
from tempfile import NamedTemporaryFile
# import fasta2pandas as f2p
from main import read_multifasta
#import shannonEntropy as se
from main import MarkerLociIdentificationStrategy
import math
import pandas as pd
#from quasiAlignmentCompatible import align_sequences

# Create an instance of the class MarkerLociIdentificationStrategy in order to use find_conserved_regions_shannon_entropy and align_sequences
analyzer = MarkerLociIdentificationStrategy()

def manhattan_distance(vec1, vec2):
    return np.sum(np.abs(vec1 - vec2))

def calculate_centroid_mean(vectors):
    return np.mean(vectors, axis=0)

# def trim_sequence(source_df,seq_id,beginning,end,continuation_len_thr,segment_size):
#     """
#     based on beginning and end points, get a new sequence fragment and add a few bases at the beginning and end (determined by continuation_len_thr) if in bound of the sequence
#     """
#     seq_beginning = beginning - round(0.5*continuation_len_thr*segment_size) if (beginning - round(0.5*continuation_len_thr*segment_size) > 0) else 0
#     seq_end = end + round(0.5*continuation_len_thr*segment_size) if (end + round(0.5*continuation_len_thr*segment_size) < len(source_df.loc[source_df["genome_id"]==seq_id]["sequence"].iloc[0])) else len(source_df.loc[source_df["genome_id"]==seq_id]["sequence"].iloc[0])
#     sequence = source_df.loc[source_df["genome_id"]==seq_id]["sequence"].iloc[0][seq_beginning:seq_end]

#     return sequence
                    
# def find_candidate_regions_new(candidate_conserved_region_df0, source_df, segment_size ,continuation_len_thr):
#     """
#     Identify the candidates for conserved regions based on the clusters.
#     In this step, several segments of the same sequence within one cluster are merged if they overlap or are close together.
#     """
#     print("candidate_conserved_region_df0:\n",candidate_conserved_region_df0)

#     # Extract the numbers of the clusters as well as the total number of different sequences
#     cluster_numbers = candidate_conserved_region_df0["cluster"].unique()

#     # Initiate an empty data frame to write the information to
#     columns = ['cluster','genome_id', 'sequence', 'beginning', 'end']
#     candidate_conserved_region_df = pd.DataFrame(columns=columns)

#     for cluster_n in cluster_numbers:
#         # For each cluster, we consider all the segments that are in the cluster
#         current_cluster = candidate_conserved_region_df0[candidate_conserved_region_df0["cluster"]==cluster_n]
#         # We calculate a consensus score that tells us how big a fraction of the different sequences has at least one segment in the cluster
#         cluster_origins = current_cluster['genome_id'].unique()

#         for seq_id in cluster_origins:
#             beginning = current_cluster[current_cluster["genome_id"]==seq_id]['beginning'].values
#             end = current_cluster[current_cluster["genome_id"]==seq_id]['end'].values

#             if(isinstance(beginning,np.ndarray)):
#                 # If more than just one segment per sequence was added, we want to either consider them as one large segment 
#                 # if they are overlapping or close we just adjust beginning and end, otherwise we ?
#                 beginning_sorted = sorted(beginning)
#                 end_sorted = sorted(end)

#                 new_beginning = beginning[0] 
#                 new_end = end[0]
#                 for i in range(1,len(beginning_sorted)):
#                     if (beginning_sorted[i]-beginning_sorted[i-1]-segment_size > continuation_len_thr*segment_size):
#                         # todo: for now i just added both of them but i might want to come up with a better method... one possibility: in this case we want to align only one of the segments, we take the one that is closer to the center
                        
#                         # the two segments are too far apart to be united so we add the first one to the df and continue looking at the next segment
#                         sequence = trim_sequence(source_df,seq_id,new_beginning,new_end,continuation_len_thr,segment_size)
#                         new_row = pd.DataFrame({'cluster' : [cluster_n],'genome_id' : [seq_id],'sequence' : [sequence], 'beginning': [new_beginning], 'end': [new_end]})
#                         candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

#                         # and set the new beginning and end to the new segment
#                         new_beginning = beginning[i] 
#                         new_end = end[i]
#                     else:
#                         # we want to unite the segments, so we overwrite the "new end" and continue looking at the next segment
#                         new_end = end_sorted[i]
#                 # after looking at each segment, we add the last one to our df
#                 sequence = trim_sequence(source_df,seq_id,new_beginning,new_end,continuation_len_thr,segment_size)
#                 new_row = pd.DataFrame({'cluster' : [cluster_n],'genome_id' : [seq_id],'sequence' : [sequence], 'beginning': [new_beginning], 'end': [new_end]})
#                 candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
            
#             else:
#                 # In cases where we only have one segment for a sequence, we sinply add the segment to the df:
#                 # We let the sequence be a bit longer than the candidate region so we can find the actual beginning of the conserved region in the alignment
#                 sequence = trim_sequence(source_df,seq_id,beginning,end,continuation_len_thr,segment_size)
#                 new_row = pd.DataFrame({'cluster' : [cluster_n],'genome_id' : [seq_id],'sequence' : [sequence], 'beginning': [beginning], 'end': [end]})
#                 candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
        
#     print("candidate_conserved_region_df:\n",candidate_conserved_region_df)
#     return candidate_conserved_region_df


def trim_sequence(seq, beginning, end, continuation_len_thr, segment_size):
    """
    Based on beginning and end points, get a new sequence fragment and add a few bases at the beginning and end
    (determined by continuation_len_thr) if in bound of the sequence.
    """
    seq_beginning = max(beginning - round(0.5 * continuation_len_thr * segment_size), 0)
    seq_end = min(end + round(0.5 * continuation_len_thr * segment_size), len(seq))
    return (seq[seq_beginning:seq_end], seq_beginning, seq_end)

def find_candidate_regions_new(candidate_conserved_region_df0, source_df, segment_size, continuation_len_thr):
    """
    Identify the candidates for conserved regions based on the clusters.
    In this step, several segments of the same sequence within one cluster are merged if they overlap or are close together.
    """
    #print("candidate_conserved_region_df0:\n", candidate_conserved_region_df0)

    # Extract the numbers of the clusters as well as the total number of different sequences
    cluster_numbers = candidate_conserved_region_df0["cluster"].unique()

    # Initiate an empty data frame to write the information to
    columns = ['cluster', 'genome_id', 'sequence', 'beginning', 'end']
    candidate_conserved_region_df = pd.DataFrame(columns=columns)

    # Cache sequences and their lengths
    sequence_cache = {}
    for seq_id in source_df["genome_id"].unique():
        seq = source_df.loc[source_df["genome_id"] == seq_id, "sequence"].iloc[0]
        sequence_cache[seq_id] = seq

    for cluster_n in cluster_numbers:
        # For each cluster, we consider all the segments that are in the cluster
        current_cluster = candidate_conserved_region_df0[candidate_conserved_region_df0["cluster"] == cluster_n]
        # We calculate a consensus score that tells us how big a fraction of the different sequences has at least one segment in the cluster
        cluster_origins = current_cluster['genome_id'].unique()

        for seq_id in cluster_origins:
            beginning = current_cluster[current_cluster["genome_id"] == seq_id]['beginning'].values
            end = current_cluster[current_cluster["genome_id"] == seq_id]['end'].values

            if isinstance(beginning, np.ndarray):
                # If more than just one segment per sequence was added, we want to either consider them as one large segment 
                # if they are overlapping or close we just adjust beginning and end, otherwise we ?
                beginning_sorted = np.sort(beginning)
                end_sorted = np.sort(end)

                new_beginning = beginning_sorted[0]
                new_end = end_sorted[0] # why not -1 ?
                for i in range(1, len(beginning_sorted)):
                    if (beginning_sorted[i] - beginning_sorted[i - 1] - segment_size > continuation_len_thr * segment_size):
                        # The two segments are too far apart to be united so we add the first one to the df and continue looking at the next segment
                        seq = sequence_cache[seq_id]
                        (sequence, trimmed_beg, trimmed_end) = trim_sequence(seq, new_beginning, new_end, continuation_len_thr, segment_size)
                        new_row = pd.DataFrame({'cluster': [cluster_n], 'genome_id': [seq_id], 'sequence': [sequence], 'beginning': [trimmed_beg], 'end': [trimmed_end]})
                        candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

                        # And set the new beginning and end to the new segment
                        new_beginning = beginning_sorted[i]
                        new_end = end_sorted[i]
                    else:
                        # We want to unite the segments, so we overwrite the "new end" and continue looking at the next segment
                        new_end = end_sorted[i]
                # After looking at each segment, we add the last one to our df
                seq = sequence_cache[seq_id]
                (sequence, trimmed_beg, trimmed_end) = trim_sequence(seq, new_beginning, new_end, continuation_len_thr, segment_size)
                new_row = pd.DataFrame({'cluster': [cluster_n], 'genome_id': [seq_id], 'sequence': [sequence], 'beginning': [trimmed_beg], 'end': [trimmed_end]})
                candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
            else:
                # In cases where we only have one segment for a sequence, we simply add the segment to the df:
                # We let the sequence be a bit longer than the candidate region so we can find the actual beginning of the conserved region in the alignment
                seq = sequence_cache[seq_id]
                (sequence, trimmed_beg, trimmed_end) = trim_sequence(seq, trimmed_beg, trimmed_end, continuation_len_thr, segment_size)
                new_row = pd.DataFrame({'cluster': [cluster_n], 'genome_id': [seq_id], 'sequence': [sequence], 'beginning': [beginning], 'end': [end]})
                candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

    #print("candidate_conserved_region_df:\n", candidate_conserved_region_df)
    return candidate_conserved_region_df


def evaluate_and_reconsider_clusters(candidate_region_df,segment_df, clusters, cluster_centers, initial_threshold, conservation_threshold, margin=0.2, new_threshold_increase=0.1):
    total_genome_ids = segment_df['genome_id'].nunique()

    #print("df before: \n",candidate_region_df)

    #print("df before: \n",candidate_region_df.shape)

    # Calculate unique genome_ids in each cluster
    clustersIDs = candidate_region_df.groupby('cluster')['genome_id'].nunique()

   #print("total genome ids: ",total_genome_ids,"\nclusterids:\n",clustersIDs)

    # Identify clusters for reconsideration
    #reconsider_clusters = clustersIDs[((clustersIDs / total_genome_ids) < conservation_threshold) & ((clustersIDs / total_genome_ids) >= conservation_threshold - margin)]
    reconsider_clusters = clustersIDs[((clustersIDs / total_genome_ids) < 1) & ((clustersIDs / total_genome_ids) >= conservation_threshold - margin)]
    
    #print("reconsider_clusters\n",reconsider_clusters)

    # Lower threshold
    new_threshold = initial_threshold * (1 + new_threshold_increase)
    

    for cluster_index in reconsider_clusters.index:
        cluster_genome_ids = set(candidate_region_df[candidate_region_df['cluster'] == cluster_index]['genome_id'])
        remaining_genome_ids = set(candidate_region_df['genome_id']) - cluster_genome_ids

        for genome_id in remaining_genome_ids:
            segment_rows = segment_df[segment_df['genome_id'] == genome_id]
            #vectors = np.stack(segment_rows['vector'].values)
            # print("segment rows\n",segment_rows)
            for vector, sequence, begin, end in zip(segment_rows['vector'], segment_rows['sequence'], segment_rows['beginning'], segment_rows['end']):
                #print("vector: ",vector)
                distance = manhattan_distance(cluster_centers[cluster_index], vector)

                if distance <= new_threshold:
                    clusters[cluster_index].append(vector)
                    cluster_centers[cluster_index] = calculate_centroid_mean(np.array(clusters[cluster_index]))
                    new_row = pd.DataFrame({'cluster':cluster_index,'genome_id':genome_id, 'sequence':sequence, 'vector':[vector], 'beginning':begin, 'end':end})
                    candidate_region_df = pd.concat([candidate_region_df, new_row], ignore_index=True)

    # Recalculate unique genome_ids in each cluster

    new_cluster_IDs = candidate_region_df.groupby('cluster')['genome_id'].nunique()
    # print("new_cluster_IDs:\n",new_cluster_IDs)
    # print("total:\n",total_genome_ids)
    # print("divided:\n",(new_cluster_IDs / total_genome_ids))

   # Filter clusters that meet or exceed conservation threshold
    valid_clusters = new_cluster_IDs[(new_cluster_IDs / total_genome_ids) >= conservation_threshold]
    #print("valid_clusters:\n",valid_clusters)

    # Filter dataframe to include only valid clusters
    valid_cluster_indices = valid_clusters.index.tolist()
    filtered_df = candidate_region_df[candidate_region_df['cluster'].isin(valid_cluster_indices)]

    #print("df after: \n",filtered_df.shape)

    return filtered_df, clusters, cluster_centers


def find_clustersNew(segment_df,metric,threshold):
    """
    Cluster the segments based on their metric distance
    """
    only_add_to_best_match = True
    # todo: maybe different way of determining which cluster to add this too if two clusters are equally good.
    # todo: more efficient clustering mechanisms

    # The cluster centers will hold an average vector for each cluster and be used to determine the distance of each segment to each cluster. 
    clusters = []
    cluster_centers = []

    unique_genome_ids = segment_df['genome_id'].unique()

    # Initiate an empty data frame to write the information to
    columns = ['cluster','genome_id', 'sequence', 'vector', 'beginning', 'end']
    candidate_conserved_region_df = pd.DataFrame(columns=columns)
    
    for genome_id in unique_genome_ids:

        # For each of the different genome IDs, look at each vector representing a sequence segment summary vector and either create a new cluster for it
        # or put it into one of the existing clusters
        segment_rows = segment_df[segment_df['genome_id'] == genome_id]
        
        for vector, sequence, begin, end in zip(segment_rows['vector'], segment_rows['sequence'], segment_rows['beginning'], segment_rows['end']):

            if len(cluster_centers) == 0:
                # Create a new cluster for the first segment
                clusters.append([vector])
                # A cluster with only one sequence has that sequence as its center
                cluster_centers.append(vector)
                # The information which cluster we added the segment to is added to the df
                #segment_df.at[rowWithId, 'cluster'] = [0]
                new_row = pd.DataFrame({'cluster':0,'genome_id':genome_id, 'sequence':sequence, 'vector':[vector], 'beginning':begin, 'end':end})
                candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
                
            else:
                if metric == "manhattan":
                        distances = np.array([manhattan_distance(center, vector) for center in cluster_centers])
                else:
                    raise Exception(f"unknown metric {metric} selected for conserved region identification based on quasi alignments")
                
                if only_add_to_best_match:
                    min_distance = distances.min()
                    if min_distance <= threshold:
                        min_dist_cluster_index = distances.argmin()
                        clusters[min_dist_cluster_index].append(vector)
                        cluster_centers[min_dist_cluster_index] = calculate_centroid_mean(np.array(clusters[min_dist_cluster_index]))
                        #segment_df.at[rowWithId, 'cluster'] = [min_dist_cluster_index]
                        new_row = pd.DataFrame({'cluster':min_dist_cluster_index,'genome_id':genome_id, 'sequence':sequence, 'vector':[vector], 'beginning':begin, 'end':end})
                        candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

                    else:
                        clusters.append([vector])
                        cluster_centers.append(vector)
                        #segment_df.at[rowWithId, 'cluster'] = [len(clusters)]
                        new_row = pd.DataFrame({'cluster':len(clusters),'genome_id':genome_id, 'sequence':sequence, 'vector':[vector], 'beginning':begin, 'end':end})
                        candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

                else:
                    if metric == "manhattan":
                        distances = np.array([manhattan_distance(center, vector) for center in cluster_centers])
                    else:
                        raise Exception(f"unknown metric {metric} selected for conserved region identification based on quasi alignments")
                        
                    under_threshold_indices = np.where(distances <= threshold)[0]

                    if under_threshold_indices.size > 0:
                        for cluster_index in under_threshold_indices:
                            clusters[cluster_index].append(vector)
                            cluster_centers[cluster_index] = calculate_centroid_mean(np.array(clusters[cluster_index]))
                            new_row = pd.DataFrame({'cluster':cluster_index,'genome_id':genome_id, 'sequence':sequence, 'vector':[vector], 'beginning':begin, 'end':end})
                            candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)
                    else:
                        clusters.append([vector])
                        cluster_centers.append(vector)
                        #segment_df.at[rowWithId, 'cluster'] = [len(clusters)]
                        new_row = pd.DataFrame({'cluster':len(clusters),'genome_id':genome_id, 'sequence':sequence, 'vector':[vector], 'beginning':begin, 'end':end})
                        candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

    #print("df:\n",candidate_conserved_region_df)

    return candidate_conserved_region_df, clusters, cluster_centers


def align_candidate_regions(candidate_conserved_region_df):
    """
    Align the sequences of the conserved region candidates
    """
    conserved_regions_positions = []
    conserved_regions_dominant = []
    cluster_nums = candidate_conserved_region_df['cluster'].unique()

    for cluster_num in cluster_nums:
        # Get sequences to consider for the current cluster
        sequences_to_consider = candidate_conserved_region_df[candidate_conserved_region_df["cluster"] == cluster_num]

        # Prepare sequence tuples for alignment
        sequence_tuples = [(row["genome_id"], row["sequence"]) for index, row in sequences_to_consider.iterrows()]

        # Write sequences to a temporary FASTA file
        with NamedTemporaryFile(mode='w', delete=False) as temp_file:
            for genome_id, sequence in sequence_tuples:
                temp_file.write(f">{genome_id} alignment_cluster_{cluster_num}_{genome_id}\n{sequence}\n")
            temp_file.flush()
            temp_file_path = temp_file.name

        # Perform sequence alignment
        alignment_path = analyzer.align_sequences(temp_file_path)

        # Read the aligned sequences
        df = read_multifasta(alignment_path)

        # Find conserved regions
        conservedRegion = analyzer.find_conserved_regions_shannon_entropy(df)
        conserved_regions_positions.extend(conservedRegion[0])
        conserved_regions_dominant.extend(conservedRegion[1])

    return conserved_regions_positions, conserved_regions_dominant

def align_candidate_regions_new(candidate_conserved_region_df):
    """
    Align the sequences of the conserved region candidates
    """

    # Initiate an empty data frame to write the information to
    columns = ['cluster',  'genome_id', 'conserved_region', 'beginning', 'end']
    aligned_candidate_conserved_region_df = pd.DataFrame(columns=columns)

    conserved_regions_positions = []
    conserved_regions_dominant = []
    cluster_nums = candidate_conserved_region_df['cluster'].unique()

    for cluster_num in cluster_nums:
        # Get sequences to consider for the current cluster
        sequences_to_consider = candidate_conserved_region_df[candidate_conserved_region_df["cluster"] == cluster_num]

        # Prepare sequence tuples for alignment
        sequence_tuples = [(row["genome_id"], row["sequence"]) for index, row in sequences_to_consider.iterrows()]

        # Write sequences to a temporary FASTA file
        with NamedTemporaryFile(mode='w', delete=False) as temp_file:
            for genome_id, sequence in sequence_tuples:
                temp_file.write(f">{genome_id} alignment_cluster_{cluster_num}_{genome_id}\n{sequence}\n")
            temp_file.flush()
            temp_file_path = temp_file.name

        # Perform sequence alignment
        alignment_path = analyzer.align_sequences(temp_file_path)

        # Read the aligned sequences
        df = read_multifasta(alignment_path)

        # Find conserved regions
        conservedRegion = analyzer.find_conserved_regions_shannon_entropy(df)

        # the position of the identified conserved region in the actual sequence is thus the following:
        seq_pos = {}  # Dictionary to store counts
        for genome_id in df['genome_id'].unique():
            # todo: when using .values[0] here, make sure there actually is only one item in the list
            begin_orig = candidate_conserved_region_df[(candidate_conserved_region_df['cluster'] == cluster_num) & (df['genome_id'] == genome_id)]['beginning'].values[0]
            end_orig = candidate_conserved_region_df[(candidate_conserved_region_df['cluster'] == cluster_num) & (df['genome_id'] == genome_id)]['end'].values[0]
            seq_pos[genome_id] = (begin_orig, end_orig)


        for j in range(0,len(conservedRegion[0])):

            # extract consenus sequence
            consensus_seq = conservedRegion[1][j]

            # extract the beginning and the end of the conserved region 
            (begin_conserved, end_conserved) = conservedRegion[0][j]
            conserved_regions_positions.append((begin_conserved, end_conserved))

            #for i in seq_pos.keys():
            for index, row in df.iterrows():

                current_genome_id = row['genome_id']

                (begin_orig, end_orig) = seq_pos[current_genome_id]

                # Check each sequence for the number of - before and in the conserved region
                sequence = row['sequence']
                offset_begin = sequence[:begin_conserved].count("-")
                offset_end = sequence[:end_conserved].count("-")

                # the beginning and end of the conserved regions in the orignal sequence can be calculated as follows:
                begin = begin_orig + (begin_conserved - offset_begin)
                end = begin_orig + (end_conserved - offset_end)

                new_row = pd.DataFrame({'cluster':cluster_num,'genome_id':current_genome_id, 'conserved_region':[j], 'beginning':[begin], 'end':[end]})
                aligned_candidate_conserved_region_df = pd.concat([aligned_candidate_conserved_region_df, new_row], ignore_index=True)

        conserved_regions_dominant.extend(consensus_seq)

    return aligned_candidate_conserved_region_df, conserved_regions_positions, conserved_regions_dominant


def find_overlapping_conserved_regions(candidate_conserved_region_df: pd.DataFrame, source_df, clusters_to_consider: list[int], clusters_in_final_df: list[int], continuation_len_thr, conservation_thr, number_of_sequences, segment_size):
    """
    Check if the candidate regions overlap or are very close to each other in any of the sequences and if so, combine them.
    """
    for cluster_num in clusters_to_consider:
        # Select the cluster of interest
        cluster_of_interest = candidate_conserved_region_df[candidate_conserved_region_df["cluster"] == cluster_num]
        edited = False

        # Consider clusters that have only one segment per sequence
        # But I am not acutally checking that here at all?
        if len(cluster_of_interest["genome_id"].unique()) == len(cluster_of_interest):
            for cluster_to_compare_to in clusters_to_consider:
                if cluster_num != cluster_to_compare_to:
                    all_seq_distances = []
                    sequences_in_cluster_to_compare_to = candidate_conserved_region_df[candidate_conserved_region_df["cluster"] == cluster_to_compare_to]["genome_id"]

                    if len(sequences_in_cluster_to_compare_to.unique()) == len(sequences_in_cluster_to_compare_to):
                        sequences_in_both_sets = set(cluster_of_interest["genome_id"]).intersection(sequences_in_cluster_to_compare_to)
                        
                        if len(sequences_in_both_sets) >= conservation_thr * number_of_sequences:
                            for seq_id in sequences_in_both_sets:
                                beginning0 = cluster_of_interest[cluster_of_interest["genome_id"] == seq_id]["beginning"].values[0]
                                beginning1 = candidate_conserved_region_df[(candidate_conserved_region_df["cluster"] == cluster_to_compare_to) & (candidate_conserved_region_df["genome_id"] == seq_id)]["beginning"].values[0]
                                distance = abs(beginning1 - beginning0)
                                all_seq_distances.append(distance)

                            average_distance = np.mean(all_seq_distances)
                            if average_distance <= continuation_len_thr:
                                for seq_id in sequences_in_both_sets:
                                    beginning0 = cluster_of_interest[cluster_of_interest["genome_id"] == seq_id]["beginning"].values[0]
                                    beginning1 = candidate_conserved_region_df[(candidate_conserved_region_df["cluster"] == cluster_to_compare_to) & (candidate_conserved_region_df["genome_id"] == seq_id)]["beginning"].values[0]
                                    new_beginning = min(beginning0, beginning1)
                                    end0 = cluster_of_interest[cluster_of_interest["genome_id"] == seq_id]["end"].values[0]
                                    end1 = candidate_conserved_region_df[(candidate_conserved_region_df["cluster"] == cluster_to_compare_to) & (candidate_conserved_region_df["genome_id"] == seq_id)]["end"].values[0]
                                    new_end = max(end0, end1)
                                    (new_sequence, trimmed_beg, trimmed_end) = trim_sequence(source_df, seq_id, new_beginning, new_end, continuation_len_thr, segment_size)
                                    new_cluster_name = candidate_conserved_region_df['cluster'].max() + 1
                                    new_row = pd.DataFrame({'cluster': [new_cluster_name], 'genome_id': [seq_id], 'sequence': [new_sequence], 'beginning': [trimmed_beg], 'end': [trimmed_end]})
                                    candidate_conserved_region_df = pd.concat([candidate_conserved_region_df, new_row], ignore_index=True)

                                clusters_to_consider = np.append(clusters_to_consider[(clusters_to_consider != cluster_num) & (clusters_to_consider != cluster_to_compare_to)], new_cluster_name)
                                edited = True
                                candidate_conserved_region_df, clusters_in_final_df = find_overlapping_conserved_regions(candidate_conserved_region_df, source_df, clusters_to_consider, clusters_in_final_df, continuation_len_thr, conservation_thr, number_of_sequences, segment_size)
                                break

        if not edited:
            clusters_to_consider = clusters_to_consider[clusters_to_consider != cluster_num]
            clusters_in_final_df = np.append(clusters_in_final_df, cluster_num)
        else:
            break

    return candidate_conserved_region_df, clusters_in_final_df

def auto_parameters(segment_df,metric,conservation_threshold, source_df, segment_size, continuation_len_thr, overlap, threshold=None):
    #consider: overlap, conservation thr, segment size
    # we want a sufficient number of candidate clusters but one that is not too big -> start with a guess calculated based on the considered parameters and then refine based on number of candidate clusters 

    #no_sequences = source_df["genome_id"].nunique()
    
    average_seg_length = source_df["sequence"].str.len().mean()/segment_size

    min_cons_regions = max(math.floor(average_seg_length*0.2),1) #((overlap*segment_df.shape[0])/(segment_size))*0.1
    max_cons_regions = max(math.ceil(average_seg_length*0.4), 2)#((overlap*segment_df.shape[0])/(segment_size))*0.5
    # todo: make it possible to edit this

    print("min: ", min_cons_regions,"max: ",max_cons_regions)

    #init_thr = segment_size / 3.5
    init_thr = segment_size / 3
    last_thr = init_thr + 0.2*init_thr

    number_of_candidates = 0

    print(min_cons_regions, max_cons_regions)

    if threshold == None:

        while (number_of_candidates < min_cons_regions) or (number_of_candidates > max_cons_regions):

            candidate_conserved_region_df,clusters,cluster_centers  = find_clustersNew(segment_df,metric,init_thr)

            # we attempt to fill up clusters if they are almost at the conservation threshold and otherwise filter them out
            candidate_conserved_region_df_new,clusters,cluster_centers = evaluate_and_reconsider_clusters(candidate_conserved_region_df,segment_df,clusters,cluster_centers,init_thr,conservation_threshold)

            # calculate stats of the candidate region df
            # number of candidate regions
            no_clusters = len(cluster_centers)
            #print(no_clusters)
            # average length 
            average_no_sequences = sum(len(cluster) for cluster in clusters) / len(clusters)
            #print(average_no_sequences)
            # candidate regions
            number_of_candidates = candidate_conserved_region_df_new["cluster"].nunique()
            print(number_of_candidates)

            # todo: make sure this can not run indefinateley
            if number_of_candidates < min_cons_regions:
                last_thr = init_thr
                init_thr = init_thr + 0.5
                print("increasing thr ",init_thr)
            if number_of_candidates > max_cons_regions:
                dif = last_thr - init_thr
                # if this is positive we are curretnly going down
                last_thr = init_thr
                init_thr = init_thr - max(0.7*abs(dif),0.2)
                print("decreasing thr ",init_thr)

    else:
        candidate_conserved_region_df,clusters,cluster_centers  = find_clustersNew(segment_df,metric,threshold)
        # we attempt to fill up clusters if they are almost at the conservation threshold and otherwise filter them out
        candidate_conserved_region_df_new,clusters,cluster_centers = evaluate_and_reconsider_clusters(candidate_conserved_region_df,segment_df,clusters,cluster_centers,threshold,conservation_threshold)

    candidate_conserved_region_df_new = find_candidate_regions_new(candidate_conserved_region_df_new, source_df, segment_size, continuation_len_thr)

    return candidate_conserved_region_df_new






def quasi_align(source_df,overlap,conservation_threshold,threshold=None,additional_sequence_symbols:list[str]=[],metric:str="manhattan",segment_size:int=300,continuation_len_thr:int=1):
    """
    Find conserved regions without aligning sequences based on a clustering method.
    The sequences are divided into segments of a given length and the number of occurrences of the different possible base triplets is used to 
    determine the distance between two segments based on a given Metric. Similar segments are then clustered to find candidates for conserved regions.

    The parameters used here are:
    - segment_size: The length of segments the sequences are divided into (default: 300)
    - threshold: The maximum allowed metric distance that segments can have and still be clustered together
    - overlap: When dividing the sequence into segments to compare triplet occurrences in, an overlap of segments can be allowed and is given as 1/n to indicate the fraction of segments at which the overlap should begin.
    - conservation_thr: Sets the threshold for the fraction of all sequences that have to be in one cluster in order to consider the cluster as a potential conserved region
    - continuation_len_thr: distance relative to segment size that two segments can have to be considered one continuous segment (default: 1)
    """
    # Todo: different distance for specific symbols? if e.g. n can be a or t then the distance of a n to these should be smaller than to g or c
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

        # Extract information from source df
        sample = source_df.iloc[row]["sequence"]
        genome_id = source_df.iloc[row]["genome_id"]

        for i in range(0, len(sample)-round(segment_size*(1-overlap)), round(segment_size*overlap)):

            # Extract the actual sequence for each segment
            if len(sample) < segment_size:
                print(f"sequence provided for sample {genome_id} is shorter than the minimum segment size set for quasi-alignments ({segment_size}). Skipping sequence.")
            else:
                if (len(sample)-i<segment_size):
                    segBegin = len(sample)-segment_size
                    segEnd = len(sample)
                else:
                    segBegin = i
                    segEnd = i + segment_size

                segment = sample[segBegin:segEnd]

                # Construct a summary vector for triplet occurrences for each of the segments
                vector = []
                for word in possible_words:
                    vector.append(segment.count(word))
                vector = np.array(vector)

                # Add the new information to the new df
                new_row = pd.DataFrame({'genome_id': genome_id, 'beginning' : segBegin, 'end' : segEnd, 'sequence' : segment, 'vector' : [vector], "cluster" : [np.nan]})
                segment_df = pd.concat([segment_df, new_row], ignore_index=True)

    # Find clusters and update the segment_df 
    #candidate_regions_df = find_clustersNew(segment_df,metric,threshold,conservation_threshold, source_df, segment_size, continuation_len_thr)
    candidate_regions_df = auto_parameters(segment_df,metric,conservation_threshold, source_df, segment_size, continuation_len_thr, overlap, threshold)

    # Merge overlapping regions:
    (extended_candidate_regions_df,clusters_in_final_df) = find_overlapping_conserved_regions(candidate_regions_df, source_df, candidate_regions_df["cluster"].unique(), [],continuation_len_thr,conservation_threshold,len(source_df),segment_size)

    candidate_regions_df = extended_candidate_regions_df[extended_candidate_regions_df["cluster"].isin(clusters_in_final_df)]

    # align the candidate regions and calculate their conservation score
    (conserved_regions_positions,conserved_regions_dominant) = align_candidate_regions(candidate_regions_df)

    return(conserved_regions_positions, conserved_regions_dominant)


def quasi_align_new(source_df,overlap,conservation_threshold,threshold=None,additional_sequence_symbols:list[str]=[],metric:str="manhattan",segment_size:int=300,continuation_len_thr:int=1):
    # new: return in a format that contains information on where in the individual sequences the conserved region begins and ends
    """
    Find conserved regions without aligning sequences based on a clustering method.
    The sequences are divided into segments of a given length and the number of occurrences of the different possible base triplets is used to 
    determine the distance between two segments based on a given Metric. Similar segments are then clustered to find candidates for conserved regions.

    The parameters used here are:
    - segment_size: The length of segments the sequences are divided into (default: 300)
    - threshold: The maximum allowed metric distance that segments can have and still be clustered together
    - overlap: When dividing the sequence into segments to compare triplet occurrences in, an overlap of segments can be allowed and is given as 1/n to indicate the fraction of segments at which the overlap should begin.
    - conservation_thr: Sets the threshold for the fraction of all sequences that have to be in one cluster in order to consider the cluster as a potential conserved region
    - continuation_len_thr: distance relative to segment size that two segments can have to be considered one continuous segment (default: 1)
    """
    # Todo: different distance for specific symbols? if e.g. n can be a or t then the distance of a n to these should be smaller than to g or c
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

        # Extract information from source df
        sample = source_df.iloc[row]["sequence"]
        genome_id = source_df.iloc[row]["genome_id"]

        for i in range(0, len(sample)-round(segment_size*(1-overlap)), round(segment_size*overlap)):

            # Extract the actual sequence for each segment
            if len(sample) < segment_size:
                print(f"sequence provided for sample {genome_id} is shorter than the minimum segment size set for quasi-alignments ({segment_size}). Skipping sequence.")
            else:
                if (len(sample)-i<segment_size):
                    segBegin = len(sample)-segment_size
                    segEnd = len(sample)
                else:
                    segBegin = i
                    segEnd = i + segment_size

                segment = sample[segBegin:segEnd]

                # Construct a summary vector for triplet occurrences for each of the segments
                vector = []
                for word in possible_words:
                    vector.append(segment.count(word))
                vector = np.array(vector)

                # Add the new information to the new df
                new_row = pd.DataFrame({'genome_id': genome_id, 'beginning' : segBegin, 'end' : segEnd, 'sequence' : segment, 'vector' : [vector], "cluster" : [np.nan]})
                segment_df = pd.concat([segment_df, new_row], ignore_index=True)

    # Find clusters and update the segment_df 
    #candidate_regions_df = find_clustersNew(segment_df,metric,threshold,conservation_threshold, source_df, segment_size, continuation_len_thr)
    candidate_regions_df = auto_parameters(segment_df,metric,conservation_threshold, source_df, segment_size, continuation_len_thr, overlap, threshold)

    # Merge overlapping regions:
    (extended_candidate_regions_df,clusters_in_final_df) = find_overlapping_conserved_regions(candidate_regions_df, source_df, candidate_regions_df["cluster"].unique(), [],continuation_len_thr,conservation_threshold,len(source_df),segment_size)

    candidate_regions_df = extended_candidate_regions_df[extended_candidate_regions_df["cluster"].isin(clusters_in_final_df)]

    # align the candidate regions and calculate their conservation score
    (conserved_regions_positions_df, conserved_regions_positions ,conserved_regions_dominant) = align_candidate_regions_new(candidate_regions_df)

    # test if the correct original conserved parts were located
    for row in range(0,num_samples):

        # Extract information from source df
        sample = source_df.iloc[row]["sequence"]
        genome_id = source_df.iloc[row]["genome_id"]

        if genome_id in conserved_regions_positions_df["genome_id"].values:
            #begin_orig = candidate_conserved_region_df[(candidate_conserved_region_df['cluster'] == cluster_num) & (df['genome_id'] == genome_id)]['beginning'].values[0]
            beg = conserved_regions_positions_df.loc[(conserved_regions_positions_df["genome_id"] == genome_id) & (conserved_regions_positions_df["conserved_region"] == 0), "beginning"].values[0]
            end = conserved_regions_positions_df.loc[(conserved_regions_positions_df["genome_id"] == genome_id) & (conserved_regions_positions_df["conserved_region"] == 0), "end"].values[0]
            print(f"identified region: {sample[beg:end]}")

    return(conserved_regions_positions_df, conserved_regions_positions,  conserved_regions_dominant)

## figure out weather the regions could encode primers that enclose primers with an appropriate length for species identification
#def identifyEnclosedRegion()

def find_primer_sets(min_len_ampicon, max_len_amplicon):
    # it probably makes more senece to eventually do this step not before but after the design of primers and after finding out which of them can be used together (compatiple melting temperatures etc.)
    # all that i will need to change is however the step where I identify the regions between the primers for the different sequences

    # 1. Identify the beginnging and the end of the potential amplicons
    pass

#############################################################################################

import time

inFileTests = "sequence-5.fasta"

ifDf = read_multifasta(inFileTests)

start_time_quasi_3 = time.time()
(conserved_regions_positions_df, conserved_regions_positions, conserved_regions_dominant_3) = quasi_align_new(ifDf,1/3,0.8,None,[],"manhattan",300,2)
end_time_quasi_3 = time.time()
duration_quasi_3 = end_time_quasi_3 - start_time_quasi_3
# length of conserved region identified this way
numberOfBases_cquasi_3 = 0
for consRegion in conserved_regions_positions:
    numberOfBases_cquasi_3 = numberOfBases_cquasi_3 + (consRegion[1]-consRegion[0])

# statistics
print(f"\n\nRESULT:")
print(f"\nQuasi aligning the sequences took {duration_quasi_3} seconds.\nIdentified conserved regions in alignment: \n{conserved_regions_positions_df} \nConserved region with a total length of {numberOfBases_cquasi_3} were identified")

