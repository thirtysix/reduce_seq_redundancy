#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 2.7.0 ###########################################################
# Libraries ####################################################################
import os
from Bio import SeqIO
from operator import itemgetter
import re
import time

################################################################################
# Functions ####################################################################
################################################################################    

def clustalo(input_seqs_filename, dist_mat_outfilename, alignment_outfilename):
    """
    Requires installation of clustal omega.  Creates a distance matrix from
    all entries in the input_seqs_filename file.  FASTA descriptions limited
    to 128 characters or less.
    """

    clustalo_args = ["clustalo -i", input_seqs_filename, "--full", "--distmat-out=" + dist_mat_outfilename, ">", alignment_outfilename]
    clustalo_str = " ".join(clustalo_args)
    os.system(clustalo_str)
    

def parse_matrix(dist_mat_outfilename):
    """
    Parse into a dictionary the distance matrix file output by Clustal Omega.
    """

    with open(dist_mat_outfilename, 'r') as dist_mat_file:
        line_count = 0
        distances_dict = {}
        readlines = dist_mat_file.readlines()
        entries_count = int(readlines[0].strip())
        

        # iterate through each line of the matrix file
        for line in readlines[1:]:
            distances = line.split(" ")
            distances_dict[line_count] = []

            distances = line[128:]
            distances = distances.split(" ")
            distances = [float(x) for x in distances]

            if len(distances) == entries_count:
                distances_dict[line_count] = distances
            else:
                print "WAAAAAAAT"

            line_count += 1

    return distances_dict


def sequence_list(input_seqs_filename):
    """
    Import sequences from input sequence file.
    """

    with open(input_seqs_filename, 'r') as input_seqs_file:
        input_seqs = list(SeqIO.parse(input_seqs_file, 'fasta'))

    return input_seqs


def cluster_similar(distances_dict, similarity_thresh):
    """
    Identify clusters of similar sequences.
    """

    similar_indexes_dict = {}

    # iterate threugh distances matrices, identify indexes which
    # are similar to the current index
    for index, distances in distances_dict.iteritems():
        distance_count = 0
        similar_indexes = []

        for distance in distances:
            if distance_count != index:
                if distance <= similarity_thresh:
                    similar_indexes.append(distance_count)

            distance_count += 1

        similar_indexes.append(index)
        similar_indexes_dict[index] = similar_indexes            

    return similar_indexes_dict


def prune_clusters(input_seqs, similar_indexes_dict, distances_dict):
    """
    Retain the sequence that is most similar to the others, ignore the rest.
    """
    cleaned_output_list = []
    to_ignore_list = []
    for index, similar_indexes in similar_indexes_dict.iteritems():
        similar_indexes = [x for x in similar_indexes if x not in to_ignore_list]
        similar_indexes = [x for x in similar_indexes if x not in cleaned_output_list]

        # ignore entries which have already been added, or which are similar to a previously added entry
        if index not in to_ignore_list and index not in cleaned_output_list:
            node, to_ignore = identify_most_similar(similar_indexes, distances_dict)
            cleaned_output_list.append(node[0])
            to_ignore_list += to_ignore

    return cleaned_output_list, to_ignore_list      
        

def identify_most_similar(similar_indexes, distances_dict):
    """
    From a group of indexes, identify the sequence most like the others.
    """

    distances_list = []
    for i in range(0, len(similar_indexes)):
        current_index = similar_indexes[i]
        other_indexes = similar_indexes[:i] + similar_indexes[i+1:]

        distances_sum = sum([distances_dict[current_index][other_index] for other_index in other_indexes])
        distances_list.append((current_index, distances_sum))

    distances_list_sorted = sorted(distances_list, key=itemgetter(1))

    node = distances_list_sorted[0]
    to_ignore = [x[0] for x in distances_list_sorted[1:]]

    return node, to_ignore
    

def write_cleaned_entries(cleaned_output_list, unique_entries_outfilename, input_seqs):
    """
    Write to fasta file the unique entries.
    """

    unique_entries = [input_seqs[x] for x in cleaned_output_list]
    with open(unique_entries_outfilename, 'w') as outfile:
        SeqIO.write(unique_entries, outfile, 'fasta')


################################################################################
# Initiating Variables #########################################################
################################################################################
# enter name of starting sequence file
input_seqs_filename = "./all.results.fasta.unique_seq"


dist_mat_outfilename = "./distance_matrix.mat"
alignment_outfilename = "./aligned_seqs.fasta"
unique_entries_outfilename = input_seqs_filename + ".unique_nodes"
similarity_thresh = 0.10


################################################################################
# Execution ####################################################################
################################################################################
start_time = time.time()
clustalo(input_seqs_filename, dist_mat_outfilename, alignment_outfilename)
distances_dict = parse_matrix(dist_mat_outfilename)
input_seqs = sequence_list(input_seqs_filename)
similar_indexes_dict = cluster_similar(distances_dict, similarity_thresh)
cleaned_output_list, to_ignore_list = prune_clusters(input_seqs, similar_indexes_dict, distances_dict)
write_cleaned_entries(cleaned_output_list, unique_entries_outfilename, input_seqs)
end_time = time.time()
print "total time", end_time-start_time, "seconds"




