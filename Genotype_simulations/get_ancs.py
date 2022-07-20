#!/usr/bin/python3

import tskit, sys, getopt
import numpy as np
import pandas as pd
import time as t
from collections import defaultdict

#Initialise arguments
def main(argv):
    global tree_file
    global sites_file
    global cutoff
    global outputfile
    global begin
    global end

    tree_file = ''
    sites_file = ''
    cutoff = ''
    outputfile = ''
    begin = ''
    end = ''

    try:
        opts, args = getopt.getopt(argv,"hi:s:c:o:b:e:",["treeFile=", "sitesFile=","cutoff=","outFile=", "begin=", "end="])
    except getopt.GetoptError:
        print('gen_trees.py -i <tree_file> -s <sites_file> -c <cutoff> -o <outputfile> -b <begin> -e <end>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('gen_trees.py -i <tree_file> -c <cutoff> -s <sites_file> -o <outputfile> -b <begin> -e <end>')
            sys.exit()
        elif opt in ("-i", "--treeFile"):
            tree_file = str(arg)
        elif opt in ("-s", "--sitesFile"):
            sites_file = str(arg)
        elif opt in ("-c", "--cutoff"):
            cutoff = int(arg)
        elif opt in ("-o", "--ofile"):
            outputfile = arg 
        elif opt in ("-b", "--begin"):
            begin = int(arg)
        elif opt in ("-e", "--end"):
            end = int(arg)
    print('Input file is', tree_file)
    print('Sites to keep', sites_file)
    print('Time cutoff is', cutoff)
    print('First individual', begin)
    print('Last individual', end)
    print('Output file is', outputfile)

def find_local_ancestry(samples, time, ts, mig_int_tree):
    t0 = t.time()
    ancestor_before_timex_of_tree = dict()
    pop_at_time_of_parent = dict()
    pop_at_time_of_tree = dict()
    intervals_of_tree = dict()
    merged_segments_from_pop = dict()
    for sample in samples:
        ancestor_before_timex_of_tree[sample] = dict()
        pop_at_time_of_tree[sample] = dict()
        intervals_of_tree[sample] = dict()
        merged_segments_from_pop[sample] = dict()
    for tree in ts.trees():
        if(tree.num_sites > 0):
            for sample in samples:
                target = sample
                node_time = tree.time(target)
                parent_node = tree.parent(target)
                if parent_node != tskit.NULL:
                    parent_time = tree.time(tree.parent(target))
                else:
                    parent_time = time+1 
                while parent_time < time:
                    node_time = parent_time
                    target = tree.parent(target)
                    parent_node = tree.parent(target)
                    if parent_node != tskit.NULL:
                        parent_time = tree.time(tree.parent(target))
                    else:
                        parent_time = time+1 
                parent_node = target
                if parent_node in mig_int_tree:
                    overlapping_migrations = list(filter(lambda x: x.left <= tree.interval[0] and x.right >= tree.interval[0], mig_int_tree[parent_node]))
                else:
                    overlapping_migrations = []
                pop_at_time_of_parent[parent_node] = tree.population(parent_node)
                if len(overlapping_migrations) > 0:
                    overlapping_migrations = sorted(overlapping_migrations, 
                                                    key = lambda x : x.time)
                    last_mig = overlapping_migrations.pop()
                    pop_at_time_of_parent[parent_node] = last_mig.dest
                for site in tree.sites():
                    pop_at_time_of_tree[sample][int(site.position)] = pop_at_time_of_parent[parent_node]
    t5 = t.time()
    print("done in ", t5-t0, "seconds.")
    return(pop_at_time_of_tree)

if __name__ == "__main__":
    main(sys.argv[1:])
 
ts_mu = tskit.load(tree_file)
ex = list(pd.read_csv(sites_file, sep="\t", header=None)[1])

ex_ids = []
for i in ts_mu.sites():
    if(int(i.position) not in ex):
        ex_ids.append(i.id)
ts_trim = ts_mu.delete_sites(ex_ids)

sample_nodes = list()
pop_nodes = list()
for node in ts_trim.nodes():
    if node.time == 0 :
        sample_nodes.append(node.id)
        if node.population == 4:
            pop_nodes.append("Nama")
        if node.population == 6:
            pop_nodes.append("MSL")
        if node.population == 7:
            pop_nodes.append("GBR")
        if node.population == 8:
            pop_nodes.append("EP")
        if node.population == 10:
            pop_nodes.append("EAS")
        if node.population == 11:
            pop_nodes.append("SAS")
        if node.population == 12:
            pop_nodes.append("SAC")

time = 100
mig_int_tree = dict()
for migration in ts_mu.migrations():
    if migration.time < time:
        if migration.node not in mig_int_tree:
            mig_int_tree[migration.node] = [migration]
        else:
            mig_int_tree[migration.node].append(migration)

LA = find_local_ancestry(sample_nodes[begin:end], time, ts_trim, mig_int_tree)

frame = pd.DataFrame(LA)
frame.to_csv(outputfile)
