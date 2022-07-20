#!/usr/bin/python3

import sys, getopt, gzip, msprime, stdpopsim, tskit, demes
from demes import convert

#Initialise arguments
def main(argv):
    global inputfile
    global chromosome
    global outputfile
    global numind

    inputfile = ''
    chromosome = ''
    outputfile = ''
    numind = 1000

    try:
        opts, args = getopt.getopt(argv,"hi:c:o:n:",["inFile=","chr=","outFile=","numInd="])
    except getopt.GetoptError:
        print('gen_trees.py -i <inputfile> -c <chromosome> -o <outputfile> -n <numind>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('gen_trees.py -i <inputfile> -c <chromosome> -o <outputfile> -n <numind>')
            sys.exit()
        elif opt in ("-i", "--inFile"):
            inputfile = arg
        elif opt in ("-c", "--chromosome"):
            chromosome = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg 
        elif opt in ("-n", "--numInd"):
            numind = int(arg)
    print('Input file is', inputfile)
    print('Chromosome is', chromosome)
    print('Number of individuals', numind)
    print('Output file is', outputfile)
   
#Write out vcf
def write_vcf(treeseq, chromosome, outputfile):
    sample_nodes = list()
    pop_nodes = list()
    for node in treeseq.nodes():
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
    n_dip_indv = int(len(sample_nodes) / 2)
    indv_names = [f"{pop_nodes[2*i]}{sample_nodes[2*i]}-{pop_nodes[2*i+1]}{sample_nodes[2*i+1]}" for i in range(n_dip_indv)]
    with gzip.open(str(outputfile + ".vcf.gz"), "wt") as vcf_file:
        treeseq.write_vcf(vcf_file, position_transform="legacy", contig_id=chromosome, individual_names=indv_names)

#Generate treesequences
def gen_trees(inputfile, numind, seed, hapmap, outputfile):
    multi_graph = demes.load(inputfile)
    ms_model = convert.to_msprime(multi_graph)
    ms_demography = msprime.Demography.from_old_style(population_configurations=ms_model[0], demographic_events=ms_model[1], migration_matrix=ms_model[2])

    gcr = 2e-7
    gctl = 125
    
    ts = msprime.sim_ancestry(samples={"Nama": numind,
                                   "SAC": numind,
                                   "GBR": numind/5,
                                   "MSL": numind/5,
                                   "EP": numind/5,
                                   "EAS": numind/5,
                                   "SAS": numind/5},
                          demography=ms_demography, 
                          recombination_rate=hapmap,
                          gene_conversion_rate=gcr,
                          gene_conversion_tract_length=gctl,
                          random_seed=seed,
                          record_migrations=True)
    ts.dump(str(outputfile + ".tree"))
    return(ts)

#Superimpose mutations on treesequences using default Jukes-Canter
def gen_mutations(treeseq, mu, seed, outputfile):
    ts_mu = msprime.sim_mutations(treeseq, rate=mu, random_seed=seed)
    ts_mu.dump(str(outputfile + "_mu.tree"))
    return(ts_mu)

if __name__ == "__main__":
    main(sys.argv[1:])
 
species = stdpopsim.get_species("HomSap")
contig = species.get_contig(str("chr" + str(chromosome)))

hapmap = species.get_genetic_map("HapMapII_GRCh37").get_chromosome_map(str("chr" + str(chromosome))).map
seed = 1111
mu = contig.mutation_rate

ts = gen_trees(inputfile, numind, seed, hapmap, outputfile)

ts_mu = gen_mutations(ts, mu, seed, outputfile)

write_vcf(ts_mu, chromosome, outputfile)
    
