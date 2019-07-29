from __future__ import print_function
from __future__ import division

from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

import pandas as pd
import os
import glob
import tqdm
import numpy as np
import urllib

import subprocess
from contextlib import closing

#directory = "/mnt/scratch3/meta_genome/combined_contamenation/sep27_coverage/"
#bam_dir = "/mnt/scratch3/meta_genome/combined_contamenation/sep27_minimap2_bam/"
#cedar_dir = "/mnt/scratch3/meta_genome/combined_contamenation/sep27_cedar_files/"
#
#taxadf = pd.read_table('/mnt/scratch3/meta_genome/combined_contamenation/dataframes/taxadf.tsv')
#df4 = pd.read_table('/mnt/scratch3/meta_genome/combined_contamenation/dataframes/coverage_score_sept27.tsv', index_col=0)
#df4_top = df4.apply(lambda s: s.abs().nlargest(10).index.tolist(), axis=1)

#result = {}
#prefix = "/mnt/scratch3/meta_genome/combined_contamenation/sep27_minimap2_bam/"
#missing_stuff = {}
#sequence_dict = {}
from Bio.SeqUtils import lcc

from multiprocessing import Pool


import re
def slugify(value):
    """
    Normalizes string, converts to lowercase, removes non-alpha characters,
    and converts spaces to hyphens.
    """
    import unicodedata
    # value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = str(re.sub('[^\w\s-]', '', value).strip().lower())
    value = str(re.sub('[-\s]+', '-', value))
    return value


def get_best_hits_from_blast_mono(arg):
    org_name, read_seq = arg[0], arg[1]

    out_dir = 'xml_files'
    xml_file = out_dir + '/' + org_name + '.xml'

    done, success, ResponseData = True, False, ""

    blast_hits = []
    if os.path.exists(xml_file):
        return (done, False, "")

    done = False
    num_tried = 0

    while (num_tried < 10):
        try:
            print("Trying blast for ...", org_name)
            result_handle = NCBIWWW.qblast('blastn', 'nt', read_seq, alignments=5000, hitlist_size=5000)
            print("Done blast for ...", org_name)
            success = True
            break
        except urllib.error.URLError as e: 
            print("Got exception, waiting ..... ")
            ResponseData = e.reason
            time.sleep(1200)

        num_tried += 1
    
    if(success):
        with open(xml_file, "w") as out_handle:
            out_handle.write(result_handle.read())
    else:
        print("failed....")

    return (done, success, ResponseData)


import time
from Bio import SeqIO
from pyfasta import Fasta

illumina_fasta_file = 'illumina_seq.fasta'
rRNA_fasta_file = 'bacterial_rRNA.fasta'

illumina_seq = Fasta(illumina_fasta_file)
rRNA_seq = Fasta(rRNA_fasta_file)

arguments = []
seq_dict = []
for name in rRNA_seq:
    arguments += [(slugify(name), str(rRNA_seq[name]))]  

for name in illumina_seq:
    arguments += [(slugify(name), str(illumina_seq[name]))]  

offset = len(rRNA_seq)

print("Try blasting",len(arguments[offset:]), " sequences")
for i in range(len(arguments[offset:])):
    repeat = True
    arg = arguments[i+offset]

    while(repeat):
        done, success, ResponseData = get_best_hits_from_blast_mono(arg)
        if done:
            repeat = False
        elif success:
            # construct dataframe

            time.sleep(10)
            repeat = False
        else:
            print("Sleeping for a long time....")
            time.sleep(3500)
