import os
import subprocess
import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbitblastnCommandline
import multiprocessing


def TA_ref_locations(df_TA, TAfile, replicons_folder, core_file, outdir, **kwargs):
    ref = kwargs.get("ref", "First")

    list_replicons = [x.rsplit(".",1)[0] for x in os.listdir(replicons_folder)]
    replicons_extension = os.listdir(replicons_folder)[0].rsplit(".",1)[1]
    if ref == "First":
        ref_genome = list_replicons[0]

    else :
        if ref in list_replicons :
            ref_genome = ref
        else :
            print(f"{ref} not in {replicons_folder}")
            raise Exception

    # doing a tblastn for each TAs to get their correct location in the reference genome
    subprocess.call(f"makeblastdb -in {replicons_folder}/{ref_genome}.{replicons_extension} -dbtype nucl -input_type fasta -out {outdir}/TATdecay/tmp/tblastn_db/ref_replicondb", shell = True)
    ref_tblastn_cline = NcbitblastnCommandline(query = TAfile, db = f"{outdir}/TATdecay/tmp/tblastn_db/ref_replicondb", evalue = 0.001,
                                               outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                                               out = f"{outdir}/TATdecay/tmp/tblastn/ref_tblastn.tsv")
    ref_tblastn_cline()

    df_tblastn_ref = pd.read_csv(f"{outdir}/TATdecay/tmp/tblastn/ref_tblastn.tsv", sep = "\t", names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
    df_tblastn_ref = df_tblastn_ref.drop_duplicates(subset = "qseqid")

    
    #updating the df_TA with the coordinates of each components in the ref genomes
    df_TA = df_TA.merge(df_tblastn_ref[df_tblastn_ref["qseqid"].isin(df_TA["Toxin_name"].tolist())][["qseqid","sstart","send"]], left_on = "Toxin_name", right_on = "qseqid")
    df_TA.drop(columns ="qseqid", inplace = True)
    df_TA.rename(columns = {"sstart":"tox_start", "send":"tox_end"}, inplace = True)
    df_TA = df_TA.merge(df_tblastn_ref[df_tblastn_ref["qseqid"].isin(df_TA["Antitoxin_name"].tolist())][["qseqid","sstart","send"]], left_on = "Antitoxin_name", right_on = "qseqid")
    df_TA.drop(columns ="qseqid", inplace = True)
    df_TA.rename(columns = {"sstart":"antitox_start", "send":"antitox_end"}, inplace = True)

    #now we get the closest core families for each TA
    df_core = pd.read_csv(core_file, sep = "\t")
    df_core = df_core[df_core["genome_name"] == ref_genome]
    df_TA[["left_core_family","ref_interval_leftc"]] = df_TA.apply(lambda x: get_closest_core(x,df_core,"left"), axis = 1)
    df_TA[["right_core_family","ref_interval_rightc"]] = df_TA.apply(lambda x: get_closest_core(x,df_core,"right"), axis = 1)
    
    df_TA.to_csv(f"{outdir}/TATdecay/tmp/TAs_index.tsv", sep= "\t", index = False)
    return df_TA




def get_closest_core(TA_row, df_core, which):

    if which == "left":
        df_tmp_core = df_core[df_core["left_coordinate"] <= min(TA_row[["tox_start","tox_end","antitox_start","antitox_end"]])]
        df_tmp_core = df_tmp_core.sort_values(by = "left_coordinate", ascending = False).reset_index(drop = True)
        return pd.Series([df_tmp_core.loc[0,"core_family"], int(df_tmp_core.loc[0,"left_coordinate"])])
    
    elif which == "right":
        df_tmp_core = df_core[df_core["right_coordinate"] >= min(TA_row[["tox_start","tox_end","antitox_start","antitox_end"]])]
        df_tmp_core = df_tmp_core.sort_values(by = "right_coordinate", ascending = True).reset_index(drop = True)
        return pd.Series([df_tmp_core.loc[0,"core_family"], int(df_tmp_core.loc[0,"right_coordinate"])])



def all_genomes_search(df_TA, outdir, **kwargs):
    n_cpu = kwargs.get("cpu",1)

    list_TA_per_jobs = np.array_split(df_TA.index.tolist(), n_cpu)
    list_TA_per_jobs = [x.tolist() for x in list_TA_per_jobs if x.size != 0]
    list_iterables4multi_cpu = []

    for i in range(0,len(list_TA_per_jobs)):
        list_iterables4multi_cpu.append([list_TA_per_jobs[i], df_TA, outdir])

    with multiprocessing.Pool(n_cpu) as pool:
        pool.starmap(func=multi_cpu_tblastn, iterable=list_iterables4multi_cpu)




def multi_cpu_tblastn(list_TA_index, df_TA, outdir):

    for TA_index in list_TA_index:
        tmp_str = ""
        tmp_str += f">{df_TA.loc[TA_index,'Toxin_name']}\n{df_TA.loc[TA_index,'Toxin_seq']}\n"
        tmp_str += f">{df_TA.loc[TA_index,'Antitoxin_name']}\n{df_TA.loc[TA_index,'Antitoxin_seq']}\n"
        with open(f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_TA_query.fna", "w") as f:
            f.write(tmp_str)

        tblastn_all_genomes_cline = NcbitblastnCommandline(query = f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_TA_query.fna", 
                                db = f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots", evalue = 0.1,
                                outfmt = "6 qseqid sseqid pident length qstart qend sstart send sframe evalue qseq sseq",
                                out = f"{outdir}/TATdecay/tmp/tblastn/{TA_index}_all_genomes_same_spot.tsv")
        tblastn_all_genomes_cline()
        

