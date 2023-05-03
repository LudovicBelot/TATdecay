import os
import subprocess
import shutil
import pandas as pd
import numpy as np
import multiprocessing
from Bio import SeqIO


def create_dir(outdir,**kwargs):

    remove_tmp = kwargs.get("remove_tmp", True)
    #remove_tmp = False
    #in case the tmp folder already exist, remove it to avoid analysis problems
    if remove_tmp == True:
        shutil.rmtree(f"{outdir}/TATdecay/tmp", ignore_errors= True)

    list_dir = ["TATdecay","TATdecay/tmp","TATdecay/results", "TATdecay/tmp/tblastn_db", "TATdecay/tmp/tblastn", "TATdecay/results/1-spot_tblastn", "TATdecay/results/1-spot_tblastn/TA_alignments"]
    list_dir.insert(0,outdir)


    for i in list_dir:
        if i == list_dir[0]:
            try :
                os.mkdir(i)
            except FileExistsError:
                continue
        else:
            try :
                os.mkdir(f"{outdir}/{i}")
            except FileExistsError:
                continue


def core2use(n):
    if int(n) != 0 and int(n) <= multiprocessing.cpu_count():
        return int(n)
    else :
        return int(multiprocessing.cpu_count())


def TA_parser(TAfile,outdir):

    n_seq = 0
    n_index_TA = 0
    d_TA = {}
    for record in SeqIO.parse(TAfile, "fasta"):
        n_seq += 1

        if (n_seq % 2) == 1 :
            n_index_TA += 1
            d_TA[n_index_TA] = {}
            d_TA[n_index_TA]["Toxin_name"] = record.id
            d_TA[n_index_TA]["Toxin_seq"] = record.seq

        elif (n_seq % 2) == 0 :
            d_TA[n_index_TA]["Antitoxin_name"] = record.id
            d_TA[n_index_TA]["Antitoxin_seq"] = record.seq
    
    df_TA = pd.DataFrame.from_dict(d_TA, orient = "index")
    df_TA.to_csv(f"{outdir}/TATdecay/tmp/TAs_index.tsv", sep= "\t", index = False)
    
    return df_TA


def generate_spot_db(df_TA, core_file, replicons_folder, outdir, **kwargs):
    n_cpu = kwargs.get("cpu", 1)

    df_core = pd.read_csv(core_file, sep = "\t")
    df_core = df_core.reset_index(drop = False)
    d_seqIO = {}
    for file in os.listdir(replicons_folder):
        d_seqIO[file.rsplit(".",1)[0]] = SeqIO.index(f"{replicons_folder}/{file}", "fasta")

    list_TA_per_jobs = np.array_split(df_TA.index.tolist(), n_cpu)
    list_TA_per_jobs = [x.tolist() for x in list_TA_per_jobs if x.size != 0]
    list_iterables4multi_cpu = []

    for i in range(0,len(list_TA_per_jobs)):
        list_iterables4multi_cpu.append([list_TA_per_jobs[i], df_TA, df_core, replicons_folder, outdir])
    
    with multiprocessing.Pool(n_cpu) as pool:
        pool.starmap(func=recreate_spot, iterable=list_iterables4multi_cpu)



def recreate_spot(list_index_TA, df_TA, df_core, replicons_folder, outdir):

    d_seqIO = {}
    str_list_genomes = ""
    for file in os.listdir(replicons_folder):
        d_seqIO[file.rsplit(".",1)[0]] = SeqIO.index(f"{replicons_folder}/{file}", "fasta")
        str_list_genomes += f"{file.rsplit('.',1)[0]}\n"
    
    with open(f"{outdir}/TATdecay/tmp/list_genomes.lst", "w") as f:
        f.write(str_list_genomes)

    for TA_index in list_index_TA:
        tmp_left_core = df_TA.loc[TA_index,"left_core_family"]
        tmp_right_core = df_TA.loc[TA_index,"right_core_family"]
        tmp_str = ""

        #now recreating the same interval in each genome
        for genome in d_seqIO.keys():
            new_left_core_coordinates = df_core[(df_core["core_family"] == tmp_left_core) & (df_core["genome_name"] == genome)][["contig", "index","left_coordinate", "right_coordinate"]].values.tolist()[0]
            new_right_core_coordinates = df_core[(df_core["core_family"] == tmp_right_core) & (df_core["genome_name"] == genome)][["contig","index","left_coordinate", "right_coordinate"]].values.tolist()[0]
            

            if new_left_core_coordinates[0] == new_right_core_coordinates[0]: #same contig
                #need to check whether there are no core rearrangement
                if check_core_rearangment(min(new_left_core_coordinates[2:]+new_right_core_coordinates[2:]), max(new_left_core_coordinates[2:]+new_right_core_coordinates[2:]), new_left_core_coordinates[0],df_core) == True:
                    tmp_str += f">{new_left_core_coordinates[0]}_{min(new_left_core_coordinates[2:]+new_right_core_coordinates[2:])}-{max(new_left_core_coordinates[2:]+new_right_core_coordinates[2:])}\n"
                    tmp_str += f"{d_seqIO[genome][new_left_core_coordinates[0]].seq[min(new_left_core_coordinates[2:]+new_right_core_coordinates[2:]):max(new_left_core_coordinates[2:]+new_right_core_coordinates[2:])]}\n"

            elif new_left_core_coordinates[0] != new_right_core_coordinates[0]:
                tmp_str += get_best_combination(new_left_core_coordinates, new_right_core_coordinates, df_core, d_seqIO, genome)

        with open(f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta", "w") as f:
            f.write(tmp_str)
        
        subprocess.call(f"makeblastdb -in {outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta -dbtype nucl -parse_seqids -input_type fasta -out {outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots", shell = True)


def check_core_rearangment(tmp_interval_leftc, tmp_interval_rightc, contig, df_core):

    df_is_rearrangment = df_core[(df_core["contig"] == contig) & (df_core["left_coordinate"] >= tmp_interval_leftc) & (df_core["right_coordinate"] <= tmp_interval_rightc)]

    if len(df_is_rearrangment) == 2:
        return True
    else:
        return False




def get_best_combination(core1, core2, df_core, d_seqIO, genome):
    
    df_tmp_core1 = df_core[df_core["contig"] == core1[0]].sort_values(by = "left_coordinate").reset_index(drop=True)
    df_tmp_core2 = df_core[df_core["contig"] == core2[0]].sort_values(by = "left_coordinate").reset_index(drop=True)
    tmp_str = ""

    if (df_tmp_core1.loc[0,"index"] == core1[1] or df_tmp_core1.loc[df_tmp_core1.index[-1],"index"] == core1[1]) and (df_tmp_core2.loc[0,"index"] == core2[1] or df_tmp_core2.loc[df_tmp_core2.index[-1],"index"] == core2[1]):
        
        if df_tmp_core1.loc[0,"index"] == core1[1]:
            tmp_str += f">{core1[0]}_1-{max(core1[2],core1[3])}\n"
            tmp_str += f"{d_seqIO[genome][core1[0]].seq[:max(core1[2],core1[3])]}\n"
        elif df_tmp_core1.loc[df_tmp_core1.index[-1],"index"] == core1[1]:
            tmp_str += f">{core1[0]}_{min(core1[2],core1[3])}-{len(d_seqIO[genome][core1[0]].seq)}\n"
            tmp_str += f"{d_seqIO[genome][core1[0]].seq[min(core1[2],core1[3]):]}\n"

        if df_tmp_core2.loc[0,"index"] == core2[1]:
            tmp_str += f">{core2[0]}_1-{max(core2[2],core2[3])}\n"
            tmp_str += f"{d_seqIO[genome][core2[0]].seq[:max(core2[2],core2[3])]}\n"
        elif df_tmp_core2.loc[df_tmp_core2.index[-1],"index"] == core2[1]:
            tmp_str += f">{core2[0]}_{min(core2[2],core2[3])}-{len(d_seqIO[genome][core2[0]].seq)}\n"
            tmp_str += f"{d_seqIO[genome][core2[0]].seq[min(core2[2],core2[3]):]}\n"
    
    return tmp_str




