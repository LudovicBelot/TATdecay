import pandas as pd
import numpy as np
import multiprocessing
import itertools
from Bio import SeqIO

def is_both_components(df_TA, outdir, **kwargs):
    n_cpu = kwargs.get("cpu",1)

    list_genomes = []
    with open(f"{outdir}/TATdecay/tmp/list_genomes.lst", "r") as f:
        list_genomes = [line.strip() for line in f if line.strip() != ""]

    list_TA_per_jobs = np.array_split(df_TA.index.tolist(), n_cpu)
    list_TA_per_jobs = [x.tolist() for x in list_TA_per_jobs if x.size != 0]
    list_iterables4multi_cpu = []

    for i in range(0,len(list_TA_per_jobs)):
        list_iterables4multi_cpu.append([list_TA_per_jobs[i], df_TA, list_genomes, outdir])

    with multiprocessing.Pool(n_cpu) as pool:
        all_d_res = pool.starmap(func=multi_is_both_components, iterable=list_iterables4multi_cpu)

    all_d_res = dict(itertools.chain(*map(dict.items, all_d_res)))
    df_all_res = pd.DataFrame.from_dict(all_d_res, orient = "index")
    df_all_res.to_csv(f"{outdir}/TATdecay/results/1-spot_tblastn/all_TAs_decays.tsv", sep = "\t")

    return df_all_res


def multi_is_both_components(list_index_TA, df_TA, list_genomes, outdir):

    d_res = {}
    for TA_index in list_index_TA:
        df_tblastn = pd.read_csv(f"{outdir}/TATdecay/tmp/tblastn/{TA_index}_all_genomes_same_spot.tsv", sep = "\t",
                                 names = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "sframe", "evalue", "qseq", "sseq"])

        #adding a new column which is the name of the genome
        df_tblastn["genome_name"] = df_tblastn.apply(lambda x: x["sseqid"].rsplit(".",1)[0], axis =1)
        df_tblastn.to_csv(f"{outdir}/TATdecay/tmp/tblastn/{TA_index}_all_genomes_same_spot.tsv", sep = "\t", index = False)

        #now checking for each genome whether 0,1 or 2 components were found for each TAs
        for genome in list_genomes:
            d_res[f"{TA_index}_{genome}"] = {"genome_name": genome,
                               "toxin_name": df_TA.loc[TA_index,"Toxin_name"],
                               "antitoxin_name": df_TA.loc[TA_index,"Antitoxin_name"],                       
                               "tblastn_toxin": "No",
                               "tblastn_antitoxin": "No",
                               "min_distance_(bp)":""
                               }

            if df_tblastn[(df_tblastn["genome_name"] == genome) & (df_tblastn["qseqid"] == df_TA.loc[TA_index,"Toxin_name"])].empty == False:
                d_res[f"{TA_index}_{genome}"]["tblastn_toxin"] = "Yes"
            if df_tblastn[(df_tblastn["genome_name"] == genome) & (df_tblastn["qseqid"] == df_TA.loc[TA_index,"Antitoxin_name"])].empty == False:
                d_res[f"{TA_index}_{genome}"]["tblastn_antitoxin"] = "Yes"
            
            if d_res[f"{TA_index}_{genome}"]["tblastn_toxin"] == "Yes" and d_res[f"{TA_index}_{genome}"]["tblastn_antitoxin"] == "Yes":
                d_res[f"{TA_index}_{genome}"]["min_distance_(bp)"] = get_TA_min_distance(df_tblastn[(df_tblastn["genome_name"] == genome) & (df_tblastn["qseqid"] == df_TA.loc[TA_index,"Toxin_name"])]["sstart"],
                                                                           df_tblastn[(df_tblastn["genome_name"] == genome) & (df_tblastn["qseqid"] == df_TA.loc[TA_index,"Antitoxin_name"])]["sstart"]
                                                                            )

    return d_res


def get_TA_min_distance(list_coordinates_1,list_coordinates_2):

    min_value = 100000000000
    for value1 in list_coordinates_1:
        for value2 in list_coordinates_2:
            if abs(value1-value2) < min_value:
                min_value = abs(value1-value2)
    
    return min_value


def check_component_integrity(df_TA, df_tblastn, outdir, **kwargs):
    n_cpu = kwargs.get("cpu",1)
    ref_genome = kwargs.get("ref", "First")

    list_genomes = []
    with open(f"{outdir}/TATdecay/tmp/list_genomes.lst", "r") as f:
        list_genomes = [line.strip() for line in f if line.strip() != ""]
    
    if ref_genome == "First":
        ref_genome = list_genomes[0]

    list_TA_per_jobs = np.array_split(df_TA.index.tolist(), n_cpu)
    list_TA_per_jobs = [x.tolist() for x in list_TA_per_jobs if x.size != 0]
    list_iterables4multi_cpu = []

    for i in range(0,len(list_TA_per_jobs)):
        list_iterables4multi_cpu.append([list_TA_per_jobs[i], list_genomes, ref_genome, df_TA, outdir])

    with multiprocessing.Pool(n_cpu) as pool:
        all_d_res = pool.starmap(func=multi_check_component_integrity, iterable=list_iterables4multi_cpu)

    all_d_res = dict(itertools.chain(*map(dict.items, all_d_res)))
    df_all_res = pd.DataFrame.from_dict(all_d_res, orient ="index")

    df_all_res = df_tblastn.merge(df_all_res, left_index = True, right_index =True)
    df_all_res.to_csv(f"{outdir}/TATdecay/results/1-spot_tblastn/all_TAs_decays.tsv", sep = "\t")
    return df_all_res


def multi_check_component_integrity(list_TA_index, list_genomes, reference_genome, df_TA, outdir):

    d_res = {}
    for TA_index in list_TA_index:
        ref_toxin_seq = df_TA.loc[TA_index,"Toxin_seq"]
        ref_antitoxin_seq = df_TA.loc[TA_index,"Antitoxin_seq"]
        df_tblastn_results = pd.read_csv(f"{outdir}/TATdecay/tmp/tblastn/{TA_index}_all_genomes_same_spot.tsv", sep = "\t")
        str_tox_ali = f"#Toxin: {df_TA.loc[TA_index,'Toxin_name']}\n{reference_genome}\t\t {ref_toxin_seq}\n"
        str_antitox_ali = f"#Antitoxin: {df_TA.loc[TA_index,'Antitoxin_name']}\n{reference_genome}\t\t {ref_antitoxin_seq}\n"

        for genome in list_genomes:
            d_res[f"{TA_index}_{genome}"] = {"tox_cov":"",
                                            "tox_identity":"",
                                            "tox_in_CDS":"",
                                            "tox_early_stop":"",
                                            "tox_frameshift":"",
                                            "antitox_cov":"",
                                            "antitox_in_CDS":"",
                                            "antitox_identity":"",
                                            "antitox_early_stop":"",
                                            "antitox_frameshift":"",
                                            }

            if df_tblastn_results[(df_tblastn_results["qseqid"] == df_TA.loc[TA_index,"Toxin_name"]) & (df_tblastn_results["genome_name"] == genome)].empty == False:
                df_tmp = df_tblastn_results[(df_tblastn_results["qseqid"] == df_TA.loc[TA_index,"Toxin_name"]) & (df_tblastn_results["genome_name"] == genome)].reset_index(drop = True)
                # only one hit
                if len(df_tmp) == 1:
                    
                    d_res[f"{TA_index}_{genome}"]["tox_cov"] = f"{len(df_tmp.loc[0,'sseq'])}/{len(ref_toxin_seq)}"
                    d_res[f"{TA_index}_{genome}"]["tox_identity"] = df_tmp.loc[0,'pident']
                    d_res[f"{TA_index}_{genome}"]["tox_early_stop"] = check_early_stop(df_tmp.loc[0,'sseq'])
                    str_tox_ali += f"{genome}\t\t{' '*int(df_tmp.loc[0,'qstart'])}{df_tmp.loc[0,'sseq']}\n"
                    d_res[f"{TA_index}_{genome}"]["tox_in_CDS"] = check_is_part_of_CDS(df_tmp.iloc[0], len(ref_toxin_seq),f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta")

                elif len(df_tmp) > 1:
                    n_tmp_index = 0
                    for row in df_tmp.iterrows():
                        n_tmp_index += 1
                        if n_tmp_index == 1 :
                            tmp_cov = f"{len(row[1]['sseq'])}/{len(ref_toxin_seq)}"
                            tmp_id = row[1]['pident']
                            tmp_early_stop = check_early_stop(row[1]['sseq'])
                            str_tox_ali += f"{genome}\t\t{' '*int(row[1]['qstart'])}{row[1]['sseq']}\n"
                            tmp_frame = row[1]['sframe']
                            tmp_qstart = row[1]['qstart']
                            tmp_qend = row[1]['qend']
                            tmp_tox_isCDS = [check_is_part_of_CDS(row[1], len(ref_toxin_seq),f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta")]

                        elif n_tmp_index > 1 :
                            #keep only if this is another part of the toxin and it is not far from the best hit (and same strand but not frame)
                            if (tmp_frame in [-1, -2, -3] and row[1]["sframe"] in [-1, -2, -3]) or (tmp_frame in [1, 2, 3] and row[1]["sframe"] in [1, 2, 3]):
                                same_strand = True
                            else :
                                same_strand = False

                            if (row[1]['qstart'] < tmp_qstart or row[1]['qend'] > tmp_qend) and same_strand == True:
                                tmp_cov = f"{tmp_cov.split('/',1)[0]}-{len(row[1]['sseq'])}/{tmp_cov.split('/',1)[1]}"
                                tmp_id = f"{tmp_id}-{row[1]['pident']}"
                                if tmp_early_stop == "No":
                                    tmp_early_stop = check_early_stop(row[1]['sseq'])
                                str_tox_ali += f"{' '*len(genome)}\t\t{' '*int(row[1]['qstart'])}{row[1]['sseq']}\n"

                                tmp_tox_isCDS.append(check_is_part_of_CDS(row[1], len(ref_toxin_seq),f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta"))

                                if tmp_frame != row[1]["sframe"]:
                                    d_res[f"{TA_index}_{genome}"]["tox_frameshift"] = "Yes"
                            
                    tmp_tox_isCDS = (',').join(tmp_tox_isCDS)


                    d_res[f"{TA_index}_{genome}"]["tox_cov"] = tmp_cov
                    d_res[f"{TA_index}_{genome}"]["tox_identity"] = tmp_id
                    d_res[f"{TA_index}_{genome}"]["tox_early_stop"] = tmp_early_stop
                    d_res[f"{TA_index}_{genome}"]["tox_in_CDS"] = tmp_tox_isCDS


            if df_tblastn_results[(df_tblastn_results["qseqid"] == df_TA.loc[TA_index,"Antitoxin_name"]) & (df_tblastn_results["genome_name"] == genome)].empty == False:
                df_tmp = df_tblastn_results[(df_tblastn_results["qseqid"] == df_TA.loc[TA_index,"Antitoxin_name"]) & (df_tblastn_results["genome_name"] == genome)].reset_index(drop = True)
                # only one hit
                if len(df_tmp) == 1:
                    
                    d_res[f"{TA_index}_{genome}"]["antitox_cov"] = f"{len(df_tmp.loc[0,'sseq'])}/{len(ref_antitoxin_seq)}"
                    d_res[f"{TA_index}_{genome}"]["antitox_identity"] = df_tmp.loc[0,'pident']
                    d_res[f"{TA_index}_{genome}"]["antitox_early_stop"] = check_early_stop(df_tmp.loc[0,'sseq'])
                    str_antitox_ali += f"{genome}\t\t{' '*int(df_tmp.loc[0,'qstart'])}{df_tmp.loc[0,'sseq']}\n"
                    d_res[f"{TA_index}_{genome}"]["antitox_in_CDS"] = check_is_part_of_CDS(df_tmp.iloc[0], len(ref_antitoxin_seq),f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta")

                elif len(df_tmp) > 1:
                    n_tmp_index = 0
                    for row in df_tmp.iterrows():
                        n_tmp_index += 1
                        if n_tmp_index == 1 :
                            tmp_cov = f"{len(row[1]['sseq'])}/{len(ref_antitoxin_seq)}"
                            tmp_id = row[1]['pident']
                            tmp_early_stop = check_early_stop(row[1]['sseq'])
                            str_antitox_ali += f"{genome}\t\t{' '*int(row[1]['qstart'])}{row[1]['sseq']}\n"
                            tmp_frame = row[1]['sframe']
                            tmp_qstart = row[1]['qstart']
                            tmp_qend = row[1]['qend']
                            tmp_antitox_isCDS = [check_is_part_of_CDS(row[1], len(ref_antitoxin_seq),f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta")]

                        elif n_tmp_index > 1 :
                            #keep only if this is another part of the toxin and it is not far from the best hit (and same strand but not frame)
                            if (tmp_frame in [-1, -2, -3] and row[1]["sframe"] in [-1, -2, -3]) or (tmp_frame in [1, 2, 3] and row[1]["sframe"] in [1, 2, 3]):
                                same_strand = True
                            else :
                                same_strand = False

                            if (row[1]['qstart'] < tmp_qstart or row[1]['qend'] > tmp_qend) and same_strand == True:
                                tmp_cov = f"{tmp_cov.split('/',1)[0]}-{len(row[1]['sseq'])}/{tmp_cov.split('/',1)[1]}"
                                tmp_id = f"{tmp_id}-{row[1]['pident']}"
                                if tmp_early_stop == "No":
                                    tmp_early_stop = check_early_stop(row[1]['sseq'])
                                str_antitox_ali += f"{' '*len(genome)}\t\t{' '*int(row[1]['qstart'])}{row[1]['sseq']}\n"
                                
                                tmp_antitox_isCDS.append(check_is_part_of_CDS(row[1], len(ref_antitoxin_seq),f"{outdir}/TATdecay/tmp/tblastn_db/{TA_index}_all_genomes_spots.fasta"))

                                if tmp_frame != row[1]["sframe"]:
                                    d_res[f"{TA_index}_{genome}"]["antitox_frameshift"] = "Yes"
                    
                    d_res[f"{TA_index}_{genome}"]["antitox_cov"] = tmp_cov
                    d_res[f"{TA_index}_{genome}"]["antitox_identity"] = tmp_id
                    d_res[f"{TA_index}_{genome}"]["antitox_early_stop"] = tmp_early_stop
                    d_res[f"{TA_index}_{genome}"]["antitox_in_CDS"] = (",").join(tmp_antitox_isCDS)


        with open(f"{outdir}/TATdecay/results/1-spot_tblastn/TA_alignments/{TA_index}_toxin_alignments.txt", "w") as f:
            f.write(str_tox_ali)
        with open(f"{outdir}/TATdecay/results/1-spot_tblastn/TA_alignments/{TA_index}_antitoxin_alignments.txt", "w") as f:
            f.write(str_antitox_ali)
    
    return d_res



def check_early_stop(sequence2test):
    if "*" in sequence2test[:len(sequence2test)-1]:
        return "Yes"
    else :
        return "No"



def check_is_part_of_CDS(row, len_query, spot_fasta_file):

    #function which return "CDS" or "Nucleic" + (nucl_start:nucl_end) also, if "CDS" add _({prt_start_in_CDS}/{prt_qstart}:{prt_end_in_CDS}/{prt_qend})
    sseqid = row["sseqid"]
    start = row["sstart"]
    end = row["send"]

    spot_seqio = SeqIO.index(spot_fasta_file, "fasta")
    interval_seq = spot_seqio[sseqid].seq
    list_codon_start = ["ATG", "GTG", "TTG"]
    list_codon_stop = ["TGA", "TAG", "TAA"]

    #depending on the query length, we will be looking for start/stop codon farthest from the tblastn hit (up to 33%)
    len_query_33percent = len_query
    if start < end: #tblastn hit located on the positive strand
        #starting to assess whether there are a start codon which could be either the first codon of the tblastn hit or before the tblastn hit
        tblastn_seq = interval_seq[start-1:end]
        new_start = 0
        CDS_hit_start = start-1
        CDS_hit_end = end

        while tblastn_seq[:3] not in list_codon_start and new_start+3 <= len_query_33percent and start-1-(new_start+3) >= 3:
            new_start += 3
            tblastn_seq = interval_seq[start-1-new_start:end]
            CDS_hit_start = start-1-new_start

            if tblastn_seq[:3] in list_codon_stop:
                break

        if tblastn_seq[:3] not in list_codon_start:
            # in this case we will check whether we can find a codon start after our hit start
            tblastn_seq = interval_seq[start-1:end]
            new_start = 0
            while tblastn_seq[:3] not in list_codon_start and new_start+3 < end :
                new_start += 3
                tblastn_seq = interval_seq[start-1+new_start:end]

                if tblastn_seq[:3] in list_codon_stop:
                    break
            
            CDS_hit_start = start-1+new_start

        if tblastn_seq[:3] not in list_codon_start:
            return f"nucleic({start}:{end})_qcover({row['qstart']}:{row['qend']})"
        
        # now we're looking for a codon stop if we detected a codon start
        #NOTE: we're looking for codon stop from the first codon in the tblastn hit
        if end+len_query_33percent <= len(interval_seq)-1:
            tblastn_seq = interval_seq[CDS_hit_start:end + len_query_33percent]
        else :
            #we check here how many codons there are before the end of the contig
            tblastn_seq = interval_seq[CDS_hit_start:end + (len(interval_seq)-1-end)//3*3]

        codon2test = [0,3]
        while tblastn_seq[codon2test[0]:codon2test[1]] not in list_codon_stop and codon2test[1] +3 < len(tblastn_seq):
            codon2test[0] += 3
            codon2test[1] += 3
            CDS_hit_end = CDS_hit_start + codon2test[1]

        if interval_seq[CDS_hit_end-3:CDS_hit_end] not in list_codon_stop or (CDS_hit_end-CDS_hit_start)//3 < 20: # we do not consider CDS under 20 codons
            return f"nucleic({start}:{end})_qcover({row['qstart']}:{row['qend']})"

        else :
            prt_in_CDS = check_query_partofCDS(CDS_hit_start,CDS_hit_end, start, end, row["qstart"], row["qend"])
            return f"CDS({CDS_hit_start}:{CDS_hit_end})_{prt_in_CDS}"


    if start > end: #tblastn hit located on the negative strand

        tblastn_seq = interval_seq[end-1:start]
        new_start = start
        CDS_hit_start = start
        CDS_hit_end = end-1

        while reverse_codon(tblastn_seq[-3:]) not in list_codon_start and new_start+3 <= len(interval_seq)-1 and new_start+3 <= start + len_query_33percent :
            new_start += 3
            tblastn_seq = interval_seq[end-1:new_start]
            CDS_hit_start = new_start

            if reverse_codon(tblastn_seq[-3:]) in list_codon_stop:
                tblastn_seq = interval_seq[end-1:start]
                break
        
        # if we didn't find a start codon, we're looking into the hit
        new_start = start
        while reverse_codon(tblastn_seq[-3:]) not in list_codon_start and new_start-3 >= end:
            new_start -= 3
            tblastn_seq = interval_seq[end-1:new_start]
            CDS_hit_start = new_start

            if reverse_codon(tblastn_seq[-3:]) in list_codon_stop:
                break
        
        if reverse_codon(tblastn_seq[-3:]) not in list_codon_start:
            return f"nucleic({start}:{end})_qcover({row['qstart']}:{row['qend']})"
        
        # if we have a start codon, we're looking for a stop codon
        if end - len_query_33percent-1 >= 0:
            tblastn_seq = interval_seq[end-len_query_33percent-1:CDS_hit_start]

        else :
            tblastn_seq = interval_seq[end-end//3*3-1:CDS_hit_start]
        
        codon2test = [len(tblastn_seq)-3,len(tblastn_seq)]
        counted_codon = 1
        while reverse_codon(tblastn_seq[codon2test[0]:codon2test[1]]) not in list_codon_stop and codon2test[0] > 0:
            codon2test[0] -=3
            codon2test[1] -=3
            counted_codon += 1
            CDS_hit_end = CDS_hit_start - counted_codon*3

        tblastn_seq = interval_seq[CDS_hit_end: CDS_hit_start]
        if reverse_codon(tblastn_seq[:3]) not in list_codon_stop or (CDS_hit_start-CDS_hit_end)//3 < 20:
            return f"nucleic({start}:{end})_qcover({row['qstart']}:{row['qend']})"
        
        else :
            prt_in_CDS = check_query_partofCDS(CDS_hit_start,CDS_hit_end, start, end, row["qstart"], row["qend"])
            return f"CDS({CDS_hit_start}:{CDS_hit_end})_{prt_in_CDS}"


def reverse_codon(codon):
    
    reverse = ""
    for base in reversed(codon):
        if base == "A":
            reverse += "T"
        elif base == "T":
            reverse += "A"
        elif base == "C":
            reverse += "G"
        elif base == "G":
            reverse += "C"
    
    return reverse


def check_query_partofCDS(nucl_CDS_start, nucl_CDS_end, tblast_hit_start, tblast_hit_end, prt_qstart, prt_qend):

    prt_start_inCDS = prt_qstart
    prt_end_inCDS = prt_qend

    if nucl_CDS_end > nucl_CDS_start: #positive strand
        
        if nucl_CDS_start > tblast_hit_start:
            prt_start_inCDS += (nucl_CDS_start-tblast_hit_start)//3
        if nucl_CDS_end < tblast_hit_end:
            prt_end_inCDS -= (tblast_hit_end-nucl_CDS_end)//3
    
    elif nucl_CDS_end < nucl_CDS_start: #negative strand
        
        if nucl_CDS_start < tblast_hit_start:
            prt_start_inCDS += (tblast_hit_start-nucl_CDS_start)//3
        if nucl_CDS_end > tblast_hit_end:
            prt_end_inCDS -= (nucl_CDS_end-tblast_hit_end)//3
    
    return f"prt({prt_start_inCDS}/{prt_qstart}:{prt_end_inCDS}/{prt_qend})"