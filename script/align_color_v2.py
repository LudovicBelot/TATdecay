import html
from operator import index
import sys
import pandas as pd
from collections import OrderedDict
import copy
import  glob


def main():

    file2color = sys.argv[1]
    df_tblastn = pd.read_csv(sys.argv[2])
    str_html_ali = ali2color(file2color)




def ali2color(align_folder, tblastn_res_file, **kwargs):

    order_file = kwargs.get("order", None) #genome file to order the alignment results if needed
    genome_name_file = kwargs.get("names", None)
    ref_genome = kwargs.get("ref", "First")
    list_antitoxins_alignments_files = [f for f in glob.glob(f"{align_folder}/*_antitoxin_alignments.txt")]
    list_toxins_alignments_files = [f for f in glob.glob(f"{align_folder}/*_toxin_alignments.txt")]
    df_tblastn_res = pd.read_csv(tblastn_res_file, sep = "\t", comment = "#")
    df_tblastn_res = df_tblastn_res.set_index("Unnamed: 0")
    list_TAs_index = df_tblastn_res.index.map(lambda x: x.split("_")[0]).drop_duplicates().tolist()


    if ref_genome == "First":
        #quickly opening an alignment file to get the name of the ref genome
        with open(list_toxins_alignments_files[0], "r") as f:
            index_line = 0
            for line in f:
                index_line +=1
                if index_line == 2:
                    ref_genome = line.split("\t")[0].strip()
                    break

    list_genomes = df_tblastn_res["genome_name"].drop_duplicates().tolist()
    if order_file != None :
        tmp_list_genomes = copy.deepcopy(list_genomes)
        list_genomes = [ref_genome]
        with open(order_file, "r") as f:
            for line in f:
                if line.strip() in tmp_list_genomes:
                    list_genomes.append(line.strip())

    od_genomes = OrderedDict()
    od_genomes["Reference"] = {"real_name":list_genomes[0]}
    if genome_name_file == None:
        for genome in list_genomes:
            od_genomes[genome] = {"real_name":genome}
    else :
        df_names_genomes = pd.read_csv(order_file, sep = "\t", comment = "#", names = ["analysis_name", "real_name"])
        for genome in list_genomes:
            od_genomes[genome] = {"real_name":df_names_genomes.loc[df_names_genomes[df_names_genomes["analysis_name"] == genome].drop_duplicates().index, "real_name"]}

    # Now working to create a html file with a coloured coded alignment of both toxin/antitoxin
    for TA_index in list_TAs_index:
        ali2coloured_html(align_folder+f"/{TA_index}_toxin_alignments.txt", TA_index, od_genomes, df_tblastn_res)
        ali2coloured_html(align_folder+f"/{TA_index}_antitoxin_alignments.txt", TA_index, od_genomes, df_tblastn_res)




def ali2coloured_html(align_file, TA_index, od_genomes_TA, df_tblastn_res):

    str_html = "<!DOCTYPE html>\n<html>\n<body>"
    last_genome = ""
    index_line = 0
    list_genomes_with_hits = ["Reference"]

    with open(align_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                str_html += f"<b>{line.strip()}</b>\n<table>\n"
                if line.startswith("#Antitoxin"):
                    TA_component = "antitox"
                else :
                    TA_component = "tox"
                index_line += 1
                continue
            
            elif line.startswith(" ") == False:
                if index_line == 1:
                    od_genomes_TA["Reference"]["seq"] = [line.split("\t\t ")[-1]]
                    od_genomes_TA["Reference"]["isCDS"] = ["CDS"]
                    index_line += 1
                
                else :
                    count_tblastn_hit = 0
                    last_genome = line.split("\t\t ")[0]
                    list_genomes_with_hits.append(last_genome)

                    od_genomes_TA[last_genome]["seq"] = [line.split("\t\t ")[-1]]
                    od_genomes_TA[last_genome]["isCDS"] = [df_tblastn_res.loc[TA_index+'_'+last_genome, TA_component+"_in_CDS"].split(",")[count_tblastn_hit]]

            elif line.startswith(" "):
                count_tblastn_hit += 1
                od_genomes_TA[last_genome]["seq"].append(line.split("\t\t ")[-1])
                od_genomes_TA[last_genome]["isCDS"].append(df_tblastn_res.loc[TA_index+'_'+last_genome, TA_component+"_in_CDS"].split(",")[count_tblastn_hit])

        # adding empty values to the dict for genomes without tblastn hits
        for k in od_genomes_TA.keys():
            if k not in list_genomes_with_hits:
                od_genomes_TA[k]["seq"] = []
                od_genomes_TA[k]["isCDS"] = []


    #Now we're creating the html file based on the data we gathered
    for genome_name, genome_values in od_genomes_TA.items():
        str_html += "<tr>\n"
        if genome_name == "Reference":
            str_html += "<td style='line-height:0.3em;'>"+ "Reference" + "</td>"+"<td><span style='background-color: white'>" + genome_values["seq"][0] +"</span></td>\n"
            str_html += "</tr>\n"
            ref_seq = genome_values["seq"][0]

        else :
            #First checking whether we have a tblastn hit within a CDS in the given genome.
            # AA conserved in a CDS will be represented by a black box, if it is conserved but not within a CDS => grey, and if not conserved in white
            # NOTE that if there are multiples tblastn hits we keep only the first/best one encoded within a CDS and represents it in black, others secondary CDS will be in grey aswell
            # NOTE 2 : in case of overlapping tblastn hit we will priorize the representation of the best CDS
            best_hit_index = None
            tmp_best_hit_index = 0

            for tested_tblastnhit4CDS in genome_values["isCDS"]:
                if tested_tblastnhit4CDS.startswith("CDS"):
                    best_hit_index = tmp_best_hit_index
                    break
                else :
                    tmp_best_hit_index += 1
            
            if best_hit_index != None:
                best_hit = CDS(genome_values["isCDS"][best_hit_index],genome_values["seq"][best_hit_index])
                #now we're adding each all others tblastn hit which do not overlap with the best CDS
                #NOTE : put this part on hold because I'm not sure I want to represent secondary results onto the figure
                #for other hit in genome_values["isCDS"]:
            
            elif best_hit_index == None and genome_values["isCDS"] != [] :
                best_hit = tblastn_hit(genome_values["isCDS"][0], genome_values["seq"][0])
            
            else :
                str_html += "<td style='line-height:0.3em;'>"+ genome_values["real_name"] +f"<td><span style='color: white'>{ref_seq}</span></td>"+"</td>\n"
                str_html += "</tr>\n"
                continue

            #html representation
            str_html += "<td style='line-height:0.3em;'>"+ genome_values["real_name"] + "</td><td>"
            n_index_aa = 0
            if best_hit.type == "CDS":
                for aa in best_hit.seq:
                    if aa != "\n":
                        if n_index_aa >= best_hit.CDS_start-1 and n_index_aa <= best_hit.CDS_end-1:
                            aa_within_CDS = True
                        else :
                            aa_within_CDS = False

                        str_html += color_aa(ref_seq[n_index_aa], best_hit.seq[n_index_aa], aa_within_CDS)
                        n_index_aa += 1
            
            elif best_hit.type != "CDS":
                for aa in best_hit.seq:
                    if aa != "\n":
                        str_html += color_aa(ref_seq[n_index_aa], best_hit.seq[n_index_aa], False)
                        n_index_aa += 1
             
            str_html += "</span></td></tr>\n"
    
    str_html += "</table>\n</html>\n</body>"

    with open(align_file.split(".txt")[0]+".html", "w") as f:
        f.write(str_html)





def color_aa(ref_aa, new_aa, within_CDS):

    if new_aa == " ":
        return "<span style='color: white'>" + ref_aa +"</span>"
    else :
        if new_aa == ref_aa and within_CDS == True:
            return "<span style='background-color: black;color: black'>" + ref_aa +"</span>"
        
        elif new_aa == ref_aa and within_CDS == False:
            return "<span style='background-color: grey;color: grey'>" + ref_aa +"</span>"
        
        elif new_aa != ref_aa:
            return "<span style='background-color: white;color: white'>" + ref_aa +"</span>"



class CDS:
    def __init__(self,CDS_prt_res,seq):
        tmp_CDS_res = CDS_prt_res.split("prt(")[-1].split(")")[0]
        self.type = "CDS"
        self.CDS_start = int(tmp_CDS_res.split("/")[0])
        self.CDS_end = int(tmp_CDS_res.split(":")[-1].split("/")[0])
        self.TB_start = int(tmp_CDS_res.split(":")[0].split("/")[-1])
        self.TB_end = int(tmp_CDS_res.split(":")[-1].split("/")[-1])
        self.seq = seq

class tblastn_hit:
    def __init__(self,tb_hit_res,seq):
        tmp_tb_res = tb_hit_res.split("qcover(")[-1].split(")")[0]
        self.type = "tblastn"
        self.TB_start = int(tmp_tb_res.split(":")[0])
        self.TB_end = int(tmp_tb_res.split(":")[-1])
        self.seq = seq