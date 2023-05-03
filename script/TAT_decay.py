#!/usr/bin/env python
import argparse
import os

import utilities
import TA_analysis
import decay_analysis
import align_color
import align_color_v2


def main():

    args = get_args()

    utilities.create_dir(args.outdir)
    n_cpu = utilities.core2use(args.cpu)
    df_TA = utilities.TA_parser(args.TAseq, args.outdir)
    df_TA = TA_analysis.TA_ref_locations(df_TA, args.TAseq, args.replicons, args.core, args.outdir, ref = args.ref if args.ref != "First" else "First")
    utilities.generate_spot_db(df_TA, args.core, args.replicons, args.outdir, cpu = n_cpu)
    TA_analysis.all_genomes_search(df_TA, args.outdir, cpu = n_cpu)
    df_tblastn = decay_analysis.is_both_components(df_TA, args.outdir, cpu = n_cpu)

    decay_analysis.check_component_integrity(df_TA, df_tblastn, args.outdir, cpu = n_cpu)

    """
    for file in os.listdir(f"{args.outdir}/TATdecay/results/1-spot_tblastn/TA_alignments"):
        if file.endswith(".txt"):
            try :
                align_color.ali2color(f"{args.outdir}/TATdecay/results/1-spot_tblastn/TA_alignments/{file}")
            except:
                print(f"Alignemnt of {file} cannot be colored yet")
                continue
    """

    align_color_v2.ali2color(args.outdir+"/TATdecay/results/1-spot_tblastn/TA_alignments", args.outdir+"/TATdecay/results/1-spot_tblastn/all_TAs_decays.tsv", 
                             order = args.order,
                             names= args.names,
                             ref = args.ref
                             )




def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("--TAseq", "-t",
                        help = " (REQUIRED) Proteins sequences of the TAs to search (Toxin seq then antitoxin seq) (fna format)", required = True)
    parser.add_argument("--core", "-c",
                        help = " (REQUIRED) Core_features_file (tsv format)", required = True)
    parser.add_argument("--replicons", "-r",
                        help = " (REQUIRED) Replicons folder (fasta format)", required = True)
    parser.add_argument("--outdir", "-o",
                        help = " (REQUIRED) Out directory folder", required = True)
    parser.add_argument("--ref",
                        help = "Reference genome to use (if not provided, the first alphabetical one will be used, default: 'First')",
                        default = "First")
    parser.add_argument("--order",
                        help = "A file specifiying the order of each genome for the final figure (one per line, default = None)",
                        default = None)
    parser.add_argument("--names",
                        help = "TSV files with two columns, the first with the analysis names and the other with the real names (no columns header, default = None)",
                        default = None)
    parser.add_argument("--cpu",
                    help = "Number of cores to use (default = 1)",
                    default = 1)

    args = parser.parse_args()
    return args




if __name__ == "__main__":
    main()














