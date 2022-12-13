#!/usr/bin/env python
### This is modified from Charles Y. Lin CRC program: https://github.com/linlabcode/CRC
### This can be normaly installed through python and only used for calling subpeaks/valleys inside enhancers
###
#######

"""CRC main function call."""

import argparse
import sys

from crc.resources import crc_utils, utils


def parse_args(args=None):
    """Argument parser."""
    parser = argparse.ArgumentParser(
        usage=(
            "crc [options]"
            " -e [ENHANCER_FILE]"
            " -c [CHROMOSOMES_FOLDER_PATH]"
            " -g [GENOME]"
            " -o [OUTPUTFOLDER]"
            " -n [NAME]"
        )
    )

    # Required flags
    parser.add_argument("-e", "--enhancer_file", dest="enhancers", default=None, type=str,
                        help=(
                            "Provide a ROSE generated enhancer table "
                            "(_AllEnhancers_ENHANCER_TO_GENE.txt)"
                        ), required=True)
    parser.add_argument("-g", "--genome", dest="genome", default=None, type=str,
                        help=(
                            "Provide the build of the genome to be used for the analysis. "
                            "Currently supports HG19, MM10, and RN6"
                        ), required=True)
    parser.add_argument("-c", "--chrom-path", dest="chrom_path", type=str,
                        help=(
                            "Provide a path to a folder with a seperate fasta file for each "
                            "chromosome"
                        ), required=True)
    parser.add_argument("-o", "--output", dest="output", default=None, type=str,
                        help="Enter an output folder", required=True)
    parser.add_argument("-n", "--name", dest="name", default=None, type=str,
                        help="Provide a name for the job", required=True)

    # Either bams for valleys or subpeaks are needed: This is Required
    parser.add_argument("-b", "--bam", dest="bam", default=None, type=str,
                        help="Enter a comma separated list of bams of valley finding",
                        required=True)

    # Additional options: Unnecessary
    parser.add_argument("-m", "--mask", dest="mask_file", type=str,
                        help="Masking file in BED format", default=None, required=False)
    parser.add_argument("-a", "--activity", dest="activity", default=None, type=str,
                        help="A table with active gene names in the first column", required=False)
    parser.add_argument("-l", "--extension-length", dest="extension", default=100, type=int,
                        help=(
                            "Enter the length to extend subpeak regions for motif finding. "
                            "default is 100"
                        ), required=False)
    parser.add_argument("-N", "--number", dest="number", default=1, type=int,
                        help=(
                            "Enter the number of non overlapping motifs in a region required "
                            "to assign a binding event. Default=1"
                        ), required=False)
    parser.add_argument("--motifs", dest="motifs", default=False, type=str,
                        help="Enter additional PWM file for the analysis", required=False)
    parser.add_argument("--config", dest="config", default='', type=str,
                        help="Enter genome configuration file to overwrite default paths",
                        required=False)

    return parser.parse_args(args)


def crc(enhancers, genome_input, chrom_path, output, analysis_name, bam,
        mask_file=None, activity_path=None, const_extension=100, number=1, motifs=False,
        config=''):
    """CRC main function."""
    # =====================================================================================
    # ===============================I. PARSING ARGUMENTS==================================
    # =====================================================================================
## CRC defined genome class
    genome = crc_utils.load_genome(
        genome_input,
        chrom_path,
        mask_file=mask_file,
        config_file=config,
    )

    motif_database_file = genome.return_feature('motif_database')
    motif_convert_file = genome.return_feature('motif_convert')

    # User input files
    enhancer_file = enhancers

    if bam is None and subpeak_file is None:
        print('ERROR: Must provide either bams for valley finding or subpeaks as a .bed')
        sys.exit()

    # bam file is necessary for valleys finding
    if bam:
        bam_file_list = [bam_path for bam_path in bam.split(',') if bam_path]
        print(bam_file_list)
    else:
        print('ERROR: Must provide bams for valley finding or subpeaks')
        sys.exit()

    # Output folder and analysis name
    print(output)
    output_folder = utils.format_folder(output, True)

    print(
        '\n\n#======================================\n#===========I AM BRINGING YOU DATA SUMMARY============\n#='
        '=====================================\n'
    )

    print('Finding valleys/subpeaks for {}'.format(analysis_name))
    print('Writing output to {}'.format(output_folder))
    print('Using {} to define active genes'.format(activity_path))

    # =====================================================================================
    # =======================II. IDENTIFYING CANDIDATE TFS AND NODES=======================
    # =====================================================================================

    print(
        '\n\n#======================================\n#===I AM MAPPING GENES AND ENHANCERS====\n#='
        '=====================================\n'
    )

    (
        gene_region_table,
        gene_tf_region_table,
        enhancer_region_table,
        enhancer_tf_region_table,
        gene_summary_table,
        candidate_tf_list,
        gene_to_enhancer_dict,
    ) = crc_utils.gene_to_enhancer(genome, enhancer_file, activity_path)


    # =====================================================================================
    # ==========================III. FINDING VALLEYS/SUBPEAKS==============================
    # =====================================================================================

    print(
        '\n\n#======================================\n#=====III. NOW I AM FINDING VALLEYS/SUBPEAKS====\n#='
        '=====================================\n'
    )

    # So here we would need to find valleys everywhere
    print('finding valleys ...... ')
    # Note: the tf_bed_path is for networks, all is for out degree finding
    all_bed_path = crc_utils.find_valleys(
        gene_to_enhancer_dict, bam_file_list, analysis_name, output_folder, cutoff=0.2
    )

    # Generate subpeaks fasta files
    all_sub_bed, all_fasta = crc_utils.generate_subpeak_fasta(
        gene_to_enhancer_dict, all_bed_path, genome, analysis_name, const_extension
    )
    # This is the case where we did valleys # only reason you would need to output the sub bed
    all_sub_out = '{}{}_all_subpeak.bed'.format(output_folder, analysis_name)
    utils.unparse_table(all_sub_bed, all_sub_out, '\t')

    # Writing the all subpeak fasta out to disk
    all_fasta_out = '{}{}_all_subpeak.fasta'.format(output_folder, analysis_name)
    utils.unparse_table(all_fasta, all_fasta_out, '')

    print(
        '\n\n#======================================\n#======FINDING VALLEYS IS DONE=======\n#='
        '=====================================\n'
    )



def main(args=None):
    """Main function call."""
    args = parse_args(args)

    print(args)

    crc(
        args.enhancers,
        args.genome,
        args.chrom_path,
        args.output,
        args.name,
        bam=args.bam,
        subpeak_file=args.subpeaks,
        mask_file=args.mask_file,
        activity_path=args.activity,
        const_extension=args.extension,
        number=args.number,
        motifs=args.motifs,
        config=args.config,
    )
