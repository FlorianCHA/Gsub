#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Florian Charriat


# Import subprocess for launch OS command from script, here it's for launch tbl2asn
import subprocess
# Import orffinder package for predict orf
from orffinder import orffinder
# For parse source file
import pandas as pd
# For path treatment
from pathlib import Path
# Import some biopython package to process fasta sequences
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# Import gooey for graphical argument parser
from gooey import Gooey
from gooey import GooeyParser
import sys
import os
import shutil
# For know the OS system
import platform
# Color package
import colored
from colored import stylize, attr, fg
# Time for know how second tools take for generate ASN
import time
# for removing warning python
import warnings
# For polymerase detection
if platform.system() == 'Linux':
    import pyhmmer

# For add color in terminal
colored.set_tty_aware(False)

# Remove warning python
warnings.filterwarnings("ignore")


def formated_error(message):
    """
    This function use stylize for format error message
    """
    error_message = stylize(message, fg('red') + attr('bold') + attr('underlined'))
    return error_message


def choose_tbl2asn():
    """
    This tools verify the OS system and choose the correct tbl2asn to launch for your computer
    """
    if platform.system() == 'Linux':
        return f'{Path(__file__).resolve().parent}/tools/table2asn.linux'
    elif platform.system() == 'Windows':
        return f'{Path(__file__).resolve().parent}/tools/table2asn.windows//table2asn.exe'
    else:
        raise TypeError(formated_error(f'You OS system are not support,'
                                       f' for now this package works only in Linux and Windows.'))


def fasta2dict(filename):
    """
    This function take a file path (fasta), and return a dictionary of sequence
    """
    with open(filename, "r") as fastaFile:
        return SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))


def create_fasta(biopython_object, output_file, description):
    """
     This function take a biopython sequence in input and create a fasta
    """
    with open(output_file, 'w') as f:
        record = SeqRecord(biopython_object.seq, id=biopython_object.id,
                           name=biopython_object.id, description=description)
        SeqIO.write(record, f, "fasta")


def search_orf_orffinder(biopython_object, length_min=300):
    """
    This function take a biopython sequence in input and create a dict object which contains all information about
    orf predict by orrfinder package (near that ncbi result).
    The minimum length of orf is by default to 75, you can change this with the length_min=XX option.
    """
    # Create the orffinder object with biopython object
    orf_result = orffinder.getORFs(biopython_object, minimum_length=length_min)
    # If the tools not detect orf, the function return a None object
    if not orf_result:
        return None
    # Create liste of gene include in other gene which must be delete
    liste_to_remove = []
    # Compare orf with previously good orf for see if orf is inclued in the previously good orf
    # For frame "+"
    if '+' in [elt["sense"] for elt in orf_result]:
        orf_result_sens_pos = [elt for elt in orf_result if elt['sense'] == '+']
        # Sort result in function of orf start
        orf_result_sens_pos = sorted(orf_result_sens_pos, key=lambda x: x['start'])
        start_orf_prev = orf_result_sens_pos[0]['start']
        end_orf_prev = orf_result_sens_pos[0]['end']
        for orf in orf_result_sens_pos[1:]:
            if orf['start'] >= start_orf_prev and orf['end'] <= end_orf_prev:
                liste_to_remove.append(orf)  # If orf is inclued in previously orf, he is remove or result
            else:  # If orf is not inclued, we change the previously start and end orf by this one
                start_orf_prev = orf['start']
                end_orf_prev = orf['end']

    # Remove all orf inclued in an another orf
    # For frame "-"
    if '-' in [elt["sense"] for elt in orf_result]:
        orf_result_sens_neg = [elt for elt in orf_result if elt['sense'] == '-']
        # Sort result in function of orf start
        orf_result_sens_neg = sorted(orf_result_sens_neg, key=lambda x: x['start'])
        start_orf_prev = orf_result_sens_neg[0]['start']
        end_orf_prev = orf_result_sens_neg[0]['end']
        for orf in orf_result_sens_neg[1:]:
            if orf['start'] >= start_orf_prev and orf['end'] <= end_orf_prev:
                liste_to_remove.append(orf)  # If orf is inclued in previously orf, he is remove or result
            else:  # If orf is not inclued, we change the previously start and end orf by this one
                start_orf_prev = orf['start']
                end_orf_prev = orf['end']
    # Remove all orf inclued in an another orf
    liste_final = [orf for orf in orf_result if orf not in liste_to_remove]


    # Create dico result for other function
    dico = {}
    for i, orf in enumerate(sorted(liste_final, key=lambda x: x['start'])):
        id_orf = f"ORF{i + 1}"
        dico[id_orf] = {"start": orf['start'],  # Start of ORF
                        "end": orf['end'],  # End of ORF
                        "end_partial": orf['trailing'],  # True if end of ORF is truncated
                        "partial": orf['trailing'],  # partial is True if start or end is truncated
                        "strand": orf['sense']}
    return dico


def remove_overlaps_gene(dico):
    """
    This function take the dict output from search_orf_orffinder and remove overlaps gene for kept only the longest
    """
    # Transform dico with length as keys :
    dico_length = {}
    liste_remove = []
    for id_orf in dico:
        length = max(dico[id_orf]['start'], dico[id_orf]['end']) - min(dico[id_orf]['start'], dico[id_orf]['end'])
        dico_length[length] = {"start": dico[id_orf]['start'], 'end': dico[id_orf]['end'], 'id': id_orf}

    # Create list for sort longest to smallest
    liste_length = [length for length in dico_length]
    for length in sorted(liste_length, reverse=True):
        # Defines the smallest position as the start and the largest position as the end
        start = min(dico_length[length]['start'], dico_length[length]['end'])
        end = max(dico_length[length]['start'], dico_length[length]['end'])
        for length_2 in dico_length:
            start_2 = min(dico_length[length_2]['start'], dico_length[length_2]['end'])
            end_2 = max(dico_length[length_2]['start'], dico_length[length_2]['end'])
            # Look if the two genes do not overlap and if one of the two genes is not already in the list of removes,
            # in this case the gene of the list is not to take into account
            if (start < start_2 < end or start < end_2 < end) and length not in liste_remove and length > length_2:
                liste_remove.append(length_2)
    # Retrieve ID orf from length id
    liste_id_to_remove = [dico_length[length]['id'] for length in liste_remove]

    i = 0
    dico_final = {}
    for id_orf in dico:
        if id_orf not in liste_id_to_remove:
            i += 1
            new_id = f'{id_orf.split("F")[0]}F{i}'  # Rename gene
            dico_final[new_id] = dico[id_orf]
    return dico_final


def remove_strand_gene(dico):
    """
    This function take the dict output from search_orf_orffinder and remove gene with different frame
    """
    # Retrieve all strand for all orf in the sequence
    liste_strand = [dico[id_orf]["strand"] for id_orf in dico]

    # Select the majority strand
    strand_majority = [strand for strand in liste_strand if liste_strand.count(strand) > 0.5 * len(liste_strand)]
    if not strand_majority:
        dico_length = {'+': [], '-': []}
        for id_orf in dico:
            length = max(dico[id_orf]['start'], dico[id_orf]['end']) - min(dico[id_orf]['start'], dico[id_orf]['end'])
            strand = dico[id_orf]['strand']
            dico_length[strand].append(length)

        if sum(dico_length['-']) > sum(dico_length['+']):
            strand_majority = '-'
        else:
            strand_majority = '+'
    else:
        strand_majority = strand_majority[0]

    i = 0
    dico_final = {}
    for id_orf in dico:
        # Select all orf with the majority strand
        if dico[id_orf]['strand'] == strand_majority:
            i += 1
            new_id = f'{id_orf.split("F")[0]}F{i}'  # Rename gene
            dico_final[new_id] = dico[id_orf]

    return dico_final


def orf_pfam(biopython_object, output_aa, dico_orf, database, score_cutoff, eval_cutoff, cov_cutoff):
    """
    This function take biopython object and orf dico from search_orf or search_orf_finder output
    and use hmmseach for search polymerase
    """
    # Create fasta with all orf find
    with open(output_aa, 'w') as fasta_file:
        for orf_id in dico_orf:
            orf = dico_orf[orf_id]
            if orf['strand'] == '+':
                sequence = biopython_object.seq[orf['start'] - 1:orf['end'] - 1]
            if orf['strand'] == '-':
                sequence = biopython_object.seq[orf['end'] - 1:orf['start'] - 1].reverse_complement()
            record = SeqRecord(sequence.translate(), id=orf_id, name=orf_id, description='')
            SeqIO.write(record, fasta_file, "fasta")

    # For know length of each orf
    seq = fasta2dict(output_aa)
    len_seq = {id_seq: len(seq[id_seq].seq) for id_seq in seq}

    # Charge orf with pyHMMER package
    with pyhmmer.easel.SequenceFile(output_aa, digital=True) as seq_file:
        sequences = list(seq_file)

    # Search polymerase with polymerase pfam database
    dico_hits = {}
    for hmm in pyhmmer.plan7.HMMFile(database):
        pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
        hits = pipeline.search_hmm(hmm, sequences)
        if hits.hits_reported != 0:
            for hit in hits:
                start_best_domaine = hit.domains[0].alignment.target_from
                end_best_domaine = hit.domains[0].alignment.target_to
                len_best_domain = end_best_domaine - start_best_domaine
                if hit.score >= score_cutoff and hit.evalue <= eval_cutoff:
                    if len_best_domain / len_seq[hit.name.decode('utf-8')] >= cov_cutoff:
                        if hit.name.decode('utf-8') in dico_hits:
                            if dico_hits[hit.name.decode('utf-8')]['score'] < hit.score:
                                dico_hits[hit.name.decode('utf-8')] = {"description": hmm.description.decode('utf-8'),
                                                                       "score": hit.score}
                        else:
                            dico_hits[hit.name.decode('utf-8')] = {"description": hmm.description.decode('utf-8'),
                                                                   "score": hit.score}
                    else:
                        polymerase = f"{hmm.description.decode('utf-8')}"
                        if hit.name.decode('utf-8') in dico_hits:
                            if dico_hits[hit.name.decode('utf-8')]['score'] < hit.score:
                                dico_hits[hit.name.decode('utf-8')] = {"description": polymerase,
                                                                       "score": hit.score}
                        else:
                            dico_hits[hit.name.decode('utf-8')] = {"description": polymerase,
                                                                   "score": hit.score}
    return dico_hits


def orf_to_feature(dico_orf, contig, output_file, dico_pol):
    """
    This function take the output of search_orf function and create a feature file for GenBank submission at this format
    >Feature Contig1670
    185	1330	gene
                gene	Contig1670_1
    185	1330	CDS
                product	Contig1670_1
                not	partial hypothetical protein
                codon_start	2
    """
    with open(output_file, 'w') as output:
        text = f">Feature {contig}\n"  # Only one feature by contigs, it's header for this file
        for gene_id in dico_orf:
            info = dico_orf[gene_id]  # Retrieve in information by orf from dictionary input
            end_symbol = '>' if info['end_partial'] else ''  # Create a > if start is truncated in this ORF
            # The start and end symbol was used before the start or end information in feature file for show the side
            # which are truncted

            # They have difference between orffinder output and tbl input, so we remove 1 position at end for only
            # complete gene
            if info["partial"]:
                remove_end = 0
            else:
                remove_end = 1
            # Create note variable for know if it's polymerase, hypothetical protein or partial protein
            note = f"{dico_pol[gene_id]['description'] if gene_id in dico_pol else 'partial hypothetical protein' if info['partial'] else 'hypothetical protein'}"
            # This 2 condition are use for remove position at end or start in function of strand
            if info['strand'] == '+':
                text = text + \
                       f"{info['start']}\t{end_symbol}{info['end'] - remove_end}\tCDS\n" \
                       f"\t\t\tproduct\t{note}\n" \
                       f"\t\t\tnote\t{gene_id}\n"
            else:
                text = text + \
                       f"{info['start'] - remove_end}\t{end_symbol}{info['end']}\tCDS\n" \
                       f"\t\t\tproduct\t{note}\n" \
                       f"\t\t\tnote\t{gene_id}\n"
        output.write(text)


def parse_src_file(src_file):
    """
    This function take the path of src file and generate a data object with all compatible feature for each seq.
    In addition create another dico with only description and comment to add at each seq.
    And create a list of seq to submit
    """
    liste_compatible_feature = ['Sequence_ID', 'Organism', 'Strain', 'Country',
                                'Host', 'Collection_date', 'Definition', 'Molecule', 'Lineage', 'reverse']
    liste_feature = ['Sequence_ID', 'Organism', 'Strain', 'Country', 'Host', 'Collection_date']

    data_src = pd.read_csv(src_file, sep='\t', header=0)
    data_src = data_src.fillna('')

    # Verify if all information are present in source file
    missing_columns = [feature for feature in liste_compatible_feature if feature not in data_src.columns]
    if len(missing_columns) != 0:
        txt_missing_columns = ', '.join([f'"{elt}"' for elt in missing_columns])
        raise TypeError(formated_error(f"The '{src_file.split('/')[-1]}' isn't correct, "
                                       f"please check that file  they missing {txt_missing_columns} columns."))

    # Keept only sequence with polymerase ( score "1" on source file)
    # data_src = data_src[data_src['Polymerase'] == 1]
    # Create dico for description & comment columns
    dico_description = {data_src.loc[index, 'Sequence_ID']:
                        {'definition': data_src.loc[index, 'Definition'],
                         'strand': f'{"+" if float(data_src.loc[index, "reverse"]) > 0 else "-"}',
                         'lineage': data_src.loc[index, 'Lineage'],
                         'molecule': str(data_src.loc[index, 'Molecule'])}  # Remove nan values
                        for index in data_src.index}

    # Create data frame for src file by sequence

    data_src = data_src[liste_feature]

    return dico_description, data_src


def create_submit_file(fasta, template, src_file, output, minlengthorf, minlengthcontig, assembler, technology,
                       overlaps, frame, database, score_cutoff, eval_cutoff, cov_cutoff, genome):
    """
    This function take fasta input which contains all contigs to submit, the template ncbi,
    the directory output and the min length for a orf and launch all other function for create the ASN file.
    """
    # liste_kept is a list which contains all contigs with polymerase
    fasta_dict = fasta2dict(fasta)
    liste_kept = list(fasta_dict.keys())
    dico_description, data_src = parse_src_file(src_file)
    count = 0  # For initiate the count (progress bar)
    len_fasta = len(dico_description) + 1  # For know how many loop the script done (progress bar)
    for seq_id in fasta_dict:
        # Only seq with polymerase are in dico_description and a contigs > min length give by user
        if seq_id not in dico_description or len(fasta_dict[seq_id].seq) <= minlengthcontig:
            liste_kept.remove(seq_id)
            continue
        count += 1
        print_progress(count, len_fasta)  # Progress bar
        sys.stdout.flush()  # For print correctly progress bar in Graphical User interface (GUI)
        directory_tmp = f'{output}/{seq_id}'
        if not Path(directory_tmp).exists():
            os.mkdir(directory_tmp)
        # Create fasta for the contigs
        fasta_tmp = f'{directory_tmp}/{seq_id}.fasta'
        if dico_description[seq_id]['strand'] == "-":
            seq = fasta_dict[seq_id].seq
            seq = seq.reverse_complement()
            fasta_dict[seq_id] = SeqRecord(seq, id=fasta_dict[seq_id].id,
                           name=fasta_dict[seq_id].id, description="")
        create_fasta(fasta_dict[seq_id], fasta_tmp,
                     dico_description[seq_id]['definition'], )
        # Create source file for the contigs
        data_sample = data_src[data_src['Sequence_ID'] == seq_id]
        for column in data_sample.columns:
            if (data_sample[column] == '').item():
                data_sample = data_sample.drop(column, axis=1)
        data_sample.to_csv(f'{directory_tmp}/{seq_id}.src', sep='\t', index=False)
        # Create feature file of contigs
        tbl_tmp = f'{directory_tmp}/{seq_id}.tbl'
        dico_orf = search_orf_orffinder(fasta_dict[seq_id], minlengthorf)
        if not overlaps:
            dico_orf = remove_overlaps_gene(dico_orf)
        if not frame:
            dico_orf = remove_strand_gene(dico_orf)
        # Use PFAM for detect polymerase
        aa_tmp = f'{directory_tmp}/{seq_id}.faa'  # Temp file which contain orf protein
        # As Hmmer is not available on windows we test if os is linux or windows. If is windows we don't use orf_pfam.
        if platform.system() == 'Linux':
            dico_pol = orf_pfam(fasta_dict[seq_id], aa_tmp, dico_orf, database, score_cutoff, eval_cutoff, cov_cutoff)
        else :
            dico_pol = {}
        # Don't create feature file if no orf find
        if dico_orf is not None:
            orf_to_feature(dico_orf, seq_id, tbl_tmp, dico_pol)
        # Add individual comment from src file to the contig
        # Src line
        tbl2asn_tool = choose_tbl2asn()
        # Comment for Assembly information
        with open(f'{directory_tmp}/{seq_id}.cmt', 'w') as comment_file:
            comment_file.write('StructuredCommentPrefix\t##Assembly-Data-START##\n')
            comment_file.write(f'Assembly Method\t{assembler}\n')
            comment_file.write(f'Sequencing Technology\t{technology}\n')
            comment_file.write('StructuredCommentSuffix\t##Assembly-Data-END##\n')
        # Lineage = Clade; order (unclassifed virus if order == NA or _order)
        feature = ''
        if genome == 'Prokaryote':
            feature = feature + "[gcode=11]"
        molecule = dico_description[seq_id]['molecule'].upper()
        if molecule == 'RNA':
            feature = feature + " [molecule=rna]"
        circular = False
        if circular:
            feature = feature + " [topology=circular]"
        lineage = dico_description[seq_id]['lineage']
        feature = feature + f" [lineage={lineage}]"
        # Use tbl2asn tools for generate from fasta & tbl the seq file in ASN.1 format
        subprocess.run(f'{tbl2asn_tool} -t {template} -i {fasta_tmp} -j "{feature}" -V bv',
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
        #######################################################
    return liste_kept


def verif_quality(output, contigs_names, list_kept):
    """
    This function check if errorsummary.val is empty, else raise a warning for user
    """
    liste_warning = []
    liste_error = []

    with open(f'{output}/Error_validation_genebank.txt', 'w') as output_qual:
        for contig_name in contigs_names:
            if contig_name not in list_kept:
                liste_warning.append(stylize(f"\t* The contig {contig_name} isn't "
                                     f"process because they haven't polymerase\n", fg('dark_orange_3a') + attr('bold')))
                continue
            contig = f"{output}/{contig_name}"
            file = f'{contig}/{contig_name}.val'
            with open(file, 'r') as file_quality:
                liste_file = list(file_quality)
                if liste_file:
                    liste_error.append(stylize(f"\t* The contig {contig_name} isn't "
                                         f"process correctly, please verify the errorsummary.val in output\n",
                                                 fg('red') + attr('bold')))
                    txt_error = "\t".join(liste_file)
                    output_qual.write(f"* {contig_name} have some problem, here you can find the GenBank error :"
                                      f"\n{txt_error}\n\n")
        if liste_warning:
            print()
            print(stylize('Warning : ', fg('red') + attr('bold') + attr('underlined')))
            print()

            for warning in liste_warning:
                print(warning)
            for error in liste_error:
                print(error)

def sort_output(output, list_kept,):
    """
    This function take output directory path and a list of contigs names for sort all output of this script.
    """
    # Create directory for result
    liste_output_directory = {"ASN_file": ".sqn", "GeneBank_file": ".gbf", "Feature_file": ".tbl",
                              "Fasta_file": ".fasta", "Error_file": ".val"}
    for output_directory in liste_output_directory:
        path_output_directory = Path(f'{output}/{output_directory}')
        if not path_output_directory.exists():
            subprocess.run(f'mkdir {path_output_directory}',
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
        for contig_name in list_kept:
            contig_path = f"{output}/{contig_name}"
            output_path = f'{path_output_directory}/{contig_name}{liste_output_directory[output_directory]}'
            if Path(output_path).exists():
                os.remove(output_path)
            os.rename(f'{contig_path}/{contig_name}{liste_output_directory[output_directory]}', output_path)

    for contig_name in list_kept:
        contig_path = f"{output}/{contig_name}"
        shutil.rmtree(contig_path)


def verify_input(fasta, template):
    """
    This function take GUI intput and verify the format and ext file
    """
    # Verify that fasta can be loaded by Biopython
    try:
        fasta_dict = fasta2dict(fasta)
    except:
        raise TypeError(formated_error(f"The {fasta.split('/')[-1]} isn't correct, "
                                       f"please check that file is at fasta format (maybe duplicates ID) "))
    if fasta_dict == {}:
        raise TypeError(formated_error(f"The {fasta.split('/')[-1]} isn't correct,"
                                       f" please check that file is at fasta format"))
    if not template.endswith('.sbt'):
        raise TypeError(formated_error("The {template} file must have a .sbt extension"))


def print_progress(index, total):
    """
    Function for Gooey API that print progress bar
    """
    print(f"Progress {int((index + 1) / total * 100)} %")
    sys.stdout.flush()


@Gooey(program_name="Submit to GenBank",
       richtext_controls=True,
       use_events=['VALIDATE_FORM'],
       advanced=True,
       program_description="\nCreate ASN file for submission to GenBank",
       progress_regex=r"^Progress (\d+) %$",
       image_dir=Path(__file__).resolve().with_name("image"),
       navigation='SIDEBAR',
       default_size=(999, 999),
       menu=[{
           'name': 'About',
           'items': [{
               'type': 'AboutDialog',
               'menuTitle': 'About',
               'name': 'Gsub - Submit to GenBank',
               'description': 'This tools is used for create ASN.1 file for facilitate submission to GeneBank.'
                              'He use orffinder package python for detect orf and tbl2asn to transformate fasta and '
                              'author template into ASN file',
               'version': '1.0.0',
               'website': 'https://github.com/FlorianCHA/Gsub',
               'developer': 'Florian CHARRIAT',
           }]
       }])
def gui_parser():
    parser = GooeyParser(description='Process some integers.')
    sub = parser.add_subparsers(dest='ssss')
    # Main TABs
    parser = sub.add_parser('Files options', prog="Options")
    parent = parser.add_argument_group('Files options', gooey_options={'columns': 2})
    parent.add_argument('Fasta', widget="FileChooser",
                        help='Path of fasta file that contains all contigs to submit to GenBank\n')
    parent.add_argument('Template', widget="FileChooser",
                        help='Path of template file generate here :\n'
                             'https://submit.ncbi.nlm.nih.gov/genbank/template/submission/\n')
    parent.add_argument('Source', widget="FileChooser",
                        help='Path of Source file')
    parent.add_argument('Output', widget="DirChooser",
                        help='Path to output directory\n')
    filter = parser.add_argument_group('Fasta Filter')
    filter.add_argument('Min_Length_Contig', type=int,
                        default=1500,
                        help='Minimum length for Contig submission\n')
    filter.add_argument('-P', '--Genome', choices=['Eukaryote', 'Prokaryote'],
                        default="Prokaryote")
    orf_parser = filter.add_argument_group('ORF - params', gooey_options={'columns': 2, 'show_border': True})
    orf_parser.add_argument('Min_Length_ORF', type=int,
                            default=300,
                            help='Minimum length for ORF prediction\n', gooey_options={"full_width": False})
    check_box = filter.add_argument_group('ORF - options', gooey_options={'columns': 1, 'show_border': True})
    check_box.add_argument('-O', '--Overlaps', widget="CheckBox",
                           action="store_true",
                           default=False,
                           help='  Keep overlaps genes ')
    check_box.add_argument('-F', '--Frame', widget="CheckBox",
                           action="store_true",
                           default=False,
                           help=' Keep gene with different frame\n')
    pfam_parser = filter.add_argument_group('PFAM option', gooey_options={'show_border': True})
    pfam_parser.add_argument('-s', '--score', type=float,
                             default=50,
                             help='Min Hmm score for predict polymerase\n')
    pfam_parser.add_argument('-e', '--evalue', type=float,
                             default=0.001,
                             help='Evalue cutoff score for predict polymerase\n')
    output_parser = parser.add_argument_group('Assembly information')
    output_parser.add_argument('-a', '--Assembler', type=str,
                               default='',
                               help='Assembler information (as exemple : megahit v. 1.2.9 & Cap3 v. 10.2011)\n')
    output_parser.add_argument('-t', '--Technology', type=str,
                               default='Illumina',
                               help='Sequencing technology (ex : Illumina)\n')

    args = parser.parse_args()
    fasta = args.Fasta
    template = args.Template
    src_file = args.Source
    genome = args.Genome
    output = args.Output
    minlengthcontig = int(args.Min_Length_Contig)
    minlengthorf = int(args.Min_Length_ORF)

    assembler = args.Assembler
    technology = args.Technology

    overlaps = args.Overlaps  # If False, we remove all overlaps gene and keept only the longest
    frame = args.Frame  # If False, we remove all gene with frame different that majority of gene
    # PFAM database for identfication of polymerase
    database = f'{Path(__file__).resolve().parent}/tools/Polymerase.hmm'
    score_cutoff = args.score
    eval_cutoff = args.evalue

    # Cut off for know if orf is only polymerase or if we must give position of polymerase in orf?
    cov_cutoff = 70 / 100

    # Verify format for fasta and the extension of template
    start = time.time()
    print(stylize(f'Gsub tools start at {time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}',
                  fg('green') + attr('bold') + attr('underlined')))
    print()
    verify_input(fasta, template)

    # Launch tbl2asn & orffinder
    list_kept = create_submit_file(fasta, template, src_file, output, minlengthorf, minlengthcontig, assembler,
                                   technology, overlaps, frame, database, score_cutoff, eval_cutoff, cov_cutoff, genome)

    # Contigs name contains all Id of seq in fasta input
    contigs_names = list(fasta2dict(fasta).keys())

    print()
    end = time.time()
    print(stylize(f'Gsub tools finish at {time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}',
                  fg('green') + attr('bold') + attr('underlined')))
    print(stylize(f'Time elapsed : {round(end - start, 3)} sec', fg('black') + attr('bold')))

    # Verify if all errorsummary.val file are empty
    verif_quality(output, contigs_names, list_kept)

    # Sort all output file
    sort_output(output, list_kept)


if __name__ == '__main__':
    gui_parser()
