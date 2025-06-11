#!/bin/env python

import sys
import os
import argparse
import numpy as np
import time

from adpred import ADpred as adp

aa_one_nat = "ACDEFGHIKLMNPQRSTVWY"

def is_aa_one_nat(sequence, additional=""):
    for aa in sequence.upper():
        if not (aa in aa_one_nat or aa in additional.upper()):
            return(False)
    return(True)


def read_fasta(filename, comment_char="#;", extra_aa="", skip_nonnatural=False, verbose=0):
    """Flexible FASTA file reader without dependencies"""
    seq_list = []
    file_id = filename.split("/")[-1].split(".")[0]
    reading_fasta = False
    with open(filename,"r") as file_handle:
        for line in file_handle:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0][0] in comment_char:
                continue
            # start reading new fasta record                                                                                                                                                                            
            if words[0][0] == ">":
                if reading_fasta:
                    if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                        seq_list.append((seq_id,seq))
                    else:
                        if verbose > 0:
                            print("Skipping protein with non-standard amino acids: %s" % (seq_id), file=sys.stderr)
                seq_id = words[0][1:]
                seq = ""
                reading_fasta = True
                # continue reading fasta record
            elif reading_fasta:
                if len(words) > 1:
                    raise ValueError("Found FASTA line with more white space separated fields: '%s'" % (line.strip()))
                seq = seq + words[0]
            # read table-format (.seq) with single-column sequences or multi-column with sequences in second column                                                                                                     
            else:
                if len(words) == 1:
                    seq = words[0]
                    if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                        seq_id = file_id+"%05d" % (len(seq_list))
                        seq_list.append((seq_id,seq))
                    else:
                        if verbose > 0:
                            print("Skipping line in non-FASTA single-column file with non-natural amino acid sequence in first column:", file=sys.stderr)
                            print("         %s" % (line.strip()), file=sys.stderr)
                elif len(words) > 1:
                    seq = words[1]
                    if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                        seq_list.append((words[0],seq))
                    else:
                        print("WARNING: Skipping line in non-FASTA multi-column file with non-natural amino acid sequence in second column:", file=sys.stderr)
                        print("         %s" % (line.strip()), file=sys.stderr)
        # append the last fasta record                                                                                                                                                                                  
        if reading_fasta:
            if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                seq_list.append((seq_id,seq))
            else:
                print("WARNING: Skipping protein with non-standard amino acids: %s" % (seq_id), file=sys.stderr)
    return(seq_list)

if __name__ == "__main__":
    # argument parser
    arg_parser = argparse.ArgumentParser(description="Predict degradation/abundance scores of proteins and peptides assuming these are solvent/PQC exposed in the cell")
    # keyword arguments
    # arg_parser.add_argument("-s", "--tile-size", type=int, default=None,
    #                         help="Number of amino acids per tile. Set to zero to keep sequence lengths as inputed. Default depends on the model used")
    # arg_parser.add_argument("-c","--ignore-ct-term", default=False, action='store_true',
    #                         help="By default, C-degron terms (if available in model) are applied to the last tile of each sequence")
    # arg_parser.add_argument("-t", "--term-method", choices=["full-tiles", "full-coverage"], default="full-tiles",
    #                         help="Method to handle terminals. Either considers tiles of full length only (default) or cover the entire sequence "+
    #                         "by evaluating tiles down to half length. If length is even the central amino acid is the one to the right of the center")
    # arg_parser.add_argument('--saturation-mut', dest='saturation', default=False, action='store_true',
    #                         help="Perform satuation mutagenesis of the central amino acid in each tile. The output will typically have 20 lines per tile"+
    #                         "with the central amino acid mutated. If tile length is even the central amino acid is the one to the right of the center")
    arg_parser.add_argument('--skip-nonnatural', default=False, action='store_true',
                            help="Ignore tiles that contains other then the 20 standard amino acids")
    # arg_parser.add_argument("--output-format", choices=["std", "prism"], default="std",
    #                         help="Output to stdout or write prism files in directory 'prism_pap'")
    arg_parser.add_argument('-v', '--verbose', action='count', default=0)
    # positional argument
    arg_parser.add_argument("seq_input", nargs="+", metavar="SEQ",
                            help="Sequence(s) in text or FASTA format. More inputs should by whitespace separated")
    args = arg_parser.parse_args()

    # extra_aa="XU"
    extra_aa=""

    # read all sequences into memory
    input_list = []
    for seq_input in args.seq_input:
        if is_aa_one_nat(seq_input.upper(), extra_aa):
            input_list.append(("seq%05d" % (len(input_list)), seq_input.upper()))
        elif seq_input[-4:].lower() in [".seq",".fas"] or seq_input[-6:].lower() == ".fasta":
            if not os.path.isfile(seq_input):
                raise ValueError("ERROR: Cannot find file %s" % (seq_input))
            input_list.extend(read_fasta(seq_input, extra_aa=extra_aa, skip_nonnatural=args.skip_nonnatural, verbose=args.verbose))
        else:
            raise ValueError("ERROR: Argument %s seems to be neither a protein sequence nor a FASTA file (.seq, .fas or .fasta)" % (seq_input))
            
    print("Loaded %d input sequences" % (len(input_list)), file=sys.stderr)

    out_fn = "adpred_scores.txt"
    completed_prot = []

    # If present, read output file for compleated predictions
    if os.path.isfile(out_fn):
        with open(out_fn, "rt") as file_handle:
            for line in file_handle:
                fields = line.split()
                # print(fields, file=sys.stderr)
                if len(fields) == 3 and fields[0][0] != "#":
                    completed_prot.append(fields[0])
        print("Found predictions of %d sequences in file %s" % (len(completed_prot),out_fn), file=sys.stderr)

    # Predict all proteins
    outfile = open(out_fn, "at")
    fails_ss = 0
    fails_adpred = 0
    for si in range(len(input_list)):
        (aaseq_name,aaseq) = input_list[si]

        # Make protein object
        adpred_prot = adp.protein(prot_id=aaseq_name, sequence=aaseq, second_struct=None)
        print("Protein %04d: %s of %d residues" % (si, adpred_prot.prot_id, len(adpred_prot.sequence)), file=sys.stderr)

        # Check for existing prediction
        if adpred_prot.prot_id in completed_prot:
            print("    protein %s has already been predicted in output file %s" % (adpred_prot.prot_id,out_fn), file=sys.stderr)
            continue
        
        # Predict secondary structure features as string of ['E','H','-']
        t0 = time.time()
        try:
            adpred_prot.predict_second_struct()
        except Exception as err:
            print("    ERROR: secondary structure prediction failed with error:", file=sys.stderr)
            print("    "+str(err), file=sys.stderr)
            fails_ss += 1
            continue
        print("    secondary structure prediction took %.1f min" % ((time.time()-t0)/60.0), file=sys.stderr)

        # Predict ADPred scores
        t0 = time.time()
        # should not predict ss if already available in adpred_prot.second_struct
        try:
            adpred_prot.predict()
        except Exception as err:
            print("    ERROR: ADpred prediction failed with error:", file=sys.stderr)
            print("    "+str(err), file=sys.stderr)
            fails_adpred += 1
            continue
        print("    adpred took %.1f min" % ((time.time()-t0)/60.0), file=sys.stderr)

        # Output
        scores_str = ";".join(["%.4f" % (score) for score in adpred_prot.predictions])
        outfile.write("%s   %s   %s\n" % (adpred_prot.prot_id, adpred_prot.second_struct, scores_str))
        outfile.flush()
        completed_prot.append(adpred_prot.prot_id)

    outfile.close()
    if fails_ss+fails_adpred > 0:
        print("Done - there were %d fails of SS prediction and %d fails of ADpred" % (fails_ss,fails_adpred), file=sys.stderr)
    else:
        print("Done - all good", file=sys.stderr)
