import os
import gzip
import sys
from Bio import Seq
import multiprocessing
import itertools
import argparse
import pandas as pd
import string

def writeToCounts(fileTup):
    print('\t'.join(fileTup[:-2]))
    sys.stdout.flush()

    mismatchDicts = fileTup[-2]
    testRun = fileTup[-1]

    statsCounts = {'A sgRNA not mapped': 0,
    'B sgRNA not mapped': 0, 
    'A sgRNA multiple mappings': 0,
    'B sgRNA multiple mappings': 0, 
    'All sgRNAs uniquely map': 0,
    'A sgRNA and B sgRNA do not match': 0,
    'A sgRNA and B sgRNA match': 0}

    pairCounts_sgRNAs = dict()
    pairCounts_double = dict()

    with gzip.open(fileTup[0]) as infile_r1:
        with gzip.open(fileTup[1]) as infile_r2:
                for i, (r1,r2) in enumerate(zip(infile_r1, infile_r2)):
                    if i%4 == 1:
                        
                        ##trimming is adjusted based on sequencing strategy
                        protospacer_a_r1 = matchBarcode(mismatchDicts['protospacer_a_r1'], r1[1:20].strip().decode("utf-8"), allowOneMismatch = False)
 
                        protospacer_b_r2 = matchBarcode(mismatchDicts['protospacer_b_r2'], r2[:19].strip().decode("utf-8"), allowOneMismatch = False)

                        if protospacer_a_r1 == 'none' \
                            or protospacer_b_r2 == 'none' \
                            or protospacer_a_r1 == 'multiple' \
                            or protospacer_b_r2 == 'multiple':

                            if protospacer_a_r1 == 'none':
                                statsCounts['A sgRNA not mapped'] += 1
                            if protospacer_b_r2 == 'none':
                                statsCounts['B sgRNA not mapped'] += 1
                            if protospacer_a_r1 == 'multiple':
                                statsCounts['A sgRNA multiple mappings'] += 1
                            if protospacer_b_r2 == 'multiple':
                                statsCounts['B sgRNA multiple mappings'] += 1

                        else:
                            statsCounts['All sgRNAs uniquely map'] += 1

                            combinedSgId = protospacer_a_r1 + '++' + protospacer_b_r2
                            if combinedSgId not in pairCounts_sgRNAs:
                                pairCounts_sgRNAs[combinedSgId] = 0

                            pairCounts_sgRNAs[combinedSgId] += 1

                            if protospacer_a_r1 != protospacer_b_r2 :
                                statsCounts['A sgRNA and B sgRNA do not match'] += 1
                                
                            else:
                                statsCounts['A sgRNA and B sgRNA match'] += 1

                                combinedSgId = protospacer_a_r1 + '++' + protospacer_b_r2

                                if combinedSgId not in pairCounts_double:
                                    pairCounts_double[combinedSgId] = 0

                                pairCounts_double[combinedSgId] += 1

                    if testRun and i == testLines:
                        break

    with open(fileTup[2] + '.all.aligned.counts.txt', 'w') as outfile:
        for pair in sorted(pairCounts_sgRNAs.keys()):
            outfile.write(pair + '\t' + str(pairCounts_sgRNAs[pair]) + '\n')

    with open(fileTup[2] + '.AB.match.counts.txt', 'w') as outfile:
        for pair in sorted(pairCounts_double.keys()):
            outfile.write(pair + '\t' + str(pairCounts_double[pair]) + '\n')

    numReads = (i+1)/4
    print(fileTup[2], numReads, 'reads')

    print('Percent A sgRNAs mapping', 100.0 - (statsCounts['A sgRNA not mapped'] * 100.0 / numReads))
    print('Percent B sgRNAs mapping', 100.0 - (statsCounts['B sgRNA not mapped'] * 100.0 / numReads))
    print('Percent all sgRNAs mapping', statsCounts['All sgRNAs uniquely map'] * 100.0 / numReads)
    print('Percent A sgRNA and B sgRNA mismatch', statsCounts['A sgRNA and B sgRNA do not match'] * 100.0 / statsCounts['All sgRNAs uniquely map'])
    print('Percent both A and B match', statsCounts['A sgRNA and B sgRNA match'] * 100.0 / statsCounts['All sgRNAs uniquely map'])
    print('Total percent matching and mapping reads', statsCounts['A sgRNA and B sgRNA match'] * 100.0 / numReads)
    
    sys.stdout.flush()


def getMismatchDict(table, column, trim_range = None, allowOneMismatch = True):
    mismatch_dict = dict()

    if trim_range:
        col = table[column].apply(lambda seq: seq[trim_range[0]:trim_range[1]])
    else:
        col = table[column]

    for sgRNA, seq in col.iteritems():
        if seq in mismatch_dict:
            print('clash with 0 mismatches', sgRNA, seq)
            mismatch_dict[seq] = 'multiple'
        else:
            mismatch_dict[seq] = sgRNA

        if allowOneMismatch:
            for position in range(len(seq)):
                mismatchSeq = seq[:position] + 'N' + seq[position + 1:]

                if mismatchSeq in mismatch_dict:
                    print('clash with 1 mismatch', sgRNA, mismatchSeq)
                    mismatch_dict[seq] = 'multiple'
                else:
                    mismatch_dict[mismatchSeq] = sgRNA

    return mismatch_dict

def matchBarcode(mismatch_dict, barcode, allowOneMismatch = True):
    if barcode in mismatch_dict:
        match = mismatch_dict[barcode]

    elif allowOneMismatch:
        match = 'none'
        for position in range(len(barcode)):
            mismatchSeq = barcode[:position] + 'N' + barcode[position + 1:]

            if mismatchSeq in mismatch_dict:
                if match == 'none':
                    match = mismatch_dict[mismatchSeq]
                else:
                    match = 'multiple'

    else:
        match = 'none'

    return match

testLines = 100000

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process raw sequencing data from screens to counts files in parallel.')
    parser.add_argument('Guide_Table', help='Table of sgRNA pairs in the library.')
    parser.add_argument('Out_File_Path', help='Directory where output files should be written.')
    parser.add_argument('Seq_File_Names', nargs='+', help='Name(s) of sequencing file(s). Unix wildcards can be used to select multiple files at once. The script will search for all *.fastq.gz, *.fastq, and *.fa(/fasta/fna) files with the given wildcard name.')
    parser.add_argument('--test', action='store_true', default=False, help='Run the entire script on only the first %d reads of each file. Be sure to delete or move all test files before re-running script as they will not be overwritten.' % testLines)

    args = parser.parse_args()

    outputDirectory = args.Out_File_Path

    inputFileList = sorted(args.Seq_File_Names)

    guideTable = pd.read_csv(args.Guide_Table, sep=',', header=0)[['sgID_AB','sgID_A','protospacer_A','sgID_B','protospacer_B']]
    guideTable = guideTable.set_index('sgID_AB').dropna()
    ##convert protospacers to mimick reads according to sequencing strategy
    ##assuming R1=19, R2=19
    guideTable['protospacer_A']=guideTable['protospacer_A'].str[1:20].str.upper()
    trans=str.maketrans('ATGC', 'TACG')
    guideTable['protospacer_B']=guideTable['protospacer_B'].str[1:20].str.upper().str.translate(trans).str[::-1]
    print('sgRNAs in library', len(guideTable))
    
    combinedMismatchDicts = {'protospacer_a_r1': getMismatchDict(guideTable, 'protospacer_A', allowOneMismatch=False),
    'protospacer_b_r2': getMismatchDict(guideTable, 'protospacer_B', allowOneMismatch=False),
    }

    fileTups = []
    for i, fastqfile in enumerate(inputFileList):
        if i%2 == 0:
            r1file = fastqfile
        elif i%2 == 1:
            r2file = fastqfile
            outputfile = os.path.join(outputDirectory, os.path.split(fastqfile)[-1].split('_R')[0])
            fileTups.append((r1file, r2file, outputfile, combinedMismatchDicts, args.test))

    pool = multiprocessing.Pool(len(fileTups))

    pool.map(writeToCounts, fileTups)

    pool.close()
    pool.join()
