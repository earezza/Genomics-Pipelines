#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:12:54 2022
Description:
    Script to run the Snakemake pipeline file without snakemake (essentially a copy of executables).
    Processes CUT&TAG sequence data to generate bigwig/bigbed files for viewing in the UCSC Browser.
    
    Reads and md5sum.txt must be saved into their own directory.
    
    Requires input files formatted as {SAMPLE_NAME}_{READ}.fastq.gz
    e.g.
        MY-SAMPLE-NAME_R1.fastq.gz
        MY-SAMPLE-NAME_R2.fastq.gz
        
    An md5 hash for each read must also be saved into a file named md5sum.txt
    e.g. md5sum.txt would look like:
        abcdefghijklmnop MY-SAMPLE-NAME_R1.fastq.gz
        qrstuvqxyz123456 MY-SAMPLE-NAME_R2.fastq.gz
    
    Output from the pipeline will generate a file tree as follows:
        logs/
        Analysis_Results/
        All_output/
    
    
@author: earezza
"""

__all__ = [ 
    'md5sum_check',
    'QCrawreads_Fastqc',
    'Compileresults_QC',
    'AdapterTrim_Cutadapt',
    'Compileresults_PosttrimQC',
    'Map_Bowtie2',
    'Collect_alignment_stats',
    'Compileresults_map',
    'Filtering_bams_PicardSamtools',
    'Compileresults_filtering',
    'GetBigwigs_BamCoverage',
    'Map2Spikein_Bowtie2',
    'Collect_Spikealignment_stats',
    'Compileresults_Spike',
    'CalcNormFactors',
    'GetNormBwsBdgs_BamCoverage',
    'Peak_Calling',
    'Clean_up',
            ]

import os
import sys
import time
import types
import argparse
import pandas as pd
import numpy as np
import logging
import subprocess
import glob
import yaml
from pathlib import Path

describe_help = 'python run_CnT_pipeline.py --reads READS_DIR/ --species Mus musculus --spikein Amp -l 50 --samples ExampleSample1 ExampleSample2 -c Target IgG'
parser = argparse.ArgumentParser(description=describe_help)
# User defined options
parser.add_argument('-logfile', '--logfile', help='Name of .log file', type=str, default="cnt_pipeline.log")
parser.add_argument('-outdir', '--outdir', help='Path of directory to store results', type=str, default="./")
parser.add_argument('-reads', '--reads', help='Path to directory containing reads (R1, R2, and md5sum.txt files)', type=str)
parser.add_argument('-species', '--species', help='Species of reads (Mus for mouse, Homo for human)', type=str, choices=['Mus', 'Homo'], default='Mus')
parser.add_argument('-spikein', '--spikein', help='Spikein type', type=str, choices=['Amp', 'Bacteria'], default='Amp')
parser.add_argument('-length', '--length', help='Read length', type=str, default='100', choices=['50', '75', '100', '150', '200'])
parser.add_argument('-controls', '--controls', help='Control reads for peaks calling', default=[], nargs='+')
parser.add_argument('-adapters', '--adapters', help='Adapter sequences, see https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/Nextera/SequencesNXTILMPrepAndPCR.htm',
                    type=int, default=1, choices=[1, 2, 3, 4, 5])
parser.add_argument('-no_spikein', '--no_spikein', help='If no spikein, skip steps for normalizing to spikein', action='store_true')
parser.add_argument('-cleanup', '--cleanup', help='If cleanup, remove all intermediate files keeping only final .bw and .bed files', action='store_true')
# Program and reference genome locations
parser.add_argument('-PicardLoc', '--PicardLoc', help='Location of picard.jar', type=str, default="java -jar /home/earezza/projects/def-jdilwort/earezza/picard.jar")
parser.add_argument('-SEACRLoc', '--SEACRLoc', help='Location of SEACR .sh', type=str, default="/home/earezza/projects/def-jdilwort/earezza/SEACR/SEACR_1.3.sh")
parser.add_argument('-bowtie2_index', '--bowtie2_index', help='Location of bowtie2 genome index files (from iGenomes)', type=str, default="/home/earezza/projects/def-jdilwort/earezza/CnT_pipeline_snakemake/Reference_files/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome")
# Usually unchanged command line options for programs
parser.add_argument('-spike_align', '--spike_align', help='Command input options for spikein alignment', type=str, default="-p 8 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700")
parser.add_argument('-bamCov_default', '--bamCov_default', help='Command input default for bamCoverage', type=str, default="--binSize 1 --ignoreForNormalization 'chrM' --extendReads")
parser.add_argument('-bamCov_min', '--bamCov_min', help='Command input options for bamCoverage', type=str, default="--binSize 1 --extendReads")
parser.add_argument('-bamCov_RPGC', '--bamCov_RPGC', help='Command input options for bamCoverage, reads per genomic content (1x normalization)', type=str, default="--binSize 1 --normalizeUsing 'RPGC'  --ignoreForNormalization 'chrM' --extendReads")
parser.add_argument('-bamCov_CPM', '--bamCov_CPM', help='Command input options for bamCoverage, number of reads per bin / number of mapped reads (in millions)', type=str, default="--binSize 10 --normalizeUsing 'CPM'  --ignoreForNormalization 'chrM' --extendReads") # Not used here
parser.add_argument('-genome_align', '--genome_align', help='Command input options for genome alignment', type=str, default="-p 8 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
parser.add_argument('-samtools_mapq', '--samtools_mapq', help='Command input option for samtools mapq', type=str, default="-q 10")
parser.add_argument('-samtools_proper_paired', '--samtools_proper_paired', help='Command input option for samtools', type=str, default="-f 2")
parser.add_argument('-mnase', '--mnase', help='If MNase-seq data, apply --MNase bam coverage option', action='store_true')

# Spikein reference index files
parser.add_argument('-spikein_index_amp', '--spikein_index_amp', help="Source of reference files for Amp spikein index file", type=str, default='/home/earezza/projects/def-jdilwort/earezza/CnT_pipeline_snakemake/Reference_files/Spikein_indices/Amp_pbluescript/Amp_index/Amp_pBlue')
parser.add_argument('-spikein_index_Ecoli', '--spikein_index_Ecoli', help="Source of reference files for Amp spikein index file", type=str, default='/home/earezza/projects/def-jdilwort/earezza/CnT_pipeline_snakemake/Reference_files/Spikein_indices/EcoliK12_index/EcoliK12Index/EcoliK12')
args = parser.parse_args()

# Define constants
ROOT_DIR = os.getcwd()
READS_DIR = args.reads
OUT_DIR = args.outdir if args.outdir[-1] == '/' else args.outdir+'/'
if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)
    
if args.mnase:
    args.bamCov_default = args.bamCov_default + ' --MNase'
    args.bamCov_min = args.bamCov_min + ' --MNase'
    args.bamCov_RPGC = args.bamCov_RPGC + ' --MNase'
    
# Assumes input reads are paired-end...
illumina_adapter_sequences = {
    1: '-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT',
    2: '-a CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA -A CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA',
    3: '-a CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG',
    4: '-g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG',
    5: '-g CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG -G AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC', # see illumina reference docs for indices...
    }
ADAPTERS = illumina_adapter_sequences[args.adapters]

global read_files
read_files = [ f for f in os.listdir(args.reads) if '.fastq.gz' in f and ('_R1' in f or '_R2' in f) and '.md5' not in f ]
read_files.sort()
global reads
reads = set([ r.replace('.fastq.gz', '').split('_')[-1] for r in read_files ])
global fastqfiles
fastqfiles = set([ r.replace('.fastq.gz', '').split('_')[0] for r in read_files ])

# Control vs non-Control files
IGGREADS = set(args.controls)
TARGETS = set([ f for f in fastqfiles if f not in IGGREADS ])

# Read lengths genome mapping
# https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
EGS_GRCh38 = {'50': '2701495761', '75': '2747877777', '100': '2805636331', '150': '2862010578', '200': '2887553303'}
EGS_GRCm38 = {'50': '2308125349', '75': '2407883318', '100': '2467481108', '150': '2494787188', '200': '2520869189'}
if args.species == 'Mus':
    EFFECTIVEGENOMESIZE = EGS_GRCm38[args.length]
if args.species == 'Homo':
    EFFECTIVEGENOMESIZE = EGS_GRCh38[args.length]

# Reference index file for spikein
if args.spikein == 'Amp':
    SPIKEINDEX = args.spikein_index_amp
if args.spikein == 'Bacteria':
    SPIKEINDEX = args.spikein_index_Ecoli


# Step 1
def md5sum_check():
    passed = True
    logger = logging.getLogger('md5sum_check')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/md5checks.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    os.chdir(READS_DIR)
    logger.info('Looking for md5sum.txt file...')
    if os.path.exists('md5sum.txt'):
        logger.info('Found!')
    else:
        logger.info('Not found, formatting reads files and creating md5sum.txt...')
        files = [ f for f in os.listdir() if 'fastq.gz' in f and ('_R1' in f or '_R2' in f)]
        files.sort()
        reads_files = files[0:len(files):2]
        md5 = files[1:len(files):2]
        if len(md5) != len(reads_files):
            logger.error('Check fastq.gz reads and corresponding .md5 files...')
            os.chdir(ROOT_DIR)
            return False
        md5sumfile = pd.DataFrame(columns=['md5', 'read'])
        # Iterate over each raw read file
        for i in range(len(reads_files)):
            f = reads_files[i]
            m = md5[i]
            r = f.split('_')[:-1][-1]
            f_formatted = '-'.join(f.split('_')[:-2])
            f_formatted = f_formatted + '_' + r + '.fastq.gz'
            # Rename reads file
            subprocess.run('mv %s %s'%(f, f_formatted), shell=True, capture_output=False, text=True)
            # Read md5 keys and write to file for pipeline input
            key = pd.read_csv(m, delim_whitespace=True).columns[0]
            md5sumfile = pd.concat([md5sumfile, pd.DataFrame.from_records([{'md5': key, 'read': f_formatted}])])
        # Write md5sum.txt
        md5sumfile.drop_duplicates(inplace=True)
        md5sumfile.to_csv('md5sum.txt', sep=' ', header=None, index=False)
        os.chdir(ROOT_DIR)
        global read_files
        read_files = [ f for f in os.listdir(args.reads) if '.fastq.gz' in f and ('_R1' in f or '_R2' in f) and '.md5' not in f ]
        read_files.sort()
        global reads
        reads = set(np.array([ [ t for t in r.replace('.fastq.gz', '').split('_') if t[0] == 'R' ] for r in read_files ]).flatten())
        global fastqfiles
        fastqfiles = set([ r.replace('.fastq.gz', '').split('_')[0] for r in read_files ])
        os.chdir(READS_DIR)
    try:
        result = subprocess.run('md5sum --check md5sum.txt', shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
        if set(result.stdout.split()[1::2]) != {'OK'}:
            passed = False
    except Exception as e:
        logging.exception(e)
        passed = False
    os.chdir(ROOT_DIR)
    
    return passed

# Step 2
def QCrawreads_Fastqc():
    if not os.path.exists(OUT_DIR+'Analysis_Results'):
        os.mkdir(OUT_DIR+'Analysis_Results')
    if not os.path.exists(OUT_DIR+'Analysis_Results/QC_Rawreads'):
        os.mkdir(OUT_DIR+'Analysis_Results/QC_Rawreads')
    if not os.path.exists(OUT_DIR+'logs/fastqc_rawreads'):
        os.mkdir(OUT_DIR+'logs/fastqc_rawreads')
    passed = True
    for f in fastqfiles:
        for r in reads:
            logger = logging.getLogger(OUT_DIR+'QCrawreads_Fastqc')
            logger.setLevel(logging.DEBUG)
            file_handler = logging.FileHandler(OUT_DIR+'logs/fastqc_rawreads/%s_%s.log'%(f, r))
            formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            if os.path.exists(OUT_DIR+'Analysis_Results/QC_Rawreads/%s_%s_fastqc.html'%(f,r)) and os.path.exists(OUT_DIR+'Analysis_Results/QC_Rawreads/%s_%s_fastqc.zip'%(f,r)):
                continue
            try:
                result = subprocess.run(('fastqc %s/%s_%s.fastq.gz --outdir=%sAnalysis_Results/QC_Rawreads'%(READS_DIR, f, r, OUT_DIR)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
            except Exception as e:
                logger.exception(e)
                passed = False
            passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/QC_Rawreads/%s_%s_fastqc.html'%(f,r))
            passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/QC_Rawreads/%s_%s_fastqc.zip'%(f,r))
    return passed

# Step 3
def Compileresults_QC():
    if not os.path.exists(OUT_DIR+'logs/compileresults'):
        os.mkdir(OUT_DIR+'logs/compileresults')
    passed = True
    logger = logging.getLogger(OUT_DIR+'Compileresults_QC')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/compileresults/QC.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists(OUT_DIR+'Analysis_Results/QC_Rawreads/Rawreads_QC_%s.html'%(args.logfile.rstrip('.log'))):
        return passed
    try:
        result = subprocess.run('multiqc %sAnalysis_Results/QC_Rawreads --force -v -o %sAnalysis_Results/QC_Rawreads -n Rawreads_QC_%s.html'%(OUT_DIR, OUT_DIR, args.logfile.rstrip('.log')), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    
    # Get read length if differs from args.length
    try:
        multiqc = open('%sAnalysis_Results/QC_Rawreads/Rawreads_QC_%s.html'%(OUT_DIR, args.logfile.rstrip('.log')), 'r')
        content = multiqc.read()
        multiqc.close()
        global READ_LENGTH
        read_length = content[content.find('length (')+8:content.find('bp).')]
        logger.info('MultiQC - Read Length: %s'%read_length)
        read_length = str(25*round(int(read_length)/25))
        READ_LENGTH = int(read_length)
        if read_length != args.length:
            global EFFECTIVEGENOMESIZE
            if args.species == 'Mus':
                EFFECTIVEGENOMESIZE = EGS_GRCm38[read_length]
            if args.species == 'Homo':
                EFFECTIVEGENOMESIZE = EGS_GRCh38[read_length]
    except Exception as e:
        logger.exception(e)
        passed = False
    
    passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/QC_Rawreads/Rawreads_QC_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 4
def AdapterTrim_Cutadapt():
    if not os.path.exists(OUT_DIR+'All_output'):
        os.mkdir(OUT_DIR+'All_output')
    if not os.path.exists(OUT_DIR+'All_output/Trimmed_reads'):
        os.mkdir(OUT_DIR+'All_output/Trimmed_reads')
    if not os.path.exists(OUT_DIR+'logs/cutadapt'):
        os.mkdir(OUT_DIR+'logs/cutadapt')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger(OUT_DIR+'AdapterTrim_Cutadapt')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/cutadapt/%s.log'%(f))
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists(OUT_DIR+'All_output/Trimmed_reads/%s_Trimmed_R1.fastq'%(f)) and os.path.exists(OUT_DIR+'All_output/Trimmed_reads/%s_Trimmed_R2.fastq'%(f)):
            continue
        try:
            result = subprocess.run(('cutadapt %s -o %sAll_output/Trimmed_reads/%s_Trimmed_R1.fastq -p %sAll_output/Trimmed_reads/%s_Trimmed_R2.fastq %s_R1.fastq.gz %s_R2.fastq.gz'%(ADAPTERS, OUT_DIR, f, OUT_DIR, f, READS_DIR+'/'+f, READS_DIR+'/'+f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists(OUT_DIR+'All_output/Trimmed_reads/%s_Trimmed_R1.fastq'%(f))
        passed = passed and os.path.exists(OUT_DIR+'All_output/Trimmed_reads/%s_Trimmed_R2.fastq'%(f))
    return passed

# Step 5
def Compileresults_PosttrimQC():
    if not os.path.exists(OUT_DIR+'Analysis_Results/Trimming'):
        os.mkdir(OUT_DIR+'Analysis_Results/Trimming')
    passed = True
    logger = logging.getLogger(OUT_DIR+'Compileresults_PosttrimQC')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/compileresults/PosttrimQC.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists(OUT_DIR+'Analysis_Results/Trimming/PostTrimming_QC_%s.html'%args.logfile.rstrip('.log')):
        return passed
    try:
        result = subprocess.run('multiqc %slogs/cutadapt --force -v -o %sAnalysis_Results/Trimming -n PostTrimming_QC_%s.html'%(OUT_DIR, OUT_DIR, args.logfile.rstrip('.log')), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Trimming/PostTrimming_QC_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 6
def Map_Bowtie2():
    if not os.path.exists(OUT_DIR+'All_output/Mapped_reads'):
        os.mkdir(OUT_DIR+'All_output/Mapped_reads')
    if not os.path.exists(OUT_DIR+'logs/primary_alignment'):
        os.mkdir(OUT_DIR+'logs/primary_alignment')
    if not os.path.exists(OUT_DIR+'logs/primary_alignment/bowtie2'):
        os.mkdir(OUT_DIR+'logs/primary_alignment/bowtie2')
    if not os.path.exists(OUT_DIR+'logs/primary_alignment/picard_sort'):
        os.mkdir(OUT_DIR+'logs/primary_alignment/picard_sort')
    passed = True
    for f in fastqfiles:
        # Bowtie2
        logger = logging.getLogger(OUT_DIR+'Map_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/primary_alignment/bowtie2/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists(OUT_DIR+'All_output/Mapped_reads/%s.bam'%f) and os.path.exists(OUT_DIR+'All_output/Mapped_reads/%s.coordsorted.bam'%f) and os.path.exists(OUT_DIR+'All_output/Mapped_reads/%s.coordsorted.bam.bai'%f):
            continue
        try:
            result = subprocess.run(('bowtie2 %s -x %s -1 %sAll_output/Trimmed_reads/%s_Trimmed_R1.fastq -2 %sAll_output/Trimmed_reads/%s_Trimmed_R2.fastq 2> %slogs/primary_alignment/bowtie2/%s.log | samtools view -bS - > %sAll_output/Mapped_reads/%s.bam'%(args.genome_align, args.bowtie2_index, OUT_DIR, f, OUT_DIR, f, OUT_DIR, f, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.error(e)
            passed = False
        # Picard
        logger = logging.getLogger(OUT_DIR+'Map_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/primary_alignment/picard_sort/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('%s SortSam -I %sAll_output/Mapped_reads/%s.bam -O %sAll_output/Mapped_reads/%s.coordsorted.bam -SORT_ORDER coordinate'%(args.PicardLoc, OUT_DIR, f, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Samtools
        try:
            result = subprocess.run(('samtools index %sAll_output/Mapped_reads/%s.coordsorted.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists(OUT_DIR+'All_output/Mapped_reads/%s.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'All_output/Mapped_reads/%s.coordsorted.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'All_output/Mapped_reads/%s.coordsorted.bam.bai'%f)
    return passed    
        
# Step 7
def Collect_alignment_stats():
    if not os.path.exists(OUT_DIR+'logs/primary_alignment/PostAlignmentStats'):
        os.mkdir(OUT_DIR+'logs/primary_alignment/PostAlignmentStats')
    if not os.path.exists(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/dupstats'):
        os.mkdir(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/dupstats')
    if not os.path.exists(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/flagstat'):
        os.mkdir(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/flagstat')
    passed = True
    for f in fastqfiles:
        # Picard
        logger = logging.getLogger(OUT_DIR+'Collect_alignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/dupstats/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/dupstats/%s.dupMarked.bam'%f) and os.path.exists(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/dupstats/%s_picard.dupMark.txt'%f):
            continue
        try:
            result = subprocess.run(('%s MarkDuplicates -I %sAll_output/Mapped_reads/%s.coordsorted.bam -O %slogs/primary_alignment/PostAlignmentStats/dupstats/%s.dupMarked.bam -METRICS_FILE %slogs/primary_alignment/PostAlignmentStats/dupstats/%s_picard.dupMark.txt'%(args.PicardLoc, OUT_DIR, f, OUT_DIR, f, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Samtools
        logger = logging.getLogger(OUT_DIR+'Collect_alignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/flagstat/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('samtools flagstat %sAll_output/Mapped_reads/%s.coordsorted.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/dupstats/%s.dupMarked.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'logs/primary_alignment/PostAlignmentStats/dupstats/%s_picard.dupMark.txt'%f)
    return passed

# Step 8
def Compileresults_map():
    passed = True
    logger = logging.getLogger(OUT_DIR+'Compileresults_map')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/compileresults/map.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists(OUT_DIR+'Analysis_Results/primary_alignment/Alignment_results_%s.html'%args.logfile.rstrip('.log')):
        return passed
    try:
        result = subprocess.run('multiqc %slogs/primary_alignment --force -v -o %sAnalysis_Results/primary_alignment -n Alignment_results_%s.html'%(OUT_DIR, OUT_DIR, args.logfile.rstrip('.log')), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/primary_alignment/Alignment_results_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 9
def Filtering_bams_PicardSamtools():
    if not os.path.exists(OUT_DIR+'logs/filtered_bams'):
        os.mkdir(OUT_DIR+'logs/filtered_bams')
    if not os.path.exists(OUT_DIR+'All_output/Processed_reads'):
        os.mkdir(OUT_DIR+'All_output/Processed_reads')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger(OUT_DIR+'Filtering_bams_PicardSamtools')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/filtered_bams/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if (OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%f) and os.path.exists(OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam.bai'%f) and os.path.exists(OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%f) and os.path.exists(OUT_DIR+'logs/filtered_bams/%s_picard.rmDup.txt'%f) and os.path.exists(OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam.bai'%f):
            continue
        try:
            result = subprocess.run(('samtools view -bu %s %sAll_output/Mapped_reads/%s.coordsorted.bam | samtools view -b %s - | samtools sort - -o %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%(args.samtools_proper_paired, OUT_DIR, f, args.samtools_mapq, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('%s MarkDuplicates -I %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.bam -O %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -REMOVE_DUPLICATES true -METRICS_FILE %slogs/filtered_bams/%s_picard.rmDup.txt'%(args.PicardLoc, OUT_DIR, f, OUT_DIR, f, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('samtools flagstat %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('samtools index %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('samtools index %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
    
        passed = passed and os.path.exists(OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam.bai'%f)
        passed = passed and os.path.exists(OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'logs/filtered_bams/%s_picard.rmDup.txt'%f)
        passed = passed and os.path.exists(OUT_DIR+'All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam.bai'%f)
    return passed

# Step 10
def Compileresults_filtering():
    passed = True
    logger = logging.getLogger(OUT_DIR+'Compileresults_filtering')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/compileresults/filtering.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists(OUT_DIR+'Analysis_Results/primary_alignment/filteringbamsStats_%s.html'%args.logfile.rstrip('.log')):
        return passed
    try:
        result = subprocess.run('multiqc %slogs/filtered_bams --force -v -o %sAnalysis_Results/primary_alignment -n filteringbamsStats_%s.html'%(OUT_DIR, OUT_DIR, args.logfile.rstrip('.log')), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/primary_alignment/filteringbamsStats_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 11
def GetBigwigs_BamCoverage():
    if not os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws'):
        os.mkdir(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws')
    if not os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws'):
        os.mkdir(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws')
    if not os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization'):
        os.mkdir(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization')
    if not os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups'):
        os.mkdir(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups')
    if not os.path.exists(OUT_DIR+'logs/RPGC_and_Unnormalized_bws'):
        os.mkdir(OUT_DIR+'logs/RPGC_and_Unnormalized_bws')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger(OUT_DIR+'GetBigwigs_BamCoverage')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/RPGC_and_Unnormalized_bws/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bw'%f) and os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo_norm.bw'%f) and os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups/%s_wo_norm_wDups.bw'%f) and os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_CPM.bw'%f):
            continue
        # RPGC normalized without duplicates
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bw %s --effectiveGenomeSize %s'%(OUT_DIR, f, OUT_DIR, f, args.bamCov_RPGC, EFFECTIVEGENOMESIZE)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # RPGC normalized without duplicated bedgraph for peaks
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam --outFileFormat bedgraph -o %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph %s --effectiveGenomeSize %s'%(OUT_DIR, f, OUT_DIR, f, args.bamCov_RPGC, EFFECTIVEGENOMESIZE)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # No normalization without duplicates
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o %sAnalysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo_norm.bw %s'%(OUT_DIR, f, OUT_DIR, f, args.bamCov_min)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # No normalization without duplicates bedgraph for peaks
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam --outFileFormat bedgraph -o %sAnalysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo_norm.bedgraph %s'%(OUT_DIR, f, OUT_DIR, f, args.bamCov_min)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # No normalization with duplicates
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.bam -o %sAnalysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups/%s_wo_norm_wDups.bw %s'%(OUT_DIR, f, OUT_DIR, f, args.bamCov_min)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # CPM normalized without duplicates
        '''
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_CPM.bw %s --effectiveGenomeSize %s'%(OUT_DIR, f, OUT_DIR, f, args.bamCov_CPM, EFFECTIVEGENOMESIZE)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        '''
        '''
        # Get bedgraph
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam --outFileFormat bedgraph -o %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph %s --effectiveGenomeSize %s'%(OUT_DIR, f, OUT_DIR, f, args.bamCov_RPGC, EFFECTIVEGENOMESIZE)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        '''
        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bw'%f)
        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo_norm.bw'%f)
        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups/%s_wo_norm_wDups.bw'%f)
        #passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_CPM.bw'%f)
    return passed
        
# Step 12
def Map2Spikein_Bowtie2():
    if not os.path.exists(OUT_DIR+'logs/Spike_Alignment'):
        os.mkdir(OUT_DIR+'logs/Spike_Alignment')
    if not os.path.exists(OUT_DIR+'logs/Spike_Alignment/bowtie2'):
        os.mkdir(OUT_DIR+'logs/Spike_Alignment/bowtie2')
    if not os.path.exists(OUT_DIR+'logs/Spike_Alignment/picard_sort'):
        os.mkdir(OUT_DIR+'logs/Spike_Alignment/picard_sort')
    if not os.path.exists(OUT_DIR+'All_output/Spike_mapped_reads'):
        os.mkdir(OUT_DIR+'All_output/Spike_mapped_reads')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger(OUT_DIR+'Map2Spikein_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/Spike_Alignment/bowtie2/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists(OUT_DIR+'All_output/Spike_mapped_reads/%s.bam'%f) and os.path.exists(OUT_DIR+'All_output/Spike_mapped_reads/%s.coordsorted.bam'%f) and os.path.exists(OUT_DIR+'All_output/Spike_mapped_reads/%s.coordsorted.bam.bai'%f):
            continue
        #if args.no_spikein:
        #    logger.info('No spikein, skipping step.')
        #    break
        # Bowtie2
        try:
            result = subprocess.run(('bowtie2 %s -x %s -1 %sAll_output/Trimmed_reads/%s_Trimmed_R1.fastq -2 %sAll_output/Trimmed_reads/%s_Trimmed_R2.fastq 2> %slogs/Spike_Alignment/bowtie2/%s.log | samtools view -Sb - > %sAll_output/Spike_mapped_reads/%s.bam'%(args.spike_align, SPIKEINDEX, OUT_DIR, f, OUT_DIR, f, OUT_DIR, f, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Picard
        logger = logging.getLogger(OUT_DIR+'Map2Spikein_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/Spike_Alignment/picard_sort/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('%s SortSam -I %sAll_output/Spike_mapped_reads/%s.bam -O %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -SORT_ORDER coordinate'%(args.PicardLoc, OUT_DIR, f, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Samtools
        try:
            result = subprocess.run(('samtools index %sAll_output/Spike_mapped_reads/%s.coordsorted.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
            
        passed = passed and os.path.exists(OUT_DIR+'All_output/Spike_mapped_reads/%s.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'All_output/Spike_mapped_reads/%s.coordsorted.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'All_output/Spike_mapped_reads/%s.coordsorted.bam.bai'%f)
    return passed

# Step 13
def Collect_Spikealignment_stats():
    if not os.path.exists(OUT_DIR+'logs/Spike_Alignment/dupstats'):
        os.mkdir(OUT_DIR+'logs/Spike_Alignment/dupstats')
    if not os.path.exists(OUT_DIR+'logs/Spike_Alignment/flagstat'):
        os.mkdir(OUT_DIR+'logs/Spike_Alignment/flagstat')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger(OUT_DIR+'Collect_Spikealignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/Spike_Alignment/dupstats/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists(OUT_DIR+'logs/Spike_Alignment/dupstats/%s.dupMarked.bam'%f) and os.path.exists(OUT_DIR+'logs/Spike_Alignment/dupstats/%s_picard.dupMark.txt'%f):
            continue
        #if args.no_spikein:
        #    logger.info('No spikein, skipping step.')
        #    break
        try:
            result = subprocess.run(('%s MarkDuplicates -I %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -O %slogs/Spike_Alignment/dupstats/%s.dupMarked.bam -METRICS_FILE %slogs/Spike_Alignment/dupstats/%s_picard.dupMark.txt'%(args.PicardLoc, OUT_DIR, f, OUT_DIR, f, OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        logger = logging.getLogger(OUT_DIR+'Collect_Spikealignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/Spike_Alignment/flagstat/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('samtools flagstat %sAll_output/Spike_mapped_reads/%s.coordsorted.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists(OUT_DIR+'logs/Spike_Alignment/dupstats/%s.dupMarked.bam'%f)
        passed = passed and os.path.exists(OUT_DIR+'logs/Spike_Alignment/dupstats/%s_picard.dupMark.txt'%f)
    return passed
    
# Step 14
def Compileresults_Spike():
    if not os.path.exists(OUT_DIR+'Analysis_Results/Spikein_alignment'):
        os.mkdir(OUT_DIR+'Analysis_Results/Spikein_alignment')
    passed = True
    logger = logging.getLogger(OUT_DIR+'Compileresults_Spike')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/compileresults/Spikealign.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists(OUT_DIR+'Analysis_Results/Spikein_alignment/Spike_alignment_%s.html'%args.logfile.rstrip('.log')):
        return passed
    #if args.no_spikein:
    #    logger.info('No spikein, skipping step.')
    #    return passed
    try:
        result = subprocess.run('multiqc %slogs/Spike_Alignment --force -v -o %sAnalysis_Results/Spikein_alignment -n Spike_alignment_%s.html'%(OUT_DIR, OUT_DIR, args.logfile.rstrip('.log')), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Spikein_alignment/Spike_alignment_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 15
def CalcNormFactors():
    if not os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs'):
        os.mkdir(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs')
    if not os.path.exists(OUT_DIR+'Analysis_Results/Spikein_alignment/Spike_alignment_%s_data'%args.logfile.rstrip('.log')):
        os.mkdir(OUT_DIR+'Analysis_Results/Spikein_alignment/Spike_alignment_%s_data'%args.logfile.rstrip('.log'))
    Allsamples = list(fastqfiles)
    passed = True
    logger = logging.getLogger(OUT_DIR+'CalcNormFactors')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/compileresults/Scalefacs.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%args.logfile.rstrip('.log')):
        return passed
    #if args.no_spikein:
    #    logger.info('No spikein, skipping step.')
    #    return passed
    try:
        filetoread = OUT_DIR+"Analysis_Results/Spikein_alignment/Spike_alignment_%s_data/multiqc_bowtie2.txt"%args.logfile.rstrip('.log')
        logger.info('Looking to read %s'%filetoread)
        Spike = pd.read_csv(filetoread, sep='\t')
        Spike.set_index('Sample', inplace=True)
        if len(Spike.index) == len(Allsamples):
            logger.info('mutiqc_bowtie2.txt processing...')
            Spike["TotalMappedFragments"] = Spike["paired_aligned_multi"] + Spike["paired_aligned_one"]
            Spike.columns = [str(col) + '_spikein' for col in Spike.columns]
            Smin = Spike['TotalMappedFragments_spikein'].min()
            Spike['ScalingFactors'] = Smin / Spike['TotalMappedFragments_spikein']
            Spike.to_csv(OUT_DIR+"Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv"%args.logfile.rstrip('.log'))
            logger.info('Saving to %sAnalysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%(OUT_DIR, args.logfile.rstrip('.log')))
        else:
            missingfileAll = (set(Allsamples)).difference(set(Spike.index))
            logger.warning('Sample is/are missing')
            logger.warning(str(missingfileAll))
            passed = False
    except:
        logger.exception("Sorry, multiqc_bowtie2.txt file or data for some files is missing.....")
        passed = False
    passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%args.logfile.rstrip('.log'))
    return passed
    
# Step 16
def GetNormBwsBdgs_BamCoverage():
    if not os.path.exists(OUT_DIR+'logs/Spikein_normalized_bws_bdgs'):
        os.mkdir(OUT_DIR+'logs/Spikein_normalized_bws_bdgs')
    if not os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs'):
        os.mkdir(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs')
    if not os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs'):
        os.mkdir(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs')
    if not os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups'):
        os.mkdir(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups')
    passed = True
    #if args.no_spikein:
    #    print('No spikein, skipping step.')
    #    return passed
    AlignStats = pd.read_csv(OUT_DIR+"Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv"%args.logfile.rstrip('.log'))
    AlignStats.set_index('Sample', inplace=True)
    
    for f in fastqfiles:
        logger = logging.getLogger(OUT_DIR+'GetNormBwsBdgs_BamCoverage')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(OUT_DIR+'logs/Spikein_normalized_bws_bdgs/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/%s_Norm_wDups.bw'%f) and os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/%s_Norm.bw'%f) and os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph'%f):
            continue
        # with dups
        Sfvalue = AlignStats.loc[f, 'ScalingFactors']
        print(Sfvalue)
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.bam -o %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/%s_Norm_wDups.bw --scaleFactor %s %s'%(OUT_DIR, f, OUT_DIR, f, Sfvalue, args.bamCov_default)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # no dups
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/%s_Norm.bw --scaleFactor %s %s'%(OUT_DIR, f, OUT_DIR, f, Sfvalue, args.bamCov_default)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # bedgraph for peak calling
        try:
            result = subprocess.run(('bamCoverage --bam %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam --outFileFormat bedgraph -o %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph --scaleFactor %s %s'%(OUT_DIR, f, OUT_DIR, f, Sfvalue, args.bamCov_default)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
            
        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/%s_Norm_wDups.bw'%f)
        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/%s_Norm.bw'%f)
        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph'%f)
    return passed

# Step 17
def Peak_Calling():
    if not os.path.exists(OUT_DIR+'logs/Peaks'):
        os.mkdir(OUT_DIR+'logs/Peaks')
    if not os.path.exists(OUT_DIR+'Analysis_Results/Peaks'):
        os.mkdir(OUT_DIR+'Analysis_Results/Peaks')
    passed = True
    logger = logging.getLogger(OUT_DIR+'Peak_Calling')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/Peaks/'+args.logfile)
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    
    # Get peaks for spike-in normalized
    if not args.no_spikein:
        # ========== SEACR peak calling spike-in ==========
        # If control reads given
        if len(args.controls) >= 1:
            for c in args.controls:
                logger.info('SEACR - Using control %s'%c)
                for r in fastqfiles:
                    if c != r:
                        if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r):
                            continue
                        try:
                            result = subprocess.run(('bash %s %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph non stringent %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, c, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        try:
                            result = subprocess.run(('bash %s %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph non relaxed %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, c, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r)
        else:
            logger.info('SEACR - No control')
            for r in fastqfiles:
                if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r):
                    continue
                try:
                    result = subprocess.run(('bash %s %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph 0.05 non stringent %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                try:
                    result = subprocess.run(('bash %s %sAnalysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph 0.05 non relaxed %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                
                passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r)
                
        # =========== MACS2 peak calling - spike-in ==========
        # If control reads given
        if len(args.controls) >= 1:
            for c in args.controls:
                logger.info('MACS2 - Using control %s'%c)
                for r in fastqfiles:
                    if c != r:
                        if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r):
                            continue
                        try:
                            result = subprocess.run(('macs2 callpeak -t %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -c %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -f BAM -g %s --outdir %sAnalysis_Results/Peaks/ -n %s'%(OUT_DIR, r, OUT_DIR, c, EFFECTIVEGENOMESIZE, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r)
        else:
            logger.info('MACS2 - No control')
            for r in fastqfiles:
                if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r):
                    continue
                try:
                    result = subprocess.run(('macs2 callpeak -t %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -f BAM -g %s --outdir %sAnalysis_Results/Peaks/ -n %s'%(OUT_DIR, r, EFFECTIVEGENOMESIZE, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r)
        
        # ========== GoPeaks peak calling - spike-in ==========
        # If control reads given
        if len(args.controls) >= 1:
            for c in args.controls:
                logger.info('GoPeaks - Using control %s'%c)
                for r in fastqfiles:
                    if c != r:
                        if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r):
                            continue
                        try:
                            result = subprocess.run(('./gopeaks -b %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -c %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -o %sAnalysis_Results/Peaks/%s_gopeaks'%(OUT_DIR, r, OUT_DIR, c, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r)
        else:
            logger.info('GoPeaks - No control')
            for r in fastqfiles:
                if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r):
                    continue
                try:
                    result = subprocess.run(('./gopeaks -b %sAll_output/Spike_mapped_reads/%s.coordsorted.bam -o %sAnalysis_Results/Peaks/%s_gopeaks'%(OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r)
        
    else:
        # Get peaks for NON-spikein-normalized bw/bedgraphs
        # ========== SEACR peak calling - NO spike-in ==========
        # If control reads given
        if len(args.controls) >= 1:
            for c in args.controls:
                logger.info('SEACR - Using control %s'%c)
                for r in fastqfiles:
                    if c != r:
                        if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r):
                            continue
                        try:
                            result = subprocess.run(('bash %s %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph non stringent %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, c, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        try:
                            result = subprocess.run(('bash %s %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph non relaxed %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, c, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r)
        else:
            logger.info('SEACR - No control')
            for r in fastqfiles:
                if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r):
                    continue
                # RPGC-norm bedgraph peaks
                try:
                    result = subprocess.run(('bash %s %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph 0.05 non stringent %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                try:
                    result = subprocess.run(('bash %s %sAnalysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph 0.05 non relaxed %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                # Non-norm, no duplicates, bedgraph peaks
                try:
                    result = subprocess.run(('bash %s %sAnalysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo_norm.bedgraph 0.05 non stringent %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                try:
                    result = subprocess.run(('bash %s %sAnalysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo_norm.bedgraph 0.05 non relaxed %sAnalysis_Results/Peaks/%s'%(args.SEACRLoc, OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s.stringent.bed'%r)
        
        # ========== MACS2 peak calling - NO spike-in ==========
        # If control reads given
        if len(args.controls) >= 1:
            for c in args.controls:
                logger.info('MACS2 - Using control %s'%c)
                for r in fastqfiles:
                    if c != r:
                        if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r):
                            continue
                        try:
                            result = subprocess.run(('macs2 callpeak -t %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -c %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -f BAM -g %s --outdir %sAnalysis_Results/Peaks/ -n %s'%(OUT_DIR, r, OUT_DIR, c, EFFECTIVEGENOMESIZE, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r)
        else:
            logger.info('MACS2 - No control')
            for r in fastqfiles:
                if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r):
                    continue
                try:
                    result = subprocess.run(('macs2 callpeak -t %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -f BAM -g %s --outdir %sAnalysis_Results/Peaks/ -n %s'%(OUT_DIR, r, EFFECTIVEGENOMESIZE, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_peaks.narrowPeak'%r)
        
        # ========== GoPeaks peak calling - NO spike-in ==========
        # If control reads given
        if len(args.controls) >= 1:
            for c in args.controls:
                logger.info('GoPeaks - Using control %s'%c)
                for r in fastqfiles:
                    if c != r:
                        if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r):
                            continue
                        try:
                            result = subprocess.run(('./gopeaks -b %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -c %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o %sAnalysis_Results/Peaks/%s_gopeaks'%(OUT_DIR, r, OUT_DIR, c, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                            logger.info(result.stdout.rstrip('\n'))
                            logger.warning(result.stderr.rstrip('\n'))
                        except Exception as e:
                            logger.exception(e)
                            passed = False
                        passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r)
        else:
            logger.info('GoPeaks - No control')
            for r in fastqfiles:
                if os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r):
                    continue
                try:
                    result = subprocess.run(('./gopeaks -b %sAll_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o %sAnalysis_Results/Peaks/%s_gopeaks'%(OUT_DIR, r, OUT_DIR, r)), shell=True, capture_output=True, text=True)
                    logger.info(result.stdout.rstrip('\n'))
                    logger.warning(result.stderr.rstrip('\n'))
                except Exception as e:
                    logger.exception(e)
                    passed = False
                passed = passed and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_gopeaks.json'%r) and os.path.exists(OUT_DIR+'Analysis_Results/Peaks/%s_gopeaks_peaks.bed'%r)
        
    return passed

# Step 18
def Clean_up():
    logger = logging.getLogger(OUT_DIR+'Clean_up')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/cleanup.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    passed = True
    if args.cleanup:
        logger.info("Removing files except for bigwigs/bigbeds")
        try:
            result = subprocess.run(('rm %sAnalysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%(OUT_DIR, args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm %sAnalysis_Results/Spikein_alignment/Spike_alignment_%s.html'%(OUT_DIR, args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm %sAnalysis_Results/primary_alignment/filteringbamsStats_%s.html'%(OUT_DIR, args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm %sAnalysis_Results/primary_alignment/Alignment_results_%s.html'%(OUT_DIR, args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm %sAnalysis_Results/Trimming/PostTrimming_QC_%s.html'%(OUT_DIR, args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm %sAnalysis_Results/QC_Rawreads/Rawreads_QC_%s.html'%(OUT_DIR, args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm %s*.out'), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        for f in fastqfiles:
            try:
                result = subprocess.run(('rm %slogs/Spike_Alignment/dupstats/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %sAll_output/Spike_mapped_reads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %sAll_output/Processed_reads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %slogs/filtered_bams/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %slogs/primary_alignment/PostAlignmentStats/dupstats/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %sAll_output/Mapped_reads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %sAll_output/Trimmed_reads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %sAnalysis_Results/QC_Rawreads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
            except Exception as e:
                logger.exception(e)
                passed = False
    else:
        logger.info("Keeping all generated files.\nCheck disk space and manually remove large reads, .bam, and .bai files if desired.")
        '''
        logger.info("Removing intermediate large reads, .bam, and .bai files")
        for f in fastqfiles:
            try:
                result = subprocess.run(('rm %sAll_output/Processed_reads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %sAll_output/Mapped_reads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %sAll_output/Trimmed_reads/%s*'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
                result = subprocess.run(('rm %slogs/primary_alignment/PostAlignmentStats/dupstats/%s*.bam'%(OUT_DIR, f)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
            except Exception as e:
                logger.exception(e)
                passed = False
        '''
    return passed

if __name__ == '__main__':
    t_begin = time.time()
    if not os.path.exists(OUT_DIR+'logs'):
        os.mkdir(OUT_DIR+'logs')
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(OUT_DIR+'logs/%s'%args.logfile)
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info('Prepping multiqc_config.yaml...')
    
    # Multiqc configuration
    home = str(Path.home())
    datadict = {'log_filesize_limit': 2000000000}
    files = glob.glob(os.path.join(home, ".*.yaml"))
    if os.path.join(home, ".multiqc_config.yaml") in files:
        logger.info('~/.multiqc_config.yaml is present')
        with open(os.path.join(home, ".multiqc_config.yaml"), 'r') as file:
            values = yaml.safe_load(file)
            if 'log_filesize_limit' in str(values):
                logger.info('log_filesize_limit is already set')
            if values is None:
                with open(os.path.join(home, ".multiqc_config.yaml"), 'w') as file:
                    docs = yaml.dump(datadict, file)
                    logger.info('added log_filesize_limit to multiqc log file')
            else:
                with open(os.path.join(home, ".multiqc_config.yaml"), 'r') as file:
                    new_yaml = yaml.safe_load(file)
                    logger.info(type(new_yaml))
                    new_yaml.update(datadict)
                with open(os.path.join(home, ".multiqc_config.yaml"),'w') as file:
                    yaml.safe_dump(new_yaml, file)
                    logger.info('log_filesize_limit: 2000000000 is set')
    else:
    	with open(os.path.join(home, ".multiqc_config.yaml"), 'w') as file:
                documents = yaml.dump(datadict, file)
                logger.info('made new .multiqc_config.yaml')    
    
    logger.info("Getting pipeline steps...")
    pipeline = [ f for f in globals().values() if type(f) == types.FunctionType ]
    if args.no_spikein:
        del pipeline[11:-2]
    #if not args.cleanup:
    #    del pipeline[-1]
        
    logger.info("Running pipeline...")
    for p in range(0, len(pipeline)):
        logger.info("Step %s/%s - %s"%(p+1, len(pipeline), pipeline[p].__name__))
        t_start = time.time()
        if pipeline[p]():
            logger.info("PASSED - duration: %s minutes"%(round((time.time() - t_start)/60, 4)))
            continue
        else:
            logger.error("FAILED - duration: %s minutes"%(round((time.time() - t_start)/60, 4)))
            exit()
    logger.info("Pipeline COMPLETE!")
    logger.info('Duration: %s minutes'%(round((time.time() - t_begin)/60, 4)))
    