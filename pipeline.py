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
    'Peaks_SEACR',
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
parser.add_argument('-reads', '--reads', help='Path to directory containing reads (R1, R2, and md5sum.txt files)', type=str)
parser.add_argument('-species', '--species', help='Species of reads (Mus for mouse, Homo for human)', type=str, choices=['Mus', 'Homo'], default='Mus')
parser.add_argument('-spikein', '--spikein', help='Spikein type', type=str, choices=['Amp', 'Bacteria'], default='Amp')
parser.add_argument('-length', '--length', help='Read length', type=str, default='50', choices=['50', '75', '100', '150', '200'])
parser.add_argument('-controls', '--controls', help='Control reads for peaks calling', default=[], nargs='+')
parser.add_argument('-no_spikein', '--no_spikein', help='If no spikein, skip steps for normalizing to spikein', action='store_true')
parser.add_argument('-cleanup', '--cleanup', help='If cleanup, remove all intermediate files keeping only final .bw and .bed files', action='store_true')
# Program and reference genome locations
parser.add_argument('-PicardLoc', '--PicardLoc', help='Location of picard.jar', type=str, default="java -jar /home/earezza/projects/def-jdilwort/earezza/picard.jar")
parser.add_argument('-SEACRLoc', '--SEACRLoc', help='Location of SEACR .sh', type=str, default="/home/earezza/projects/def-jdilwort/earezza/SEACR/SEACR_1.3.sh")
parser.add_argument('-bowtie2_index', '--bowtie2_index', help='Location of bowtie2 genome index files (from iGenomes)', type=str, default="/home/earezza/projects/def-jdilwort/earezza/CnT_pipeline_snakemake/Reference_files/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome")
# Usually unchanged command line options for programs
parser.add_argument('-spike_align', '--spike_align', help='Command input options for spikein alignment', type=str, default="-p 8 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700")
parser.add_argument('-bamCov_default', '--bamCov_default', help='Command input default for bamCoverage', type=str, default="--binSize 10 --ignoreForNormalization 'chrM' --extendReads")
parser.add_argument('-bamCov_min', '--bamCov_min', help='Command input options for bamCoverage', type=str, default="--binSize 10 --extendReads")
parser.add_argument('-bamCov_RPGC', '--bamCov_RPGC', help='Command input options for bamCoverage, reads per genomic content (1x normalization)', type=str, default="--binSize 10 --normalizeUsing 'RPGC'  --ignoreForNormalization 'chrM' --extendReads")
parser.add_argument('-bamCov_CPM', '--bamCov_CPM', help='Command input options for bamCoverage, number of reads per bin / number of mapped reads (in millions)', type=str, default="--binSize 10 --normalizeUsing 'CPM'  --ignoreForNormalization 'chrM' --extendReads")
parser.add_argument('-genome_align', '--genome_align', help='Command input options for genome alignment', type=str, default="-p 8 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
parser.add_argument('-samtools_mapq', '--samtools_mapq', help='Command input option for samtools mapq', type=str, default="-q 10")
parser.add_argument('-samtools_proper_paired', '--samtools_proper_paired', help='Command input option for samtools', type=str, default="-f 2")
# Spikein reference index files
parser.add_argument('-spikein_index_amp', '--spikein_index_amp', help="Source of reference files for Amp spikein index file", type=str, default='/home/earezza/projects/def-jdilwort/earezza/CnT_pipeline_snakemake/Reference_files/Spikein_indices/Amp_pbluescript/Amp_index/Amp_pBlue')
parser.add_argument('-spikein_index_Ecoli', '--spikein_index_Ecoli', help="Source of reference files for Amp spikein index file", type=str, default='/home/earezza/projects/def-jdilwort/earezza/CnT_pipeline_snakemake/Reference_files/Spikein_indices/EcoliK12_index/EcoliK12Index/EcoliK12')
args = parser.parse_args()

# Define constants
ROOT_DIR = os.getcwd()
READS_DIR = args.reads


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
    file_handler = logging.FileHandler('logs/md5checks.log')
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
    if not os.path.exists('Analysis_Results'):
        os.mkdir('Analysis_Results')
    if not os.path.exists('Analysis_Results/QC_Rawreads'):
        os.mkdir('Analysis_Results/QC_Rawreads')
    if not os.path.exists('logs/fastqc_rawreads'):
        os.mkdir('logs/fastqc_rawreads')
    passed = True
    for f in fastqfiles:
        for r in reads:
            logger = logging.getLogger('QCrawreads_Fastqc')
            logger.setLevel(logging.DEBUG)
            file_handler = logging.FileHandler('logs/fastqc_rawreads/%s_%s.log'%(f, r))
            formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            if os.path.exists('Analysis_Results/QC_Rawreads/%s_%s_fastqc.html'%(f,r)) and os.path.exists('Analysis_Results/QC_Rawreads/%s_%s_fastqc.zip'%(f,r)):
                continue
            try:
                result = subprocess.run(('fastqc %s/%s_%s.fastq.gz --outdir=./Analysis_Results/QC_Rawreads'%(READS_DIR, f, r)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
            except Exception as e:
                logger.exception(e)
                passed = False
            passed = passed and os.path.exists('Analysis_Results/QC_Rawreads/%s_%s_fastqc.html'%(f,r))
            passed = passed and os.path.exists('Analysis_Results/QC_Rawreads/%s_%s_fastqc.zip'%(f,r))
    return passed

# Step 3
def Compileresults_QC():
    if not os.path.exists('logs/compileresults'):
        os.mkdir('logs/compileresults')
    passed = True
    logger = logging.getLogger('Compileresults_QC')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/compileresults/QC.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists('Analysis_Results/QC_Rawreads/Rawreads_QC_%s.html'%(args.logfile.rstrip('.log'))):
        return passed
    try:
        result = subprocess.run('multiqc ./Analysis_Results/QC_Rawreads --force -v -o ./Analysis_Results/QC_Rawreads -n Rawreads_QC_%s.html'%args.logfile.rstrip('.log'), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists('Analysis_Results/QC_Rawreads/Rawreads_QC_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 4
def AdapterTrim_Cutadapt():
    if not os.path.exists('All_output'):
        os.mkdir('All_output')
    if not os.path.exists('All_output/Trimmed_reads'):
        os.mkdir('All_output/Trimmed_reads')
    if not os.path.exists('logs/cutadapt'):
        os.mkdir('logs/cutadapt')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger('AdapterTrim_Cutadapt')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/cutadapt/%s.log'%(f))
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists('All_output/Trimmed_reads/%s_Trimmed_R1.fastq'%(f)) and os.path.exists('All_output/Trimmed_reads/%s_Trimmed_R2.fastq'%(f)):
            continue
        try:
            result = subprocess.run(('cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT  -o All_output/Trimmed_reads/%s_Trimmed_R1.fastq -p All_output/Trimmed_reads/%s_Trimmed_R2.fastq %s_R1.fastq.gz %s_R2.fastq.gz'%(f, f, READS_DIR+'/'+f, READS_DIR+'/'+f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists('All_output/Trimmed_reads/%s_Trimmed_R1.fastq'%(f))
        passed = passed and os.path.exists('All_output/Trimmed_reads/%s_Trimmed_R2.fastq'%(f))
    return passed

# Step 5
def Compileresults_PosttrimQC():
    if not os.path.exists('Analysis_Results/Trimming'):
        os.mkdir('Analysis_Results/Trimming')
    passed = True
    logger = logging.getLogger('Compileresults_PosttrimQC')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/compileresults/PosttrimQC.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists('Analysis_Results/Trimming/PostTrimming_QC_%s.html'%args.logfile.rstrip('.log')):
        return passed
    try:
        result = subprocess.run('multiqc ./logs/cutadapt --force -v -o ./Analysis_Results/Trimming -n PostTrimming_QC_%s.html'%args.logfile.rstrip('.log'), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists('Analysis_Results/Trimming/PostTrimming_QC_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 6
def Map_Bowtie2():
    if not os.path.exists('All_output/Mapped_reads'):
        os.mkdir('All_output/Mapped_reads')
    if not os.path.exists('logs/primary_alignment'):
        os.mkdir('logs/primary_alignment')
    if not os.path.exists('logs/primary_alignment/bowtie2'):
        os.mkdir('logs/primary_alignment/bowtie2')
    if not os.path.exists('logs/primary_alignment/picard_sort'):
        os.mkdir('logs/primary_alignment/picard_sort')
    passed = True
    for f in fastqfiles:
        # Bowtie2
        logger = logging.getLogger('Map_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/primary_alignment/bowtie2/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists('All_output/Mapped_reads/%s.bam'%f) and os.path.exists('All_output/Mapped_reads/%s.coordsorted.bam'%f) and os.path.exists('All_output/Mapped_reads/%s.coordsorted.bam.bai'%f):
            continue
        try:
            result = subprocess.run(('bowtie2 %s -x %s -1 All_output/Trimmed_reads/%s_Trimmed_R1.fastq -2 All_output/Trimmed_reads/%s_Trimmed_R2.fastq 2> logs/primary_alignment/bowtie2/%s.log | samtools view -bS - > All_output/Mapped_reads/%s.bam'%(args.genome_align, args.bowtie2_index, f, f, f, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.error(e)
            passed = False
        # Picard
        logger = logging.getLogger('Map_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/primary_alignment/picard_sort/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('%s SortSam -I All_output/Mapped_reads/%s.bam -O All_output/Mapped_reads/%s.coordsorted.bam -SORT_ORDER coordinate'%(args.PicardLoc, f, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Samtools
        try:
            result = subprocess.run(('samtools index All_output/Mapped_reads/%s.coordsorted.bam'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists('All_output/Mapped_reads/%s.bam'%f)
        passed = passed and os.path.exists('All_output/Mapped_reads/%s.coordsorted.bam'%f)
        passed = passed and os.path.exists('All_output/Mapped_reads/%s.coordsorted.bam.bai'%f)
    return passed    
        
# Step 7
def Collect_alignment_stats():
    if not os.path.exists('logs/primary_alignment/PostAlignmentStats'):
        os.mkdir('logs/primary_alignment/PostAlignmentStats')
    if not os.path.exists('logs/primary_alignment/PostAlignmentStats/dupstats'):
        os.mkdir('logs/primary_alignment/PostAlignmentStats/dupstats')
    if not os.path.exists('logs/primary_alignment/PostAlignmentStats/flagstat'):
        os.mkdir('logs/primary_alignment/PostAlignmentStats/flagstat')
    passed = True
    for f in fastqfiles:
        # Picard
        logger = logging.getLogger('Collect_alignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/primary_alignment/PostAlignmentStats/dupstats/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists('logs/primary_alignment/PostAlignmentStats/dupstats/%s.dupMarked.bam'%f) and os.path.exists('logs/primary_alignment/PostAlignmentStats/dupstats/%s_picard.dupMark.txt'%f):
            continue
        try:
            result = subprocess.run(('%s MarkDuplicates -I All_output/Mapped_reads/%s.coordsorted.bam -O logs/primary_alignment/PostAlignmentStats/dupstats/%s.dupMarked.bam -METRICS_FILE logs/primary_alignment/PostAlignmentStats/dupstats/%s_picard.dupMark.txt'%(args.PicardLoc, f, f, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Samtools
        logger = logging.getLogger('Collect_alignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/primary_alignment/PostAlignmentStats/flagstat/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('samtools flagstat All_output/Mapped_reads/%s.coordsorted.bam'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists('logs/primary_alignment/PostAlignmentStats/dupstats/%s.dupMarked.bam'%f)
        passed = passed and os.path.exists('logs/primary_alignment/PostAlignmentStats/dupstats/%s_picard.dupMark.txt'%f)
    return passed

# Step 8
def Compileresults_map():
    passed = True
    logger = logging.getLogger('Compileresults_map')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/compileresults/map.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists('Analysis_Results/primary_alignment/Alignment_results_%s.html'%args.logfile.rstrip('.log')):
        return passed
    try:
        result = subprocess.run('multiqc ./logs/primary_alignment --force -v -o ./Analysis_Results/primary_alignment -n Alignment_results_%s.html'%args.logfile.rstrip('.log'), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists('Analysis_Results/primary_alignment/Alignment_results_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 9
def Filtering_bams_PicardSamtools():
    if not os.path.exists('logs/filtered_bams'):
        os.mkdir('logs/filtered_bams')
    if not os.path.exists('All_output/Processed_reads'):
        os.mkdir('All_output/Processed_reads')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger('Filtering_bams_PicardSamtools')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/filtered_bams/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if ('All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%f) and os.path.exists('All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam.bai'%f) and os.path.exists('All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%f) and os.path.exists('logs/filtered_bams/%s_picard.rmDup.txt'%f) and os.path.exists('All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam.bai'%f):
            continue
        try:
            result = subprocess.run(('samtools view -bu %s All_output/Mapped_reads/%s.coordsorted.bam | samtools view -b %s - | samtools sort - -o All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%(args.samtools_proper_paired, f, args.samtools_mapq, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('%s MarkDuplicates -I All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam -O All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -REMOVE_DUPLICATES true -METRICS_FILE logs/filtered_bams/%s_picard.rmDup.txt'%(args.PicardLoc, f, f, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('samtools flagstat All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('samtools index All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('samtools index All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
    
        passed = passed and os.path.exists('All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam'%f)
        passed = passed and os.path.exists('All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam.bai'%f)
        passed = passed and os.path.exists('All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam'%f)
        passed = passed and os.path.exists('logs/filtered_bams/%s_picard.rmDup.txt'%f)
        passed = passed and os.path.exists('All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam.bai'%f)
    return passed

# Step 10
def Compileresults_filtering():
    passed = True
    logger = logging.getLogger('Compileresults_filtering')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/compileresults/filtering.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists('Analysis_Results/primary_alignment/filteringbamsStats_%s.html'%args.logfile.rstrip('.log')):
        return passed
    try:
        result = subprocess.run('multiqc ./logs/filtered_bams --force -v -o ./Analysis_Results/primary_alignment -n filteringbamsStats_%s.html'%args.logfile.rstrip('.log'), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists('Analysis_Results/primary_alignment/filteringbamsStats_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 11
def GetBigwigs_BamCoverage():
    if not os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws'):
        os.mkdir('Analysis_Results/RPGC_and_Unnormalized_bws')
    if not os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws'):
        os.mkdir('Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws')
    if not os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization'):
        os.mkdir('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization')
    if not os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups'):
        os.mkdir('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups')
    if not os.path.exists('logs/RPGC_and_Unnormalized_bws'):
        os.mkdir('logs/RPGC_and_Unnormalized_bws')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger('GetBigwigs_BamCoverage')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/RPGC_and_Unnormalized_bws/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bw'%f) and os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo.norm.bw'%f) and os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups/%s_wo.norm_wDups.bw'%f) and os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_CPM.bw'%f):
            continue
        # RPGC normalized without duplicates
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bw %s --effectiveGenomeSize %s'%(f, f, args.bamCov_RPGC, EFFECTIVEGENOMESIZE)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # No normalization without duplicates
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo.norm.bw %s'%(f, f, args.bamCov_min)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # No normalization with duplicates
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam -o Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups/%s_wo.norm_wDups.bw %s'%(f, f, args.bamCov_min)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # CPM normalized without duplicates
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_CPM.bw %s --effectiveGenomeSize %s'%(f, f, args.bamCov_CPM, EFFECTIVEGENOMESIZE)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        '''
        # Get bedgraph
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam --outFileFormat bedgraph -o Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bedgraph %s --effectiveGenomeSize %s'%(f, f, args.bamCov_RPGC, EFFECTIVEGENOMESIZE)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        '''
        passed = passed and os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_RPGC.bw'%f)
        passed = passed and os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_normalization/%s_wo.norm.bw'%f)
        passed = passed and os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/bws_wo_norm_wDups/%s_wo.norm_wDups.bw'%f)
        passed = passed and os.path.exists('Analysis_Results/RPGC_and_Unnormalized_bws/RPGC_normalized_bws/%s_CPM.bw'%f)
    return passed
        
# Step 12
def Map2Spikein_Bowtie2():
    if not os.path.exists('logs/Spike_Alignment'):
        os.mkdir('logs/Spike_Alignment')
    if not os.path.exists('logs/Spike_Alignment/bowtie2'):
        os.mkdir('logs/Spike_Alignment/bowtie2')
    if not os.path.exists('logs/Spike_Alignment/picard_sort'):
        os.mkdir('logs/Spike_Alignment/picard_sort')
    if not os.path.exists('All_output/Spike_mapped_reads'):
        os.mkdir('All_output/Spike_mapped_reads')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger('Map2Spikein_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/Spike_Alignment/bowtie2/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists('All_output/Spike_mapped_reads/%s.bam'%f) and os.path.exists('All_output/Spike_mapped_reads/%s.coordsorted.bam'%f) and os.path.exists('All_output/Spike_mapped_reads/%s.coordsorted.bam.bai'%f):
            continue
        #if args.no_spikein:
        #    logger.info('No spikein, skipping step.')
        #    break
        # Bowtie2
        try:
            result = subprocess.run(('bowtie2 %s -x %s -1 All_output/Trimmed_reads/%s_Trimmed_R1.fastq -2 All_output/Trimmed_reads/%s_Trimmed_R2.fastq 2> logs/Spike_Alignment/bowtie2/%s.log | samtools view -Sb - > All_output/Spike_mapped_reads/%s.bam'%(args.spike_align, SPIKEINDEX, f, f, f, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Picard
        logger = logging.getLogger('Map2Spikein_Bowtie2')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/Spike_Alignment/picard_sort/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('%s SortSam -I All_output/Spike_mapped_reads/%s.bam -O All_output/Spike_mapped_reads/%s.coordsorted.bam -SORT_ORDER coordinate'%(args.PicardLoc, f, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # Samtools
        try:
            result = subprocess.run(('samtools index All_output/Spike_mapped_reads/%s.coordsorted.bam'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
            
        passed = passed and os.path.exists('All_output/Spike_mapped_reads/%s.bam'%f)
        passed = passed and os.path.exists('All_output/Spike_mapped_reads/%s.coordsorted.bam'%f)
        passed = passed and os.path.exists('All_output/Spike_mapped_reads/%s.coordsorted.bam.bai'%f)
    return passed

# Step 13
def Collect_Spikealignment_stats():
    if not os.path.exists('logs/Spike_Alignment/dupstats'):
        os.mkdir('logs/Spike_Alignment/dupstats')
    if not os.path.exists('logs/Spike_Alignment/flagstat'):
        os.mkdir('logs/Spike_Alignment/flagstat')
    passed = True
    for f in fastqfiles:
        logger = logging.getLogger('Collect_Spikealignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/Spike_Alignment/dupstats/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists('logs/Spike_Alignment/dupstats/%s.dupMarked.bam'%f) and os.path.exists('logs/Spike_Alignment/dupstats/%s_picard.dupMark.txt'%f):
            continue
        #if args.no_spikein:
        #    logger.info('No spikein, skipping step.')
        #    break
        try:
            result = subprocess.run(('%s MarkDuplicates -I All_output/Spike_mapped_reads/%s.coordsorted.bam -O logs/Spike_Alignment/dupstats/%s.dupMarked.bam -METRICS_FILE logs/Spike_Alignment/dupstats/%s_picard.dupMark.txt'%(args.PicardLoc, f, f, f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        logger = logging.getLogger('Collect_Spikealignment_stats')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/Spike_Alignment/flagstat/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        try:
            result = subprocess.run(('samtools flagstat All_output/Spike_mapped_reads/%s.coordsorted.bam'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        passed = passed and os.path.exists('logs/Spike_Alignment/dupstats/%s.dupMarked.bam'%f)
        passed = passed and os.path.exists('logs/Spike_Alignment/dupstats/%s_picard.dupMark.txt'%f)
    return passed
    
# Step 14
def Compileresults_Spike():
    if not os.path.exists('Analysis_Results/Spikein_alignment'):
        os.mkdir('Analysis_Results/Spikein_alignment')
    passed = True
    logger = logging.getLogger('Compileresults_Spike')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/compileresults/Spikealign.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists('Analysis_Results/Spikein_alignment/Spike_alignment_%s.html'%args.logfile.rstrip('.log')):
        return passed
    #if args.no_spikein:
    #    logger.info('No spikein, skipping step.')
    #    return passed
    try:
        result = subprocess.run('multiqc ./logs/Spike_Alignment --force -v -o ./Analysis_Results/Spikein_alignment -n Spike_alignment_%s.html'%args.logfile.rstrip('.log'), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    passed = passed and os.path.exists('Analysis_Results/Spikein_alignment/Spike_alignment_%s.html'%args.logfile.rstrip('.log'))
    return passed

# Step 15
def CalcNormFactors():
    if not os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs'):
        os.mkdir('Analysis_Results/Spikein_normalized_bws_bdgs')
    if not os.path.exists('Analysis_Results/Spikein_alignment/Spike_alignment_%s_data'%args.logfile.rstrip('.log')):
        os.mkdir('Analysis_Results/Spikein_alignment/Spike_alignment_%s_data'%args.logfile.rstrip('.log'))
    Allsamples = list(fastqfiles)
    passed = True
    logger = logging.getLogger('CalcNormFactors')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/compileresults/Scalefacs.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%args.logfile.rstrip('.log')):
        return passed
    #if args.no_spikein:
    #    logger.info('No spikein, skipping step.')
    #    return passed
    try:
        filetoread = "Analysis_Results/Spikein_alignment/Spike_alignment_%s_data/multiqc_bowtie2.txt"%args.logfile.rstrip('.log')
        logger.info('Looking to read %s'%filetoread)
        Spike = pd.read_csv(filetoread, sep='\t')
        Spike.set_index('Sample', inplace=True)
        if len(Spike.index) == len(Allsamples):
            logger.info('mutiqc_bowtie2.txt processing...')
            Spike["TotalMappedFragments"] = Spike["paired_aligned_multi"] + Spike["paired_aligned_one"]
            Spike.columns = [str(col) + '_spikein' for col in Spike.columns]
            Smin = Spike['TotalMappedFragments_spikein'].min()
            Spike['ScalingFactors'] = Smin / Spike['TotalMappedFragments_spikein']
            Spike.to_csv("Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv"%args.logfile.rstrip('.log'))
            logger.info('Saving to Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%args.logfile.rstrip('.log'))
        else:
            missingfileAll = (set(Allsamples)).difference(set(Spike.index))
            logger.warning('Sample is/are missing')
            passed = False
    except:
        logger.exception("Sorry, multiqc_bowtie2.txt file or data for some files is missing.....")
        passed = False
    passed = passed and os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%args.logfile.rstrip('.log'))
    return passed
    
# Step 16
def GetNormBwsBdgs_BamCoverage():
    if not os.path.exists('logs/Spikein_normalized_bws_bdgs'):
        os.mkdir('logs/Spikein_normalized_bws_bdgs')
    if not os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs'):
        os.mkdir('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs')
    if not os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs'):
        os.mkdir('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs')
    if not os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups'):
        os.mkdir('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups')
    passed = True
    #if args.no_spikein:
    #    print('No spikein, skipping step.')
    #    return passed
    AlignStats = pd.read_csv("Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv"%args.logfile.rstrip('.log'))
    AlignStats.set_index('Sample', inplace=True)
    
    for f in fastqfiles:
        logger = logging.getLogger('GetNormBwsBdgs_BamCoverage')
        logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler('logs/Spikein_normalized_bws_bdgs/%s.log'%f)
        formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        if os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/%s_Norm_wDups.bw'%f) and os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/%s_Norm.bw'%f) and os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph'%f):
            continue
        # with dups
        Sfvalue = AlignStats.loc[f, 'ScalingFactors']
        print(Sfvalue)
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.bam -o Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/%s_Norm_wDups.bw --scaleFactor %s %s'%(f, f, Sfvalue, args.bamCov_default)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        # no dups
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam -o Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/%s_Norm.bw --scaleFactor %s %s'%(f, f, Sfvalue, args.bamCov_default)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
        try:
            result = subprocess.run(('bamCoverage --bam All_output/Processed_reads/%s.MappedPaired.MAPQ10.NoDups.bam --outFileFormat bedgraph -o Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph --scaleFactor %s %s'%(f, f, Sfvalue, args.bamCov_default)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
            
        passed = passed and os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs_wDups/%s_Norm_wDups.bw'%f)
        passed = passed and os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bigwigs/%s_Norm.bw'%f)
        passed = passed and os.path.exists('Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph'%f)
    return passed

# Step 17
def Peaks_SEACR():
    if not os.path.exists('logs/Peaks'):
        os.mkdir('logs/Peaks')
    if not os.path.exists('Analysis_Results/Peaks'):
        os.mkdir('Analysis_Results/Peaks')
    passed = True
    logger = logging.getLogger('Peaks_SEACR')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/Peaks/'+args.logfile)
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    #if args.no_spikein:
    #    logger.info('No spikein, skipping step.')
    #    return passed
    
    # If conrol reads given
    if len(args.controls) >= 1:
        for c in args.controls:
            logger.info('Using control %s'%c)
            for r in fastqfiles:
                if c != r:
                    if os.path.exists('Analysis_Results/Peaks/%s.stringent.bed'%r):
                        continue
                    try:
                        result = subprocess.run(('bash %s Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph non stringent Analysis_Results/Peaks/%s'%(args.SEACRLoc, r, c, r)), shell=True, capture_output=True, text=True)
                        logger.info(result.stdout.rstrip('\n'))
                        logger.warning(result.stderr.rstrip('\n'))
                    except Exception as e:
                        logger.exception(e)
                        passed = False
                    try:
                        result = subprocess.run(('bash %s Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph non relaxed Analysis_Results/Peaks/%s'%(args.SEACRLoc, r, c, r)), shell=True, capture_output=True, text=True)
                        logger.info(result.stdout.rstrip('\n'))
                        logger.warning(result.stderr.rstrip('\n'))
                    except Exception as e:
                        logger.exception(e)
                        passed = False
                    passed = passed and os.path.exists('Analysis_Results/Peaks/%s.stringent.bed'%r)
    else:
        logger.info('No control')
        for r in fastqfiles:
            if os.path.exists('Analysis_Results/Peaks/%s.stringent.bed'%r):
                continue
            try:
                result = subprocess.run(('bash %s Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph 0.05 non stringent Analysis_Results/Peaks/%s'%(args.SEACRLoc, r, r)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
            except Exception as e:
                logger.exception(e)
                passed = False
            try:
                result = subprocess.run(('bash %s Analysis_Results/Spikein_normalized_bws_bdgs/Normalized_bedgraphs/%s_Norm.bedgraph 0.05 non relaxed Analysis_Results/Peaks/%s'%(args.SEACRLoc, r, r)), shell=True, capture_output=True, text=True)
                logger.info(result.stdout.rstrip('\n'))
                logger.warning(result.stderr.rstrip('\n'))
            except Exception as e:
                logger.exception(e)
                passed = False
            passed = passed and os.path.exists('Analysis_Results/Peaks/%s.stringent.bed'%r)
    return passed

# Step 18
def Clean_up():
    logger = logging.getLogger('Clean_up')
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/cleanup.log')
    formatter = logging.Formatter('%(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info("Removing files except for bigwigs/bigbeds")
    try:
        result = subprocess.run(('rm Analysis_Results/Spikein_normalized_bws_bdgs/Spike_align_stats_%s.csv'%(args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
        result = subprocess.run(('rm Analysis_Results/Spikein_alignment/Spike_alignment_%s.html'%(args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
        result = subprocess.run(('rm Analysis_Results/primary_alignment/filteringbamsStats_%s.html'%(args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
        result = subprocess.run(('rm Analysis_Results/primary_alignment/Alignment_results_%s.html'%(args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
        result = subprocess.run(('rm Analysis_Results/Trimming/PostTrimming_QC_%s.html'%(args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
        result = subprocess.run(('rm Analysis_Results/QC_Rawreads/Rawreads_QC_%s.html'%(args.logfile.rstrip('.log'))), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
        result = subprocess.run(('rm *.out'), shell=True, capture_output=True, text=True)
        logger.info(result.stdout.rstrip('\n'))
        logger.warning(result.stderr.rstrip('\n'))
    except Exception as e:
        logger.exception(e)
        passed = False
    for f in fastqfiles:
        try:
            result = subprocess.run(('rm logs/Spike_Alignment/dupstats/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm All_output/Spike_mapped_reads/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm All_output/Processed_reads/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm logs/filtered_bams/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm logs/primary_alignment/PostAlignmentStats/dupstats/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm All_output/Mapped_reads/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm All_output/Trimmed_reads/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
            result = subprocess.run(('rm Analysis_Results/QC_Rawreads/%s*'%(f)), shell=True, capture_output=True, text=True)
            logger.info(result.stdout.rstrip('\n'))
            logger.warning(result.stderr.rstrip('\n'))
        except Exception as e:
            logger.exception(e)
            passed = False
    return True

if __name__ == '__main__':
    t_begin = time.time()
    if not os.path.exists('logs'):
        os.mkdir('logs')
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler('logs/%s'%args.logfile)
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
        del pipeline[11:-1]
    if not args.cleanup:
        del pipeline[-1]
        
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
    