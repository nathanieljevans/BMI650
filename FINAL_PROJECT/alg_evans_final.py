# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 23:40:59 2018

@author: natha


Instructions
Parameter selection is a critical feature in algorithm optimization and protocol standardization / best practices. We will be examining this in the context of alignment for a sample being aligned to a related reference genome. B6 is the reference genome for Mus Musculus. However, a number of other strains are utilized for research. We will be examining sequence data from one of the wild strains, PWK. Our concern is the number of mismatches to select for alignment given the genomic differences between B6 and PWK.  

With regard to the deliverables, you must provide a naive/brute-force analysis of the data 
***
to set a threshold to determine the number of mismatches for PWK based on a 100 base pair read. 
***
You will provide your code and a brief 1-2 page report summarizing your approach, assumptions and limitations.  After you have determined your mismatch threshold, please choose an appropriate aligner (but it must require the user to set the number/threshold of mismatches) to align the provided data. Your written report should also include a discussion of the final alignment percentage.  This is worth 80% of the total project. 

An alternative solution is to design an approach to estimate this from the data. For this path, one strategy would be to design an approach specific to one strain (e.g., PWK). However, the more realistic scenario is to allow this to be generalizable - particularly given the use of “mosaic” crosses such as the Collaborative Cross which utilizes 8 founder strains.  This should be compared to the brute-force approach above. For the remaining 20%, you must provide pseudo-code for either the strain-specific or generalizable strategy.  

The files needed for the project can be found on

state: /home/courses/BMI550/FinalProject

chr1_PWK_PhJ.mgp.v5.snps.dbSNP142.vcf  - SNP annotation for PWK Chromosome 1 (Assembly 38) 

GRCm38_chr1.fa - FASTA file for B6 Chromosome 1 (Assembly 38) 

Mus_musculus.GRCm38.gtf - GTF file for B6 (Assembly 38) 

PWK_R1.fastq  - FastQ file for PWK (RNA-seq data from Illumina Hi-Seq2500)

Each student will give a 3-min presentation on December 5th summarizing their approach and issues they faced (3 slides max). Slides must be submitted by 9am on December 5th so we can load them to allow everyone time to present. 

Project is due December 5th by 9am. 

tldr 
Parameter selection in the context of alignment for PWK (wild strain of mus musculus) being aligned to B6, mus musculus ref. genome. The goal is to find a threshold for the number of mismatches allowed between the genomes. 



Questions for Christina: 
    1) If we choose to follow part 2, is part 1 still necessary? 
    2) what does a number of mismatches for pwk mean 
    
Some Aligners (https://sakai.ohsu.edu/portal/site/BMI-550-1-AS-F18/tool/621a05f7-e013-4047-adb3-b061469b8bc0/ShowPage?returnView=&studentItemId=0&backPath=&errorMessage=&clearAttr=&messageId=&source=&title=&sendingPage=79570&newTopLevel=false&postedComment=false&addBefore=&itemId=377376&path=push&topicId=&addTool=-1&recheck=&id=&forumId=)

Bowtie, Bowtie2, BWA-SW, BWA-MEM

        Aligners: Tophat2, STAR

        HTSeq

        Cufflinks and Cuffdiff

        GTF Format 
    
    
    
Naive Approach 

Use the vcf:
    sliding window to calculated # SNPS in a 100 k-mer window (record in array)
    plot a distribution, choose a 0.05 alpha level and calculate TP,FN->acc and then use gradient descent to best acc. Weight TP,FN for preference. 
    __acc__ = (TP - FN) / (TP + FN)
    
    # Good Mouse entrez info
    https://www.ncbi.nlm.nih.gov/genome?term=mus%20musculus


Questions for Christina: 
    - using entrez, my gene converage is almost 20% of chromosome 1 (B6), are my search results pulling other builds? 
    """
    
import sys
from Bio import Entrez
import pandas as pd
import numpy as np 
import time
import pickle
from os import listdir
from matplotlib import pyplot as plt 
import seaborn as sns
import matplotlib.patches as mpatches

class vcf: 
    B6_chr1_length = 195471971
    
    def __init__(self, raw, init=True): 
        ''' 
        
        '''
        self.variants = np.array([0]*vcf.B6_chr1_length)
        self.full = []
        if(init):
            self.parse(raw)
        
    def parse(self, raw): 
        
        #length of mouse chr 1 = 195.47e6 bases 
        _lines = raw.count('\n')
        for i,variant in enumerate(raw.strip().split('\n')):
            if (i%100000 == 0): 
                print('loading vcf... progress: %.2f' % (i / _lines))
            attr = variant.split('\t')
            __chr,__start,__geneID,__seq1,__seq2 = attr
            
            self.full.append({'chr':__chr, 'start':__start, 's1':__seq1, 's2':__seq2})
            self.variants[int(__start)] = 1 
            
        print('num of snps: %d' %np.sum(self.variants))
            
    def get_kmer_snp_distribution(self, window, re_calculate = False): 
        self.last_window = window
        data = self.__check_prev_and_load(window)
        
        if (not re_calculate and data is not None):
            self.snps_count = data
            
        else: 
            self.snps_count = np.zeros(vcf.B6_chr1_length - window, dtype='byte') 
            _l = len(self.snps_count)
            
            for i in range(0, int(len(self.snps_count))):
                if (i % 10**6 == 0): 
                    print('calculating snps dist... progress: %.3f' % (i/_l))
                    
                self.snps_count[i] = np.sum(self.variants[i:i+window])
                
            with open('./data/snps_distribution- ' + str(window) + '-window.pkl', 'wb') as f: 
                pickle.dump(self.snps_count, f)
            
        self.snps_mean = np.mean(self.snps_count)
        self.snps_std = np.std(self.snps_count)
        self.snps_max = np.max(self.snps_count)
        
        print('average: %f' % (self.snps_mean))
        print('standard dev: %f' %(self.snps_std))
        print('max snps in %d window: %d' %(window, self.snps_max))
        
    
    def __check_prev_and_load(self, window, mypath = './data/'): 
        
        onlypkls = [f for f in listdir(mypath) if (f[-3:] == 'pkl' and f[0:4] == 'snps')]
        
        for f in onlypkls: 
            _window = f.split('-')[1]
            print(_window)
            if (int(_window) == window): 
                print ('loading pre-calculated snps distribution')
                with open('./data/%s' %f, 'rb') as f2: 
                    return pickle.load(f2)
                
        print('no pre-calculated distribution for this window')
        return None 
            
    def plot_snp_distribution(self): 
        
        fig, ax = plt.subplots(figsize=(11,9))
        
        bw_smooth = 1
        sns.distplot(self.snps_count, kde=True, kde_kws={'bw':bw_smooth}, bins=25, norm_hist=True,color='red', hist_kws={'alpha':0.3}, ax=ax)
        
        sns.distplot(self.exonic_snps_dist, kde=True, kde_kws={'bw':bw_smooth}, bins=25, norm_hist=True,color='blue', hist_kws={'alpha':0.3}, ax=ax)
        
        p1 = mpatches.Patch(color='red', label='ALL SNPS')
        p2 = mpatches.Patch(color='blue', label='EXON SNPS')

        plt.xlim(0,30)
        ax.legend(handles=[p1,p2])
        
        #plt.plot([mean, mean], [0, 3], linewidth=2, color='red')
        plt.title('snps probability distribution')
        
    def calculate_accuracy(self, T, dist, l): 
        
        ## Use l to choose between full chromosome or exonic length 
        ## use dist to specify exonic dist or full dist
        l = float(l)
        FN = 0
        FN = np.sum(dist > T)
        TP = len(dist) - FN
        
        ## these are theoretical FP and TN calculations
        k = float(self.last_window)
        TN = l-100 - ((l-100)*(l+100)*(np.power(4.0,T)+k))/(np.power(4.0,k))
        FP = ((l-100)*(np.power(4.0,T) + k)) / (np.power(4.0,k))
        
        acc = (TP-FN) / (TP + FN)
        
        print('T=%.3f --> acc=%.3f (TP: %f, TN: %f, FP: %f, FN: %f)' %(T,acc,TP,TN, FP, FN))
        
        return acc, TP, TN, FP, FN

    def load_exonic_intervals(self): 
        try: 
            with open('./data/ch1_exonic_intervals.pkl', 'rb') as f: 
                self.exonic_intervals = pickle.load(f)
                
        except: 
            raise('could not load exonic intervals (pickled) data, has it been curated yet?')
    
    def get_ch1_exonic_regions(self):
        '''
        This funtion retrieves all the recorded genes on chr1 of Mus Musculus from the Entrez database (Gene) and saves the start, stop and strand for each to an attribute: self.exonic_intervals
        
        inputs
            NONE
        
        outputs
            NONE
        '''
        self.exonic_intervals = []
        Entrez.email = "evansna@ohsu.edu"     # Always tell NCBI who you are
        
        search_handle = Entrez.esearch(db='Gene',term= 'GRCm38.p4[All Fields] AND "Mus musculus"[porgn] AND alive[prop] AND 1[Chr]', retmax=3000)
        
        record = Entrez.read(search_handle)
                
        ids = record["IdList"]
        print('number search results: %d' %len(ids))
                #print('ids', ids)
           
        results = Entrez.efetch(db="Gene", id=','.join(ids), rettype="fasta", retmode="xml")
        all_data = Entrez.read(results)
        _l = len(ids)
        entrez_fails = 0
        for i,data in enumerate(all_data):
            
            try:
                
                _from = data['Entrezgene_comments'][-2]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
                
                _to = data['Entrezgene_comments'][-2]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']
                
                strand = data['Entrezgene_comments'][-2]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_strand']['Na-strand'].attributes['value'] == 'plus' 
                
                gene_int = (_from, _to, strand)
                print(str(gene_int) + '...prog...%f' %(i/_l))
                self.exonic_intervals.append(gene_int)
                
            except: 
                entrez_fails += 1
                print('entrez parsing failure')

        with open('./data/ch1_exonic_intervals.pkl', 'wb') as f: 
            pickle.dump(self.exonic_intervals, f)
                
        print('failures: %d' %entrez_fails)
        print('genes parsed: %d' %(_l-entrez_fails))
    
    def create_chr1_exon_mask(self): 
        
        self.exon_mask = np.full((vcf.B6_chr1_length,2),False, dtype=bool) # row index is strand directionality: 0 -> minus, 1->plus
        
        # each strand has its own row in the array, 
        for _from, _to, _strand in self.exonic_intervals: 
            self.exon_mask[int(_from):int(_to),int(_strand)] = 1 
            
        print(np.mean(self.exon_mask)) # chr1 is almost 20% exon! That seems like too much... it's usually 2% ....
        
        
    def count_exon_snps(self): 
        
        exon_snps = np.transpose(self.exon_mask[:,0] + self.exon_mask[:,1])*( self.variants) 
        print(exon_snps.shape)
        count = np.sum(exon_snps)
        mean = count / np.sum(self.exon_mask)
        
        
        print('exon snps: %d, mean snps / 100 : %.4f' %(count, mean*100))
    
    def get_exon_snp_dist(self, window): 
        
        self.intronic_snps_dist = []
        self.boundary_snps_dist = []
        self.exonic_snps_dist = []
        exonic = self.exon_mask[:,0] + self.exon_mask[:,1]
        for i in range(vcf.B6_chr1_length-window): 
            if (i % int(vcf.B6_chr1_length/1000) == 0): 
                progress(i, vcf.B6_chr1_length, status='calculating exon snp dist...')
                
            c =np.sum(self.variants[i:i+window])
            ex_cnt = np.sum(exonic[i:i+window])
            if (ex_cnt == 0): 
                self.intronic_snps_dist.append(c)
            elif (ex_cnt < 100):
                self.boundary_snps_dist.append(c)
            else:
                self.exonic_snps_dist.append(c)
        
        # THIS ISN"T WORKING, whyyy lisahhhh
        print('saving exon,intron,boundary data to file...')
        with open('./data/exonic_snps_dist_array.pkl', 'wb') as f: 
            pickle.dump(self.exonic_snps_dist, f)
        
        with open('./data/intronic_snps_dist_array.pkl', 'wb') as f: 
            pickle.dump(self.intronic_snps_dist, f)
            
        with open('./data/boundary_snps_dist_array.pkl', 'wb') as f: 
            pickle.dump(self.boundary_snps_dist, f)
            
        print('save complete.')    
        
    def load_exon_snps_dist(self): 
        
        with open('./data/exonic_snps_dist_array.pkl', 'rb') as f: 
           data = pickle.load(f)
           self.exonic_snps_dist = data['exon']
           self.intronic_snps_dist = data['intron']
           self.boundary_snps_dist = data['boundary']
    
    def print_stats(self, dist): 
        
        print('mean: %.3f' %np.mean(dist))
        print('std: %.3f' %np.std(dist))
        print('max: %d' %np.max(dist))
 

# From: https://gist.github.com/vladignatyev/06860ec2040cb497f0f3#file-progress-py 
def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)

       
if __name__ == '__main__' : 
    
    tic = time.clock()
    print('starting...')
    with open('./data/no_header.vcf', 'r') as f: 
        raw = f.read() 
        
    variants = vcf(raw, init=True)
    
    # This can be recurated or loaded from a pickled file 
    variants.load_exonic_intervals()
    #variants.get_ch1_exonic_regions() # 1916 genes curated / 809 parse fails - should be only ~1600 , almost 20% exon coverage... 10x expected! Should add some QC for the genes
    
    variants.create_chr1_exon_mask()
    
    variants.count_exon_snps()
    
    # to recalculate snps/window distribution
    variants.get_exon_snp_dist(window=100)
    # to load from disk the variants 
    #variants.load_exon_snps_dist()
    
    # gets the full snps distribution, no distinction between exon/intron
    variants.get_kmer_snp_distribution(100, re_calculate=False)
    
    print('exon dist')
    variants.print_stats(variants.exonic_snps_dist)
    print('\n\n intron dist')
    variants.print_stats(variants.intronic_snps_dist)
    print('\n\n boundary dist')
    variants.print_stats(variants.boundary_snps_dist)
    print('\n\n full dist')
    variants.print_stats(variants.snps_count)


    '''    
    variants.plot_snp_distribution()
    
    X = np.arange(1,87,1)
    acc = np.zeros(len(X))
    TP = np.zeros(len(X))
    TN = np.zeros(len(X))
    FP = np.zeros(len(X))
    FN = np.zeros(len(X))
    
    exonic_length = np.sum(variants.exon_mask)
    for i,x in enumerate(X): 
        acc[i], TP[i], TN[i], FP[i], FN[i] = variants.calculate_accuracy(T=x, dist=variants.exonic_snps_dist, l=exonic_length)
        
    df = pd.DataFrame({'X':X, 'acc':acc, 'TP':TP, 'TN':TN, 'FP':FP, 'FN':FN, 'sensitivity':TP/(TP+FN), 'specificity':TN/(TN+FP), 'precision':TP/(TP+FP)})
    
    df.to_csv('./outputs/threshold_comparison-exon.csv')
    
    fig, ax = plt.subplots(figsize=(15,9))
    plt.gcf()
    ax.plot('X','TP', data = df, color='green')
    ax.plot('X', 'TN', data = df, color='blue')
    ax.plot('X', 'FP', data = df, color='purple')
    ax.plot('X', 'FN', data = df, color='red')
    p1 = mpatches.Patch(color='green', label='True Positive')
    p2 = mpatches.Patch(color='blue', label='True Negative')
    p3 = mpatches.Patch(color='purple', label='False Positive')
    p4 = mpatches.Patch(color='red', label='False Negative')
    ax.legend(handles=[p1,p2,p3,p4])
    plt.title('Threshold vs outcomes for Exonic Regions')
    plt.show()
    
    fig2, ax2 = plt.subplots(figsize=(15,9))
    plt.xlim(0,15)
    plt.gcf()
    ax2.plot('X', 'acc', data = df, color='red')
    ax2.plot('X', 'sensitivity', data = df, color='blue')
    ax2.plot('X', 'specificity', data = df, color='green')
    ax2.plot('X', 'precision', data = df, color='purple')
    p1 = mpatches.Patch(color='red', label='accuracy')
    p2 = mpatches.Patch(color='blue', label='sensitivity')
    p3 = mpatches.Patch(color='green', label='specificity')
    p4 = mpatches.Patch(color='purple', label='precision')
    ax2.legend(handles=[p1,p2,p3,p4])
    plt.title('Threshold vs outcome metric for Exonic Regions')
    plt.show()
    
    '''  

    toc = time.clock()
    print('finished. Time to execute: %.3f seconds' %(toc-tic))





















