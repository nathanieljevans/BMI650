# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 23:40:59 2018

@author: nathaniel evans
@class: BMI650 
@HW: FINAL PROJECT

Instructions
This program can be run from a terminal but must be in the following file structure: 
    
    alg_evans_final.py
    /data/ 
        no_header.vcf
    /outputs/
    
Additionally, if there the script has not been run previously it is necessary to change commenting in main to run full script. 
Future runs can be sped up significantly by loading pickled files from disk. 

This algorithm focuses on using the variation of the exonic regions of chromosome 1 of the ref genome B6 (mus musculus) to predict 
an acceptable number of mismatches (T) to align RNA-seq data from PWK strain to B6 genome. 

Keep in mind this is RNA-seq data, so it is processed mRNA chunks. Only exon regions will be included. Fortunately, aligners 
such as STAR can include splicing into the alignment so T mismatches only refers to the mismatches within exonic regions. 
    """
    
from sklearn.externals import joblib
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
    B6_chr1_length = 195471971 # length of chr1 B6 nucleotide sequence 
    
    def __init__(self, raw, init=True): 
        ''' 
        initialize the vcf object (name is a little deprecated now) 
        
        input 
            raw <str> .vcf file to parse 
            init <boolean> option to parse, useful to speed up testing if the vcf parsing isn't necessary 
        '''
        self.variants = np.array([0]*vcf.B6_chr1_length)
        self.full = []
        if(init):
            self.parse(raw)
        
    def parse(self, raw): 
        '''
        parse the .vcf data representing SNPs in chromosome 1 of mus musculus B6 reference genome. 
        
        inputs 
            raw <str> the vcf file to parse 
            
        outputs 
            None 
        '''
        
        #length of mouse chr 1 = 195.47e6 bases 
        _lines = raw.count('\n')
        for i,variant in enumerate(raw.strip().split('\n')):
            if (i%100000 == 0): 
                progress(i,_lines,'loading vcf')
                #print('loading vcf... progress: %.2f' % (i / _lines))
            attr = variant.split('\t')
            __chr,__start,__geneID,__seq1,__seq2 = attr
            
            self.full.append({'chr':__chr, 'start':__start, 's1':__seq1, 's2':__seq2})
            self.variants[int(__start)] = 1 
            
        print('\nnum of snps in chr1: %d' %np.sum(self.variants))
            
    def get_kmer_snp_distribution(self, window, re_calculate = False): 
        '''
        calculate the k-mer SNP distribution for B6 chromosome 1 using the full region. 
        
        inputs 
            window <int> size of k-mer 
            re_calculate <boolean> option to load from memory if available
            
        outputs 
            None 
        '''
        
        self.last_window = window
        data = self.__check_prev_and_load(window)
        
        if (not re_calculate and data is not None):
            self.snps_count = data
            
        else: 
            self.snps_count = np.zeros(vcf.B6_chr1_length - window, dtype='byte') 
            _l = len(self.snps_count)
            
            for i in range(0, int(len(self.snps_count))):
                if (i % 10**6 == 0): 
                    progress(i, _l, 'calculating full snps')

                    
                self.snps_count[i] = np.sum(self.variants[i:i+window])
                
            with open('./data/snps_distribution- ' + str(window) + '-window.pkl', 'wb') as f: 
                pickle.dump(self.snps_count, f)
            
        self.snps_mean = np.mean(self.snps_count)
        self.snps_std = np.std(self.snps_count)
        self.snps_max = np.max(self.snps_count)
        
        #print('\naverage: %f' % (self.snps_mean))
        #print('standard dev: %f' %(self.snps_std))
        #print('max snps in %d window: %d' %(window, self.snps_max))
        
    
    def __check_prev_and_load(self, window, mypath = './data/'): 
        '''
        load the ALL regions SNPs distribution into memory if it's available. 
        
        inputs
            window <int> length of k-mer for sliding window
            mypath <string> where to look for the file 
            
        outputs 
            SNPs distribution <np.array> 
        '''
        
        onlypkls = [f for f in listdir(mypath) if (f[-3:] == 'pkl' and f[0:4] == 'snps')]
        
        for f in onlypkls: 
            _window = f.split('-')[1]
            if (int(_window) == window): 
                print ('loading pre-calculated snps distribution')
                with open('./data/%s' %f, 'rb') as f2: 
                    return pickle.load(f2)
                
        print('no pre-calculated distribution for this window')
        return None 
            
    def plot_snp_distribution(self): 
        '''
        plot the cumulative distribution function of the SNPs / k-mer
        
        inputs 
            None
        
        outputs 
            None
        '''
        
        fig, ax = plt.subplots(figsize=(11,9))
        
        bw_smooth = .5
        sns.distplot(self.snps_count, kde=True, hist=False, kde_kws={'bw':bw_smooth}, bins=25, norm_hist=True,color='red', hist_kws={'alpha':0.3}, ax=ax)
        
        sns.distplot(self.exonic_snps_dist, kde=True, hist=True, kde_kws={'bw':bw_smooth}, bins=20, norm_hist=True,color='blue', hist_kws={'alpha':0.3}, ax=ax)
        
        sns.distplot(self.intronic_snps_dist, kde=True, hist=False, kde_kws={'bw':bw_smooth}, bins=25, norm_hist=True,color='green', hist_kws={'alpha':0.3}, ax=ax)

        sns.distplot(self.boundary_snps_dist, kde=True, hist=False, kde_kws={'bw':bw_smooth}, bins=25, norm_hist=True,color='purple', hist_kws={'alpha':0.3}, ax=ax)

        p1 = mpatches.Patch(color='red', label='ALL SNPS')
        p2 = mpatches.Patch(color='blue', label='EXON SNPS')
        p3 = mpatches.Patch(color='green', label='INTRON SNPS')
        p4 = mpatches.Patch(color='purple', label='BOUNDARY SNPS')

        plt.xlim(0,10)
        ax.legend(handles=[p1,p2,p3,p4])
        
        plt.xlabel('# SNPs in 100-mer')
        plt.ylabel('probability')
        plt.title('cumulative distribution of variation in 100-mers of B6 chromosome 1')
        fig.savefig('./outputs/SNPs_distribution.png')
        plt.show()
        
    def calculate_accuracy(self, T, dist, l): 
        '''
        calculate the True Positives and False Negatives given T and our SNPs / k-mer distributon. 
        Estimate the True Negatives and False Positives using a probability model. 
        
        inputs
            T <int> number of acceptable mis-matches / SNPs 
            dist <np.array> each index represents a k-mer and the value is the # of SNPs in that window
            l <int> length of the region of interest (eg chr1, exon-chr1...etc )
            
        outputs 
            None 
        '''
        ## Use l to choose between full chromosome or exonic length 
        ## use dist to specify exonic dist or full dist
        l = float(l)
        FN = np.sum(dist > T)
        TP = len(dist) - FN
        #print ('TP %.3f' %TP)
        
        ## these are theoretical FP and TN calculations
        k = float(self.last_window)
        
        # Rough estimate model
        #P = np.power(4.0, T) / np.power(4.0, k)
        TN = l-100 - ((l-100)*(l+100)*(np.power(4.0,T)+k))/(np.power(4.0,k))
        FP = ((l-100)*(np.power(4.0,T) + k)) / (np.power(4.0,k))
        
        # Factorial model
        # np.math.factorial(100) - np.math.factorial(100-T) -> 
        #partial_factorial = np.arange(k, k-T-1, -1)
        #P = (np.power(4.0,T) * (np.product(partial_factorial))) / np.power(4.0,k)
        #print('P %.3f' %P)
        
        #TN = (1.0-P)*(l-k)
        #FP = (P) * (l-k)
        
        #print('T=%.3f --> acc=%.3f (TP: %f, TN: %f, FP: %f, FN: %f)' %(T,acc,TP,TN, FP, FN))
        
        return TP, TN, FP, FN

    def load_exonic_intervals(self): 
        '''
        loads the exon intervals from disk. Saves a ton of time not to re-calculate. 
        
        inputs
            None
            
        outputs 
            None
        '''
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
        
        search_handle = Entrez.esearch(db='Gene',term= '((GRCm38.p4[Assembly Name]) AND 1[Chromosome]) AND Mus Musculus[Organism] AND srcdb_refseq_known[PROP]', retmax=3000)
        
        record = Entrez.read(search_handle)
                
        ids = record["IdList"]
        print('number search results: %d' %len(ids))
                #print('ids', ids)
           
        results = Entrez.efetch(db="Gene", id=','.join(ids), rettype="fasta", retmode="xml")
        all_data = Entrez.read(results,validate=False)
        _l = len(ids)
        entrez_fails = 0
        for i,data in enumerate(all_data):
            
            try:
                progress(i, _l, 'extracting exonic intervals [fail: %s]' %entrez_fails)
                _from = data['Entrezgene_comments'][-2]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
                
                _to = data['Entrezgene_comments'][-2]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']
                
                strand = data['Entrezgene_comments'][-2]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_comment'][-1]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_strand']['Na-strand'].attributes['value'] == 'plus' 
                
                gene_int = (_from, _to, strand)
                #print(str(gene_int) + '...prog...%f' %(i/_l))
                self.exonic_intervals.append(gene_int)
                
            except: 
                entrez_fails += 1

        with open('./data/ch1_exonic_intervals.pkl', 'wb') as f: 
            pickle.dump(self.exonic_intervals, f)
                
        print('\nfailures: %d' %entrez_fails)
        print('genes parsed: %d' %(_l-entrez_fails))
    
    def create_chr1_exon_mask(self):
        '''
        Create the exon mask representing B6 chromosome 1
        
        inputs
            None
        
        outputs 
            None
        '''
        
        self.exon_mask = np.full((vcf.B6_chr1_length,2),False, dtype=bool) # row index is strand directionality: 0 -> minus, 1->plus
        
        # each strand has its own row in the array, 
        for _from, _to, _strand in self.exonic_intervals: 
            self.exon_mask[int(_from):int(_to),int(_strand)] = 1 
            
        #print(np.mean(self.exon_mask)) # chr1 is almost 20% exon! That seems like too much... it's usually 2% .... ITS BECAUSE THE INTRONS WERN"T SPLICED OUT
        
        
    def count_exon_snps(self): 
        '''
        counts the SNPs in the exon regions
        
        inputs
            None
        
        outputs
            None
        '''
        
        exon_snps = np.transpose(self.exon_mask[:,0] + self.exon_mask[:,1])*( self.variants) 
        print(exon_snps.shape)
        count = np.sum(exon_snps)
        mean = count / np.sum(self.exon_mask)
        
        print('exon snps: %d, mean snps / 100 : %.4f' %(count, mean*100))
    
    def get_exon_snp_dist(self, window): 
        '''
        Calculate the exon SNP distribution. Loop through the SNPs mask and count SNPs per window. 
        
        input 
            window <int> size of k-mer to use for the sliding window
            
        output
            pickled file to disk
        '''
        self.intronic_snps_dist = []
        self.boundary_snps_dist = []
        self.exonic_snps_dist = []
        exonic = self.exon_mask[:,0] + self.exon_mask[:,1]
        for i in range(vcf.B6_chr1_length-window): 
            if (i % int(vcf.B6_chr1_length/1000) == 0): 
                progress(i, vcf.B6_chr1_length, status='calculating exon snp dist')
                
            c =np.sum(self.variants[i:i+window])
            ex_cnt = np.sum(exonic[i:i+window])
            if (ex_cnt == 0): 
                self.intronic_snps_dist.append(c)
            elif (ex_cnt < 100):
                self.boundary_snps_dist.append(c)
            else:
                self.exonic_snps_dist.append(c)
        
        self.intronic_snps_dist = np.array(self.intronic_snps_dist, dtype='byte')
        self.boundary_snps_dist = np.array(self.boundary_snps_dist, dtype='byte')
        self.exonic_snps_dist = np.array(self.exonic_snps_dist, dtype='byte')
        
        # THIS ISN"T WORKING, whyyy lisahhhh
        
        try:
            print('\nsaving exon,intron,boundary data to file...')
            with open('./data/exonic_snps_dist_array.pkl', 'wb') as f: 
                pickle.dump(self.exonic_snps_dist, f)
            print('first!')
           
            with open('./data/intronic_snps_dist_array.pkl', 'wb') as f: 
                pickle.dump(self.intronic_snps_dist, f)
             
            print('second!')
            with open('./data/boundary_snps_dist_array.pkl', 'wb') as f: 
                pickle.dump(self.boundary_snps_dist, f)
            print('sweet relief!')
        except:
            try: 
                joblib.dump(self.exonic_snps_dist, './data/exonic_snps_dist_array.pkl')
                joblib.dump(self.intronic_snps_dist, './data/intronic_snps_dist_array.pkl')
                joblib.dump(self.boundary_snps_dist, './data/boundary_snps_dist_array.pkl')
            except:
                print(':(')
                raise
        
        print('save complete.')    
        
    def load_exon_snps_dist(self): 
        '''
        loads exon SNPs distribution from disk (pickled file)
        
        inputs 
            None
        
        outputs 
            None
        '''
        with open('./data/exonic_snps_dist_array.pkl', 'rb') as f: 
           self.exonic_snps_dist = pickle.load(f)
           
        with open('./data/intronic_snps_dist_array.pkl', 'rb') as f: 
           self.intronic_snps_dist = pickle.load(f)
           
        with open('./data/boundary_snps_dist_array.pkl', 'rb') as f: 
           self.boundary_snps_dist = pickle.load(f)
           
    
    def print_stats(self, dist): 
        '''
        prints the stats of a given distribution 
        
        input 
            dist <numpy array> each index represents a 100-mer window and the value is a count of SNPs
        
        output
            None
        '''
        print('mean: %.3f' %np.mean(dist))
        print('std: %.3f' %np.std(dist))
        print('max: %d' %np.max(dist))
 

def progress(count, total, status=''):
    '''
    progress bar for terminal excecution, taken from: https://gist.github.com/vladignatyev/06860ec2040cb497f0f3#file-progress-py 
    '''
    
    bar_len = 20
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
    #variants.get_ch1_exonic_regions() # <-- do this if no pickled files.    1916 genes curated / 809 parse fails - should be only ~1600 , almost 20% exon coverage... 10x expected! Should add some QC for the genes... NOPE it's the gene introns that haven't been spliced out
    
    variants.create_chr1_exon_mask()
    
    variants.count_exon_snps()
    
    # to recalculate snps/window distribution
    #variants.get_exon_snp_dist(window=100) # <---- do this if no pickled files 
    # to load from disk the variants 
    variants.load_exon_snps_dist()
    
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

    variants.plot_snp_distribution()
    
    X = np.arange(1,90,1)
    TP = np.zeros(len(X))
    TN = np.zeros(len(X))
    FP = np.zeros(len(X))
    FN = np.zeros(len(X))
    
    exonic_length = np.sum(variants.exon_mask)
    for i,x in enumerate(X): 
        progress(i, len(X), 'calculating sensitivity and specificity')
        TP[i], TN[i], FP[i], FN[i] = variants.calculate_accuracy(T=x, dist=variants.exonic_snps_dist, l=exonic_length)

    df = pd.DataFrame({'X':X,'TP':TP, 'TN':TN, 'FP':FP, 'FN':FN, 'sensitivity':TP/(TP+FN), 'specificity':TN/(TN+FP), 'precision':TP/(TP+FP)})
    
    # this is what I use to look at sensitivity 
    df.to_csv('./outputs/threshold_comparison-exon.csv')
    
    #plotting graphics 
    fig, ax = plt.subplots(figsize=(15,9))
    plt.ylim(-10, 0.6e8)
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
    plt.xlabel('# mismatches (T)')
    plt.ylabel('count')
    plt.savefig('./outputs/TP_TN_plot.png')
    plt.show()
    
    fig2, ax2 = plt.subplots(figsize=(15,9))
    plt.xlim(0,90)
    plt.ylim(-.25,1.25)
    plt.gcf()
    ax2.plot('X', 'sensitivity', data = df, color='blue')
    ax2.plot('X', 'specificity', data = df, color='green')
    #ax2.plot('X', 'precision', data = df, color='purple')
    p2 = mpatches.Patch(color='blue', label='sensitivity')
    p3 = mpatches.Patch(color='green', label='specificity')
    #p4 = mpatches.Patch(color='purple', label='precision')
    ax2.legend(handles=[p3,p4])
    plt.xlabel('Threshold (T)')
    plt.ylabel('Probability')
    plt.title('Threshold vs Outcome metric for Exonic Regions')
    plt.savefig('./outputs/sens_spec_T.png')
    plt.show()

    toc = time.clock()
    print('finished. Time to execute: %.3f seconds' %(toc-tic))





















