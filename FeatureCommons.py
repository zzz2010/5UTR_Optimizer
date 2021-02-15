####extract feature for latter prediction use
###RNA folding
###Condon
### k-mer
### motif
import sys,os
from Bio import SeqIO
import Bio.SeqUtils.CodonUsage
import subprocess
from multiprocessing import Pool
import gzip
from Bio.Seq import Seq

cds_length=15 ##assumption last portion of sequence is cds



def codonFreq(seq):
        codon_str=seq.translate().tostring()        
        tot=len(codon_str)
        feature_map=dict()
        for a in codon_str:
                a="codon_"+a
                if a not in feature_map:
                        feature_map[a]=0
                feature_map[a]+=1.0/tot
        feature_map['uAUG']=codon_str.count("M") #number of start codon
        feature_map['uORF']=codon_str.count("*") #number of stop codon
        return feature_map

def singleNucleotide_composition(seq):
        dna_str=seq.tostring().upper()
        N_count=dict()  #add one pseudo count
        N_count['C']=1
        N_count['G']=1
        N_count['A']=1
        N_count['T']=1
        for a in dna_str:
                if a not in N_count:
                        N_count[a]=0
                N_count[a]+=1
        feature_map=dict()
        feature_map["CGperc"]=float(N_count['C']+N_count['G'])/len(dna_str)
        feature_map['CGratio']=abs(float(N_count['C'])/N_count['G']-1)
        feature_map['ATratio']=abs(float(N_count['A'])/N_count['T']-1)
        feature_map['utrlen_m80']=abs(len(dna_str)-80-cds_length)        
        return feature_map


def RNAfold_energy(sequence, *args):
        rnaf = subprocess.Popen(["RNAfold","--noPS"] + list(args),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                # Universal Newlines effectively allows string IO.
                                universal_newlines=True)
        rnafold_output, folderr = rnaf.communicate(sequence)
        output_lines = rnafold_output.strip().splitlines()
        sequence = output_lines[0]
        structure = output_lines[1].split(None,1)[0].strip()
        energy = float(output_lines[1].rsplit("(",1)[1].strip("()").strip())
        return energy

def RNAfold_energy_Gquad(sequence, *args):
        rnaf = subprocess.Popen(["RNAfold","--noPS"] + list(args),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                # Universal Newlines effectively allows string IO.
                                universal_newlines=True)
        rnafold_output, folderr = rnaf.communicate(sequence)
        output_lines = rnafold_output.strip().splitlines()
        sequence = output_lines[0]
        structure = output_lines[1].split(None,1)[0].strip()
        energy = float(output_lines[1].rsplit("(",1)[1].strip("()").strip())
        return energy


def foldenergy_feature(seq):
        dna_str=seq.tostring()
        feature_map=dict()
        feature_map['energy_5cap']=RNAfold_energy(dna_str[:100])
        feature_map['energy_whole']=RNAfold_energy(dna_str)
        feature_map['energy_last30bp']=RNAfold_energy(dna_str[(len(dna_str)-30):len(dna_str)])
        feature_map['energy_Gquad_5utr']=RNAfold_energy_Gquad(dna_str[:(len(dna_str)-15)])
        feature_map['energy_Gquad_5cap']=RNAfold_energy_Gquad(dna_str[:50])
        feature_map['energy_Gquad_last50bp']=RNAfold_energy_Gquad(dna_str[(len(dna_str)-50):len(dna_str)])        
        return feature_map

def Kmer_feature(seq,klen=6):
        feature_map=dict()
        seq=seq.upper()
        for k in range(1,klen+1):
                for st in range(len(seq)-klen):
                        kmer=seq[st:(st+k)]
                        featname="kmer_"+kmer.tostring()
                        if featname not in feature_map:
                                feature_map[featname]=0
                        feature_map[featname]+=1.0/(len(seq)-k+1)
        return feature_map

def oss(cmd):
        print(cmd)
        os.system(cmd)


##seq is seq object from bio.python
def Seq2Feature(seq):
        ##codon
        ret=list(codonFreq(seq).items())
        ##DNA CG composition
        ret+=list(singleNucleotide_composition(seq).items())
        ##RNA folding
        ret+=list(foldenergy_feature(seq).items())
        ret+=list(Kmer_feature(seq).items())
        return ret




