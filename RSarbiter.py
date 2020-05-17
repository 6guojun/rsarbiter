#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 14:42:04 2020

@author: Liu
"""

import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import neatbio.sequtils as utils
matplotlib.use("Agg")


def delta(x,y):
    return 0 if x == y else 1

def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))

def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)

def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice

def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()
# def gc_content(seq):
#     result = float(str(seq).count('G') + str(seq).count('C'))/len(seq) * 100
#     return result
# def at_content(seq):
#     result = float(str(seq).count('A') + str(seq).count('T'))/len(seq) * 100
#     return result






#
def main():
    """A Bioinformatics App used to visuallize recombination site"""
    st.title("RSarbiter: a machine learning supported expert in prediction and analysis of recombination spots in S. cerevisiae")
    menu = ["Intro","DNA Sequence","Dot Plot","About"]
    choice = st.sidebar.selectbox("Select Activity", menu)
    if choice == "Intro":
        st.subheader("Intro to BioInformatics")
    elif choice == "DNA Sequence":
        st.subheader("DNA Sequence Analysis")
        seq_file = st.file_uploader("Upload FASTA File",type = ["fasta","fa","txt"])
        if seq_file is not None:
            dna_record = SeqIO.read(seq_file,"fasta")
            # dna_record = list(SeqIO.parse(seq_file,"fasta"))
            dna_seq = dna_record.seq
            st.write(dna_record)
            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record.description)
            else:
                st.write(dna_record.seq)
            # Nucleotide Frequency
            st.subheader("Nucleotide Frequency")
            dan_freq = Counter(dna_seq)
            st.write(dan_freq)
            adenine_color = st.beta_color_picker("Adenine Color")
            thymine_color = st.beta_color_picker("thymine Color")
            cytosine_color = st.beta_color_picker("cytosine Color")
            guanine_color = st.beta_color_picker("guanine Color")
            if st.button("Plot Freq"):
                barlist = plt.bar(dan_freq.keys(),dan_freq.values())
                barlist[0].set_color(adenine_color)
                barlist[1].set_color(thymine_color)
                barlist[2].set_color(guanine_color)
                barlist[3].set_color(cytosine_color)
                st.pyplot()
            st.subheader("DNA composition")
            gc_content = utils.gc_content(dna_seq)
            at_content = utils.at_content(dna_seq)
            st.write({"GC Content":gc_content,"AT Content":at_content}) # st.write = st.json 

            # Nucleotide Count
            nt_count = st.text_input("Enter Nucleotide Here","Type Nucleotide Alphabet")
            st.write("Number of {} Nucleotide is {}" .format((nt_count),str(dna_seq).count(nt_count)))

            # Protein Synthesis
            st.subheader("Protein Synthesis")
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))
            if st.checkbox("Transcription"):
                st.write(dna_seq.transcribe())
            # Top Most Common Amino
            elif st.checkbox("Translation"):
                st.write(dna_seq.translate())
            elif st.checkbox("Complement"):
                st.write(dna_seq.complement())
            elif st.checkbox("AA Frequency"):
                st.write(aa_freq)
            elif st.checkbox("Plot AA Frequency"):
                aa_color = st.beta_color_picker("Pick An Amino Acid Color")
                # barlist = plt.bar(aa_freq.keys(),aa_freq.values())
                # barlist[2].set_color(aa_color)
                plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                st.pyplot()
            elif st.checkbox("Full Amino Acid Name"):
                aa_name = str(p1).replace("*","")
                aa3=utils.convert_1to3(aa_name)
                st.write(aa_name)
                st.write("======================")
                st.write(aa3)
                st.write("======================")
                st.write(utils.get_acid_name(aa3))





    elif choice == "Dot Plot":
        st.subheader("Generate a Dot Plot For Two Sequences")
        # seq_file1 = st.file_uploader("Upload 1st FASTA File",type=["fasta","fa"])
        # seq_file2 = st.file_uploader("Upload 2nd FASTA File",type=["fasta","fa"])

        # if seq_file1 and seq_file2 is not None:
        #     dna_record1 = SeqIO.read(seq_file1,"fasta")
        #     dna_record2 = SeqIO.read(seq_file2,"fasta")
        #     # st.write(dna_record)
        #     dna_seq1 = dna_record1.seq
        #     dna_seq2 = dna_record2.seq

        #     details = st.radio("Details",("Description","Sequence"))
        #     if details == "Description":
        #         st.write(dna_record1.description)
        #         st.write("=====================")
        #         st.write(dna_record2.description)
        #     elif details == "Sequence":
        #         st.write(dna_record1.seq)
        #         st.write("=====================")
        #         st.write(dna_record2.seq)


        #     cus_limit = st.number_input("Select Max number of Nucleotide",10,200,50)
        #     if st.button("Dot Plot"):
        #         st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
        #         dotplotx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit])

        #         st.pyplot()




        seq_file_1 = st.file_uploader("Upload 1st FASTA File",type = ["fasta","fa","txt"])
        seq_file_2 = st.file_uploader("Upload 2nd FASTA File",type = ["fasta","fa","txt"])

        if seq_file_1 and seq_file_2 is not None:
            dna_record_1 = SeqIO.read(seq_file_1,"fasta")
            dna_record_2 = SeqIO.read(seq_file_2,"fasta")
            dna_seq_1 = dna_record_1.seq
            dna_seq_2 = dna_record_2.seq

            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record_1.description)
                st.write("======================")
                st.write(dna_record_2.description)
            elif details == "Sequence":
                st.write(dna_record_1.seq)
                st.write("==============")
                st.write(dna_record_2.seq)
            cus_limit = st.number_input("Select Max Number of Nucleotide", 10, 200, 50)
            if st.button("Dot Plot"):
                st.write("Comparing The First {} Nucleotide of The Two Sequence".format(cus_limit))
                dotplotx(dna_seq_1[0:cus_limit],dna_seq_2[0:cus_limit])
                st.pyplot()


    elif choice == "About":
        st.subheader("About")

if __name__ == '__main__':
    main()
    
#help(SeqIO.read)
