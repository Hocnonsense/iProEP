# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 15:51:56 2017

@author: laihongyan
"""
import copy
import math
import os
import sys


def Kmers(k):
    kmers = [] 
    oligos = ['A','G','C','T']
    if k == 1:
        kmers = oligos
        return kmers  
    k_1mers = Kmers(k-1) 
    for oligo in oligos:
        for k_1mer in k_1mers:
            kmer = k_1mer + oligo
            kmers.append(kmer)
    return kmers

#构建整个位置权重矩阵（allsiteNum=82-k个位点）
def PWM(sequenceFile,kmers,allsiteNum,p0):
    pwm = []
    k = len(kmers[0])
    realCounts = []
    for i in range(allsiteNum):
        name = {}
        for kmer in kmers:
            name[kmer] = 0
        realCounts.append(name)
    pwm = copy.deepcopy(realCounts)
    totalCounts = []
    psetotalCounts = []
    ##############count the number of kmers at the corresponding allsites
    fobjr = open(sequenceFile,'r')
    for eachline in fobjr:
        if eachline[0] == '>':
            pass
        else:
            eachline = eachline.strip()
            length = len(eachline)
            for i in range(length+1-k):
                kmer = eachline[i:i+k]
                realCounts[i][kmer] += 1                          
    fobjr.close()
    for i in range(allsiteNum):
        totalCounts.append(sum(realCounts[i].values()))    ###直接对字典的所有values进行求和
        #print("the total number of real counts at the %d site is %d" % (i,sum(realCounts[i].values())))
        psetotalCounts.append(math.sqrt(totalCounts[i]))   ###math.sqrt()开平方
        #print("psetotalCounts at the %d site is %d" % (i,psetotalCounts[i]))
    for i in range(allsiteNum):
        for kmer in kmers:
            pwm[i][kmer] = (realCounts[i][kmer] + p0 * psetotalCounts[i])/(totalCounts[i]+ psetotalCounts[i])
    #print('the size of pwm is %d' % (len(pwm)))
    return pwm

def PCSF(sequences,kmers,p0,pwm,sites):
    pcsf = []
    k = len(kmers[0])
    for eachline in sequences:
        if eachline[0] == '>':
            pass
        else:
            eachline = eachline.strip()
            pcsfI = []
            length = len(eachline)
            for i in range(length+1-k):
                if i in sites:
                    kmer = eachline[i:i+k]
                    pcsfI.append(math.log(pwm[i][kmer]/float(p0),math.e))# 更改 float
                else:
                    pass             
            pcsf.append(pcsfI)
    #print('the size of pcsf is %d' % (len(pcsf)))
    return pcsf

def writePcsfFea(pcsf,num,pcsfFeaFile):    
    fobjw = open('.\\5model\\ecmodel\\'+pcsfFeaFile,'w')
    for i in range(len(pcsf)):     
        for j in range(len(pcsf[i])):
            temp = num + j        
            fobjw.write('%d:%f\t' % (temp,pcsf[i][j]))
        fobjw.write('\n')		
    fobjw.close()
        
def getPCSF(sequences,pcsfFeaFile):
    k = 3
    kmers = Kmers(k)
    allsiteNum = 81-k+1 
    p0 = 1/float(4**k)#define float
    sites = [9,23,24,25,26,44,45,46,47,48,49,50,51,52,53,58,59]
    pwmP = PWM('.\\5model\\ecmodel\\'+'promoter.txt',kmers,allsiteNum,p0)
    pcsf = PCSF(sequences,kmers,p0,pwmP,sites) 
    writePcsfFea(pcsf,329,pcsfFeaFile)

"""
sequences = ['AGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAATAAAGAGAGA','CGGGCCTATAAGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAA','TATGTAACATAATGCGACCAATAATCGTAATGAATATGAGAAGTGTGATATTATAACATTTCATGACTACTGCAAGACTAA','TCGCACGGGTGGATAAGCGTTTACAGTTTTCGCAAGCTCGTAAAAGCAGTACAGTGCACCGTAAGAAAATTACAAGTATAC']
getPCSF(sequences,'pcsfFeaFile.txt')
"""


