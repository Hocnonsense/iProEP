# -*- coding: utf-8 -*-
"""
 * @Date: Wed Mar 29 15:51:56 2017
 * @author: laihongyan
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-12 16:54:07
 * @FilePath: /iProEP_localtool/cemodel/PCSF.py
 * @Description:
"""

import copy
import math
from typing import Final
from collections import UserList
from Bio import SeqIO


class Kmers(UserList):
    oligos: Final = ["A", "G", "C", "T"]

    def __init__(self, k):
        self.k = k
        if k == 1:
            self.data = self.oligos
        else:
            self.data = [i + j for j in type(self)(k - 1) for i in self.oligos]


# 构建整个位置权重矩阵（allsiteNum=82-k个位点）
def PWM(sequenceFile: str, sequenceLen: int, kmers: Kmers, p0: float):
    allsiteNum = sequenceLen - kmers.k + 1

    realCounts: list[dict[str, int]] = [
        {kmer: 0 for kmer in kmers} for _ in range(allsiteNum)
    ]

    ##############count the number of kmers at the corresponding allsites
    for seq in (str(i.seq) for i in SeqIO.parse(sequenceFile, "fasta")):
        for i in range(len(seq) + 1 - kmers.k):
            realCounts[i][seq[i : i + kmers.k]] += 1

    ###直接对字典的所有values进行求和
    totalCounts = [sum(site.values()) for site in realCounts]
    # print("the total number of real counts at the %d site is %d" % (i,sum(realCounts[i].values())))
    ###math.sqrt()开平方
    psetotalCounts = [math.sqrt(site) for site in totalCounts]
    # print("psetotalCounts at the %d site is %d" % (i,psetotalCounts[i]))

    pwm = [
        {rk: (rv + p0 * p) / (t + p) for rk, rv in r.items()}
        for r, p, t in zip(realCounts, psetotalCounts, totalCounts)
    ]
    # print('the size of pwm is %d' % (len(pwm)))
    return pwm


def PCSF(sequences: list[str], kmers: Kmers, p0, pwm, sites):
    pcsf = []
    for eachline_ in sequences:
        if eachline_.startswith(">"):
            continue
        eachline = eachline_.strip()
        pcsfI = []
        length = len(eachline)
        for i in range(length + 1 - kmers.k):
            if i in sites:
                kmer = eachline[i : i + kmers.k]
                pcsfI.append(math.log(pwm[i][kmer] / float(p0), math.e))  # 更改float
            else:
                pass
        pcsf.append(pcsfI)
    # print('the size of pcsf is %d' % (len(pcsf)))
    return pcsf


def writePcsfFea(pcsf, num, pcsfFeaFile):
    with open(pcsfFeaFile, "w") as fobjw:
        for i in range(len(pcsf)):
            for j in range(len(pcsf[i])):
                temp = num + j
                fobjw.write(f"{temp}:{pcsf[i][j]:f}\t")
            fobjw.write("\n")


def getPCSF(sequences: list[str], pcsfFeaFile: str):
    k = 3
    kmers = Kmers(k)
    p0 = 1 / float(4**k)  # define float
    sites = [
        141,
        162,
        195,
        200,
        207,
        208,
        209,
        210,
        212,
        238,
        246,
        247,
        248,
        249,
        250,
        251,
        252,
    ]
    pwmP = PWM("promoter.txt", 300, kmers, p0)
    pcsf = PCSF(sequences, kmers, p0, pwmP, sites)
    writePcsfFea(pcsf, 389, pcsfFeaFile)


"""
sequences = ['AGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAATAAAGAGAGA','CGGGCCTATAAGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAA','TATGTAACATAATGCGACCAATAATCGTAATGAATATGAGAAGTGTGATATTATAACATTTCATGACTACTGCAAGACTAA','TCGCACGGGTGGATAAGCGTTTACAGTTTTCGCAAGCTCGTAAAAGCAGTACAGTGCACCGTAAGAAAATTACAAGTATAC']
getPCSF(sequences,'pcsfFeaFile.txt')
"""
