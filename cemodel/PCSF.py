# -*- coding: utf-8 -*-
"""
 * @Date: Wed Mar 29 15:51:56 2017
 * @author: laihongyan
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-12 22:38:12
 * @FilePath: /iProEP_localtool/cemodel/PCSF.py
 * @Description:
"""

import math
from collections import UserList
from functools import cache
from pathlib import Path
from typing import Final

from Bio import SeqIO


class Kmers(UserList):
    oligos: Final = ["A", "G", "C", "T"]

    def __init__(self, k):
        self.k = k
        self.p0 = 1 / 4**k
        if k == 1:
            self.data = self.oligos
        else:
            self.data = [i + j for j in type(self)(k - 1) for i in self.oligos]


class PromoterTrainData:
    def __init__(self, promoter: Path, seqlen: int, k: int, sites, num: int) -> None:
        self.promoter = promoter
        self.seqlen = seqlen
        self.kmers = Kmers(k)
        self.sites = sites
        self.num = num

    @cache
    def PWM(self):
        """构建整个位置权重矩阵（allsiteNum=82-k个位点）"""
        allsiteNum = self.seqlen - self.kmers.k + 1

        realCounts: list[dict[str, int]] = [
            {kmer: 0 for kmer in self.kmers} for _ in range(allsiteNum)
        ]

        ##############count the number of kmers at the corresponding allsites
        for seq in (str(i.seq) for i in SeqIO.parse(self.promoter, "fasta")):
            for i in range(len(seq) + 1 - self.kmers.k):
                realCounts[i][seq[i : i + self.kmers.k]] += 1

        ###直接对字典的所有values进行求和
        totalCounts = [sum(site.values()) for site in realCounts]
        # print("the total number of real counts at the %d site is %d" % (i,sum(realCounts[i].values())))
        ###math.sqrt()开平方
        psetotalCounts = [math.sqrt(site) for site in totalCounts]
        # print("psetotalCounts at the %d site is %d" % (i,psetotalCounts[i]))

        pwm = [
            {rk: (rv + p * self.kmers.p0) / (t + p) for rk, rv in r.items()}
            for r, p, t in zip(realCounts, psetotalCounts, totalCounts)
        ]
        # print('the size of pwm is %d' % (len(pwm)))
        return pwm

    def PCSF(self, sequences: list[str]):
        pwm = self.PWM()
        pcsf: list[list[float]] = []
        for eachline_ in sequences:
            if eachline_.startswith(">"):
                continue
            eachline = eachline_.strip()
            pcsfI = [
                math.log(pwm[i][eachline[i : i + self.kmers.k]] / self.kmers.p0, math.e)
                for i in range(len(eachline) + 1 - self.kmers.k)
                if i in self.sites
            ]
            pcsf.append(pcsfI)
        # print('the size of pcsf is %d' % (len(pcsf)))
        return pcsf

    def writePcsfFea(self, pcsf: list[list[float]], pcsfFeaFile: str):
        with open(pcsfFeaFile, "w") as fobjw:
            for i in range(len(pcsf)):
                for j in range(len(pcsf[i])):
                    temp = self.num + j
                    fobjw.write(f"{temp}:{pcsf[i][j]:f}\t")
                fobjw.write("\n")

    def getPCSF(self, sequences: list[str], pcsfFeaFile: str):
        self.writePcsfFea(self.PCSF(sequences), pcsfFeaFile)


ptd = PromoterTrainData(
    Path(__file__).parent / "promoter.txt",
    300,
    3,
    [
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
    ],
    389,
)


def getPCSF(sequences: list[str], pcsfFeaFile: str):
    ptd.getPCSF(sequences, pcsfFeaFile)


"""
sequences = ['AGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAATAAAGAGAGA','CGGGCCTATAAGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAA','TATGTAACATAATGCGACCAATAATCGTAATGAATATGAGAAGTGTGATATTATAACATTTCATGACTACTGCAAGACTAA','TCGCACGGGTGGATAAGCGTTTACAGTTTTCGCAAGCTCGTAAAAGCAGTACAGTGCACCGTAAGAAAATTACAAGTATAC']
getPCSF(sequences,'pcsfFeaFile.txt')
"""
