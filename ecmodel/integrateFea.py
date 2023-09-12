# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 15:14:03 2018

@author: laihy
"""
import pseKNC
import PCSF
import os
import sys


def merge2svmfile(newsvmfile, svmfile1, svmfile2):
    fobjr1 = open(svmfile1, "r")
    fobjr2 = open(svmfile2, "r")
    fobjw = open(newsvmfile, "w")
    for eachline1 in fobjr1:
        eachline1 = eachline1.strip()
        eachline2 = fobjr2.readline().strip()
        eachline = eachline1 + "\t" + eachline2
        fobjw.writelines(eachline + "\n")
    fobjw.close()
    fobjr1.close()
    fobjr2.close()


def getOptimalFea(mrmrOrderFile, allFeaFile, feaNum, optimalFeaFile):
    with open(mrmrOrderFile, "r") as orderFile:
        order = [eachline.strip().split("\t")[1] for eachline in orderFile]

    fobjr = open(allFeaFile, "r")
    fobjw = open(optimalFeaFile, "w")
    for eachline in fobjr:
        eachvalue = eachline.strip().split("\t")
        fobjw.write(eachvalue[0] + "\t")
        for k in range(feaNum):
            index = int(order[k])
            content = "%d:%s\t" % (k + 1, eachvalue[index].split(":")[1])
            fobjw.write(content)
        fobjw.write("\n")
    fobjw.close()
    fobjr.close()


if __name__ == "__main__":
    # sequences = ['AGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAATAAAGAGAGA','CGGGCCTATAAGCCAGGCGAGATATGATCTATATCAATTTCTCATCTATAATGCTTTGTTAGTATCTCGTCGCCGACTTAA','TATGTAACATAATGCGACCAATAATCGTAATGAATATGAGAAGTGTGATATTATAACATTTCATGACTACTGCAAGACTAA','TCGCACGGGTGGATAAGCGTTTACAGTTTTCGCAAGCTCGTAAAAGCAGTACAGTGCACCGTAAGAAAATTACAAGTATAC']
    # 读取文件并存列表
    inputFile = open(sys.argv[1])
    sequences = []

    for eachline in inputFile:
        eachline = eachline.strip("\n")
        if eachline[0] == ">":
            pass
        else:
            sequences.append(eachline)

    inputFile.close()

    ##指定输入输出文件名
    inputFilename = os.path.basename(sys.argv[1])
    outputFilename = os.path.basename(sys.argv[2])

    pseKNC.pseKNC(sequences, "pseFea_" + inputFilename)
    PCSF.getPCSF(sequences, "pcsfFea_" + inputFilename)
    merge2svmfile(
        "pse&pcsfFea_" + inputFilename,
        "pseFea_" + inputFilename,
        "pcsfFea_" + inputFilename,
    )
    getOptimalFea(
        "345ResultMRMR.txt",
        "pse&pcsfFea_" + inputFilename,
        82,
        "optimalFea_" + outputFilename,
    )
