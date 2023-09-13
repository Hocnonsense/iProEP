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
    with (
        open(svmfile1, "r") as fobjr1,
        open(svmfile2, "r") as fobjr2,
        open(newsvmfile, "w") as fobjw,
    ):
        for eachline1, eachline2 in zip(fobjr1, fobjr2):
            eachline = eachline1.strip() + "\t" + eachline2.strip()
            fobjw.writelines(eachline + "\n")


def getOptimalFea(
    mrmrOrderFile: str, allFeaFile: str, feaNum: int, optimalFeaFile: str
):
    with open(mrmrOrderFile, "r") as orderFile:
        order = [eachline.strip().split("\t")[1] for eachline in orderFile]

    with open(allFeaFile, "r") as fobjr, open(optimalFeaFile, "w") as fobjw:
        for eachline in fobjr:
            eachvalue = eachline.strip().split("\t")
            fobjw.write(eachvalue[0] + "\t")
            for k in range(feaNum):
                index = int(order[k])
                content = "%d:%s\t" % (k + 1, eachvalue[index].split(":")[1])
                fobjw.write(content)
            fobjw.write("\n")


if __name__ == "__main__":
    # sequences = ['GCGCAAAAAATTTAACAATAGTCTAGCCATTTGAACTTTTTAATGCAAAATTAAACAAAAATTATTTGAAAATTTAGAAAAATTCTATAAATTATTTAACAGTTATATTCGGTCATTACGGCCTATTTTCTGTCCTTTTAGGCAAAGAGACCCCACAATAATGGAGCTTTTTGCCTTAAATGACAGGAAATAGGCCGTAATGACCGAATATAACTGTTAAAAAAATTATAGAAATTTTCTAAATTTTTTATGATTTTTTTTCAAAATTTTGATCACACGAAATTGAAGATCCCAAAAAAT','GAGCCGAAAGAATAGGGATCATGGCGCCAGTGAAGAGACCCAAAGAAGCACACAAATCAGAGCAGTCGGCTCAAGAAAAGAGAGGAATCCGGGTGCCAACATCATTATTGGGCCATTTTTGTCTGCGTGTGTGATCTAATTTGATCAATCTCTTACTCCCATCAGCCAGCTTGAACATTAAAAAACAAGTTGACCTTCTCCTCGATAGTGCCCTCATAAGGCGCTTAAACCCACCTTACCCTTACCATCATGGCTAGTCGACGCCAAAAGCAGTTCGATCGGAAGTACAGCTCCTATCGG','TTTCATAATATTGTGCTATTTGTTTTGTGACCGGGTTGTTCTAAAGATGATGTAGTATATAAATATTTACATTTGAAAACTTTGAGCTCAAAAATCAAAACACAAGTAAAACAAAAACAAGAAAAACATTGAAAGGAGATTATTCAAGGGTAAGGGAGGATCATAGAAACAGGGAAGCGAACCGGTGGAAGGGAAAGAAATTGCAACAAAAAACAGAAGAAAATGCATAAAGTTATGAATAAATTATATAGAAATTACAATGAGATCTCACTATTTGATGGAGATGGTGGAGATCGTCTA','TTCAATCTGTGATACTGGAAATTAACTATTCAACATTAGAATTATAGCGTGGGGACCGGTGCCTCAACTATTTAAGAATGTAGCTGCAAAATTCCAATATTGAAACTATTGGAATTTTGCAAAAAAAAAAACAAGTACAATACTTTTAAGTTAAATATTGCGGACGTGTTTTCTATTTTCACGGCCATCTCGCCAATGTCAACATTTTGCTTTCATTATTCGATTTTTTTTTCGAAATTTGACTAAAATTTCTGACTTGAAGCATCAAAGTTCAAAGCATACACGAAAAGAGCAACAAAT']
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
        "405ResultMRMR.txt",
        "pse&pcsfFea_" + inputFilename,
        67,
        "optimalFea_" + outputFilename,
    )
