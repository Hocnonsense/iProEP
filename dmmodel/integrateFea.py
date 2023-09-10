# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 15:14:03 2018

@author: laihy
"""
import pseKNC
import PCSF
import os
import sys

def merge2svmfile(newsvmfile,svmfile1,svmfile2):
	fobjr1 = open('.\\5model\\dmmodel\\'+svmfile1,'r')
	fobjr2 = open('.\\5model\\dmmodel\\'+svmfile2,'r')
	fobjw = open('.\\5model\\dmmodel\\'+newsvmfile,'w')
	for eachline1 in fobjr1:
		eachline1= eachline1.strip()
		eachline2 = fobjr2.readline().strip()
		eachline = eachline1 +'\t'+ eachline2
		fobjw.writelines(eachline + '\n')
	fobjw.close()
	fobjr1.close()
	fobjr2.close()
    
def getOptimalFea(mrmrOrderFile,allFeaFile,feaNum,optimalFeaFile):
    orderFile = open('.\\5model\\dmmodel\\'+mrmrOrderFile,'r')
    order = []
    for eachline in orderFile:
            order.append(eachline.strip().split('\t')[1])       
    orderFile.close()
    
    fobjr = open('.\\5model\\dmmodel\\'+allFeaFile,'r')
    fobjw = open('.\\5model\\dmmodel\\'+optimalFeaFile,'w')
    for eachline in fobjr:
        eachline = eachline.strip().split('\t')
        fobjw.write(eachline[0] + '\t')
        for k in range(feaNum):
            index = int(order[k])
            content = "%d:%s\t" % (k+1,eachline[index].split(':')[1])
            fobjw.write(content)
        fobjw.write('\n')        
    fobjw.close()
    fobjr.close()


if __name__ == '__main__':
    #sequences = ['CGCGATCCCATCGGTTATCTGCAGAAGTTTTGGCTTGACGACAAGCACAAGTGTGACATTCCATTGCCGCAGAATATTCTTTAACATATAGCTCTGATGTAACCCAATCAAATTTCAAATTGTCAGCTAATTAAATGTTTTTGATTCGACGTTAAGCTTACGCAATGATGCATTTCTTTGATCAGTCAGGCGAAGAGATTCTGTACCGATTTAAGGGCTATAAAAGCAGATATTCCAATCTCAAACTGCATCATTCATTTCACAGTTTGCGGATCTAAAGGAAGACTGCTATAATCATGT','TCAAAAATTCATCATACCCTGCTTTCTCAAGGGTAAGGACATGCTCCTGATTCCGAGACCGCTTTCTCTTGGCTCTCTTTTGTCTCCCGCTCTCACTGTTTGTGCAACCACCTGCAGGCAGTCTGCCTTCCCATTGGCTGGGTCCATCGATCGAGAGAACCGCCGCTCTGTGCGAGAATTCTCCCCCGCCGAGAGAAACCAAAGAGAGCAAATCGTATATAAGGAGTCGGTCGTGCAGCAGGAAAACCAGAATACTACCAACTCTCAACCGAGTGCAACCAACTAAGATCCCAAGTTTTA','ATCTGGTGCAATTTGGGCTGCCAACTGTGGTCAAGTTAATCTATAAAATTCTTGATTCAAAATGCAATTTTAAGAACATTTTATTTCAAATTATATCAAATCAGATAAAAAACAAATGTCCATGAAATGATTAGGTTAGTTAAGGATATTTTCAATCAAAACACTATAAGAAATAGGACAAGGGAAGCCAAACATAAAATAAAGCCAAATTTAACTATAAACTAGCTTAAACGGCAGTGCAAACAGTTCATCTGGTCTGGCAACTCTCCGCCGAATTGTTGTTACGCGAGCGATTTCTTT','GCATTTTAAATAGCGCGCCAATTTTTTAAGCTGTTCACACTGCGCTTGCATGGAAGTAGGGCTTAAGTTATGGCAACTCGTTCTTTACTCAATGGTTTCAGTCCTTAGGTCTTATATAAGATATAAAAAGTATTTAGAAGTTATTATACATATTTTTATTAACATTTTTCTACATTCCTAGTCACATCGATAGCTAATCGCTGTCATCGATAGTCCCGTGATGTTTTTTGAGCTGCTTCCACTGCGGTCACACTAATCGAAGCATTTTCGTCCGCTCGTCACTTTTCTGCGTAATTCTTT']
    #读取文件并存列表
    inputFile = open(sys.argv[1])
    sequences = []

    for eachline in inputFile:
        eachline = eachline.strip('\n')
        if eachline[0] == '>':
                pass
        else:
               sequences.append(eachline)
               
    inputFile.close()


    ##指定输入输出文件名
    inputFilename = os.path.basename(sys.argv[1])
    outputFilename = os.path.basename(sys.argv[2])
    

    pseKNC.pseKNC(sequences, 'pseFea_'+inputFilename)
    PCSF.getPCSF(sequences, 'pcsfFea_'+inputFilename)
    merge2svmfile('pse&pcsfFea_'+inputFilename,'pseFea_'+inputFilename,'pcsfFea_'+inputFilename)
    getOptimalFea('1097ResultMRMR.txt','pse&pcsfFea_'+inputFilename,893,'optimalFea_'+outputFilename) 
