# coding=utf-8

import sys
import os
import shutil


def printHelpInfo():
    print(
        """
Useful:
    python PseKNC.py -s H -i inputFilename

  parameter description:

    -s: The species type for prediction. (required). "H" for H. sapiens, "D" for D. melanogaster, "C" for C. elegans, "B" for B. subtilis and "E" for E. coli, respectively.
    -i: The input filename. (required)
    """
    )
    exit(0)


def obtainExternalParameters():
    try:
        if sys.argv[1] == "-s":
            species = sys.argv[2]
        else:
            assert 0
        if sys.argv[3] == "-i":
            in_filename = sys.argv[4]
        else:
            assert 0
    except:
        printHelpInfo()
    return species, in_filename


def getmodelpara(preType):
    if preType == "H":
        pretypelength = 300
        modelType = "hsmodel"
        bestn = 55
    if preType == "D":
        pretypelength = 300
        modelType = "dmmodel"
        bestn = 893
    if preType == "C":
        pretypelength = 300
        modelType = "cemodel"
        bestn = 67
    if preType == "B":
        pretypelength = 81
        modelType = "bsmodel"
        bestn = 82
    if preType == "E":
        pretypelength = 81
        modelType = "ecmodel"
        bestn = 410
    return pretypelength, modelType, bestn


def slide(seqfile, pretypelength):
    f = open(seqfile, "r")
    g = open("slide_seq.fasta", "w")
    annotation = ""

    seq = ""
    for eachline in f:
        if eachline[0] == ">":
            annotation = eachline.strip()
        elif eachline.strip() != "":
            seq = eachline.strip().upper()
            if len(seq) >= pretypelength:
                for i in range(len(seq) - pretypelength + 1):
                    g.write(
                        "%s,site:%s-%s\n%s\n"
                        % (annotation, i, pretypelength + i, seq[i : pretypelength + i])
                    )
            else:
                print(
                    "The length of the %s is less than %s-bp, it has no biological meaning for promoter."
                    % (annotation, pretypelength)
                )
    g.close()
    f.close()


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
    orderFile = open(mrmrOrderFile, "r")
    order = []
    for eachline in orderFile:
        order.append(eachline.strip().split("\t")[1])
    orderFile.close()

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


def getFea(bestn):
    f = open("slide_seq.fasta", "r")
    annotation = []
    sequences = []

    for eachline in f:
        templine = eachline.strip()
        if eachline[0] == ">":
            annotation.append(templine)
        else:
            sequences.append(eachline)
    f.close()
    if sequences == []:
        print("There is no sequence that meets the requirements. Program Terminated.")
        sys.exit()
    print("Feature Encoding(time-consuming)...\n...")
    pseKNC.pseKNC(sequences, "pseFea.txt")
    PCSF.getPCSF(sequences, "pcsfFea.txt")
    merge2svmfile("pse&pcsfFea.txt", "pseFea.txt", "pcsfFea.txt")
    getOptimalFea("ResultMRMR.txt", "pse&pcsfFea.txt", bestn, "optimalFea.txt")
    return annotation


def runSVM(pathPrefix):
    commands = []
    if "linux" in sys.platform:
        commands.append(
            r"%slibsvm-3.22/svm-scale -r scale.rule optimalFea.txt > seq.scale"
            % pathPrefix
        )
        commands.append(
            r"%slibsvm-3.22/svm-predict -b 1 seq.scale predict.model res.txt >out.txt"
            % pathPrefix
        )
    if "win" in sys.platform:
        commands.append(
            r"%ssvm-scale.exe -r scale.rule optimalFea.txt > seq.scale" % pathPrefix
        )
        commands.append(
            r"%ssvm-predict.exe -b 1 seq.scale predict.model res.txt >out.txt"
            % pathPrefix
        )

    for eachCmd in commands:
        os.system(eachCmd)


def generateResult(annotation, outfile):
    outResult = "res.txt"
    f = open(outResult, "r")
    countn = 0
    g = open(outfile, "w")
    g.write("***All slide window sequence number prediction probability***\n")
    g.write("Number\tannotation\tPromoter\tPromoter probability\n")
    g1 = open("../Promoter.txt", "w")
    countp = 0
    for eachline in f:
        temp = eachline.split(" ")
        if temp[0] == "labels":
            pass
        else:
            countn += 1
            temp = eachline.split(" ")
            if temp[0] == "1":
                countp += 1
                g.write("%s\t%s\tYes\t%s\n" % (countn, annotation[countn - 1], temp[1]))
                g1.write(
                    "%s is predicted as Promoter. Its probability score of being a promoter is %s.\n"
                    % (annotation[countn - 1], temp[1])
                )
            if temp[0] == "2":
                g.write("%s\t%s\tNo\t%s\n" % (countn, annotation[countn - 1], temp[1]))
    g.close()
    f.close()
    g1.close()
    return countp


if __name__ == "__main__":
    pathPrefix = os.getcwd() + os.path.sep
    [species, seqfile] = obtainExternalParameters()
    [pretypelength, modelType, bestn] = getmodelpara(species)
    Typepath = pathPrefix + modelType + os.path.sep
    os.chdir(Typepath)

    if not Typepath in sys.path:
        sys.path.append(Typepath)
    import pseKNC
    import PCSF

    slide(pathPrefix + seqfile, pretypelength)
    print("Sliding Window Finished.")
    annotation = getFea(bestn)

    runSVM(pathPrefix)
    countp = generateResult(annotation, pathPrefix + "All_Result.txt")
    filelist = [
        "pse&pcsfFea.txt",
        "pseFea.txt",
        "pcsfFea.txt",
        "optimalFea.txt",
        "seq.scale",
        "out.txt",
    ]
    for each in filelist:
        os.remove(each)
    print(
        "Prediction Finished!\n***There are %s promoter in the query sequences.***\nFind All_Result.csv and Promoter.txt get the results."
        % countp
    )
