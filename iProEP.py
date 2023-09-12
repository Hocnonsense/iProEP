# coding=utf-8

import sys
import os
from typing import NamedTuple
import click
from pathlib import Path


class ModelPara(NamedTuple):
    pretypelength: int
    modelType: str
    bestn: int


def getmodelpara(preType):
    if preType == "H":
        return ModelPara(300, "hsmodel", 55)
    if preType == "D":
        return ModelPara(300, "dmmodel", 893)
    if preType == "C":
        return ModelPara(300, "cemodel", 67)
    if preType == "B":
        return ModelPara(81, "bsmodel", 82)
    if preType == "E":
        return ModelPara(81, "ecmodel", 410)
    raise KeyError("Unsupported Type Name")


def slide(seqfile: Path, pretypelength: int):
    annotation = ""
    seq = ""

    with open(seqfile, "r") as f, open("slide_seq.fasta", "w") as g:
        for eachline in f:
            if eachline[0] == ">":
                annotation = eachline.strip()
            elif eachline.strip() != "":
                seq = eachline.strip().upper()
                if len(seq) < pretypelength:
                    print(
                        "The length of the %s is less than %s-bp, it has no biological meaning for promoter."
                        % (annotation, pretypelength)
                    )
                    continue
                for i in range(len(seq) - pretypelength + 1):
                    g.write(
                        "%s,site:%s-%s\n%s\n"
                        % (
                            annotation,
                            i,
                            pretypelength + i,
                            seq[i : pretypelength + i],
                        )
                    )


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
    import pseKNC
    import PCSF

    annotation = []
    sequences = []

    with open("slide_seq.fasta", "r") as f:
        for eachline in f:
            templine = eachline.strip()
            if eachline[0] == ">":
                annotation.append(templine)
            else:
                sequences.append(eachline)
    if sequences == []:
        print("There is no sequence that meets the requirements. Program Terminated.")
        sys.exit()
    print("Feature Encoding(time-consuming)...\n...")
    pseKNC.pseKNC(sequences, "pseFea.txt")
    PCSF.getPCSF(sequences, "pcsfFea.txt")
    merge2svmfile("pse&pcsfFea.txt", "pseFea.txt", "pcsfFea.txt")
    getOptimalFea("ResultMRMR.txt", "pse&pcsfFea.txt", bestn, "optimalFea.txt")
    return annotation


def runSVM(pathPrefix: Path):
    commands = []
    commands.append(
        f"{pathPrefix}/libsvm-3.22/svm-scale -r scale.rule optimalFea.txt > seq.scale"
    )
    commands.append(
        f"{pathPrefix}/libsvm-3.22/svm-predict -b 1 seq.scale predict.model res.txt >out.txt"
    )

    for eachCmd in commands:
        os.system(eachCmd)


def generateResult(annotation: list, outfile: Path):
    outResult = "res.txt"
    countn = 0
    countp = 0
    with (
        open(outResult, "r") as f,
        open(outfile, "w") as g,
        open("../Promoter.txt", "w") as g1,
    ):
        g.write(
            "***All slide window sequence number prediction probability***\n"
            "Number\tannotation\tPromoter\tPromoter probability\n"
        )
        for eachline in f:
            temp = eachline.split(" ")
            if temp[0] == "labels":
                pass
            else:
                countn += 1
                temp = eachline.split(" ")
                if temp[0] == "1":
                    countp += 1
                    g.write(
                        "%s\t%s\tYes\t%s\n" % (countn, annotation[countn - 1], temp[1])
                    )
                    g1.write(
                        "%s is predicted as Promoter. Its probability score of being a promoter is %s.\n"
                        % (annotation[countn - 1], temp[1])
                    )
                if temp[0] == "2":
                    g.write(
                        "%s\t%s\tNo\t%s\n" % (countn, annotation[countn - 1], temp[1])
                    )
    return countp


@click.command(
    help="""
Useful:
    python PseKNC.py -s H -i inputFilename

  parameter description:

    -s: The species type for prediction. (required). "H" for H. sapiens, "D" for D. melanogaster, "C" for C. elegans, "B" for B. subtilis and "E" for E. coli, respectively.
    -i: The input filename. (required)
    """
)
@click.option("-s", "--species", type=str, help="species")
@click.option("-i", "--seqfile", type=str, help="input filename")
@click.option(
    "--debug",
    is_flag=True,
    flag_value=False,
    help="set this flag to keep immediate files",
)
def main(species: str, seqfile: str, debug: bool):
    pathPrefix = Path(__file__).parent
    mp = getmodelpara(species)
    Typepath = str(pathPrefix / mp.modelType)
    os.chdir(Typepath)

    if Typepath not in sys.path:
        sys.path.append(Typepath)

    slide(pathPrefix / seqfile, mp.pretypelength)
    print("Sliding Window Finished.")
    annotation = getFea(mp.bestn)

    runSVM(pathPrefix)
    countp = generateResult(annotation, pathPrefix / "All_Result.txt")

    if not debug:
        for each in [
            "pse&pcsfFea.txt",
            "pseFea.txt",
            "pcsfFea.txt",
            "optimalFea.txt",
            "seq.scale",
            "out.txt",
        ]:
            os.remove(each)
    print(
        f"Prediction Finished!\n"
        f"***There are {countp} promoter in the query sequences.***\n"
        f"Find All_Result.csv and Promoter.txt get the results."
    )


if __name__ == "__main__":
    main()
