# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-10 17:42:30
 * @Editors: Hong-Yan Lai et al.
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-12 17:11:13
 * @FilePath: /iProEP_localtool/iProEP.py
 * @Description:
"""
# coding=utf-8

import os
import sys
from pathlib import Path
from typing import NamedTuple

import click
from Bio import SeqIO, SeqRecord


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
    for seq in (i.upper() for i in SeqIO.parse(seqfile, "fasta")):
        if len(seq) < pretypelength:
            print(
                f"The length of the {seq.description} is less than {pretypelength}-bp, "
                f"it has no biological meaning for promoter."
            )
            continue
        for i in range(len(seq) - pretypelength + 1):
            seq1 = SeqRecord.SeqRecord(
                seq.seq[i : pretypelength + i],
                id=seq.id,
                description=seq.description.split(maxsplit=1)[1]
                + f",site:{i}-{pretypelength + i}",
            )
            yield seq1


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


def getFea(bestn, silce_fa):
    import PCSF
    import pseKNC

    annotation = []
    sequences = []

    for seq in SeqIO.parse(silce_fa, "fasta"):
        annotation.append(seq.description)
        sequences.append(str(seq.seq))

    assert (
        sequences
    ), "There is no sequence that meets the requirements. Program Terminated."
    print("Feature Encoding(time-consuming)...\n...")

    pseKNC.pseKNC(sequences, "pseFea.txt")
    PCSF.getPCSF(sequences, "pcsfFea.txt")
    merge2svmfile("pse&pcsfFea.txt", "pseFea.txt", "pcsfFea.txt")
    getOptimalFea("ResultMRMR.txt", "pse&pcsfFea.txt", bestn, "optimalFea.txt")
    return annotation


def runSVM(pathPrefix: Path, Typepath: str):
    os.system(
        f"{pathPrefix}/libsvm-3.22/svm-scale "
        f"-r {Typepath}/scale.rule optimalFea.txt > seq.scale"
    )
    os.system(
        f"{pathPrefix}/libsvm-3.22/svm-predict "
        f"-b 1 seq.scale {Typepath}/predict.model res.txt > out.txt"
    )


def generateResult(annotation: list, outfile: Path, outseq: Path):
    outResult = "res.txt"
    countn = 0
    countp = 0
    with open(outResult, "r") as f, open(outfile, "w") as g, open(outseq, "w") as g1:
        g.write("Number\tannotation\tPromoter\tPromoter probability\n")
        for eachline in f:
            temp = eachline.split(" ")
            if temp[0] == "labels":
                continue
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
    flag_value=True,
    help="set this flag to keep immediate files",
)
def main(species: str, seqfile: str, debug: bool):
    pathPrefix = Path(__file__).parent
    mp = getmodelpara(species)
    Typepath = str(pathPrefix / mp.modelType)
    os.chdir(Typepath)

    if Typepath not in sys.path:
        sys.path.append(Typepath)

    SeqIO.write(
        slide(pathPrefix / seqfile, mp.pretypelength), "slide_seq.fasta", "fasta-2line"
    )
    print("Sliding Window Finished.")
    annotation = getFea(mp.bestn, "slide_seq.fasta")

    runSVM(pathPrefix, Typepath)
    countp = generateResult(
        annotation, pathPrefix / "All_Result.txt", pathPrefix / "Promoter.txt"
    )

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
