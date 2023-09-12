<!--
 * @Date: 2023-09-10 17:42:30
 * @Editors: Hong-Yan Lai et al.
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-10 19:25:28
 * @FilePath: /iProEP_localtool/README.md
 * @Description:
-->
Introduction
===

- this ia modify version of iProEP edited by Hwrn

---

The iProEP suite provides online tools for predicting promoter of mRNA on Windows System.

You can download the iProEP local tool from
http://lin-group.cn/server/iProEP-mRNA/iProEP_localtool.rar

You can also use the iProEP via its website at
http://lin-group.cn/server/iProEP/


## Users Notes

To prepare a new release.

1. Download the iProEP local tool from
    http://lin-group.cn/server/iProEP/iProEP_localtool.rar

2. Prepare your query sequence in fasta format (ref example.txt)
    Each sequence must start with a greater-than symbol (">") in the first column.
    The words right after the ">" symbol in the single initial line are optional and only used for
    the purpose of identification and description.

3. Once you are ready to create a release candidate tag the current revision
    ```powershell
    # for windows
    cd iProEP_localtool/
    python3 iProEP.py -s C -i example.fasta
    ```
    ```bash
    # for Linux
    cd iProEP_localtool/
    cd libsvm-3.22 && make clean && make && cd ..
    python3 iProEP.py -s C -i example.fasta
    ```

4. Open All_Result.csv and Promoter.txt get the results.


## Citing

To cite the full MEME suite, please cite:
Hong-Yan Lai &, Zhao-Yue Zhang &, Zhen-Dong Su, Wei Su, Hui Ding, Wei Chen*, Hao Lin*.
(2019) iProEP: a computational predictor for predicting promoter. Molecular Therapy - Nucleic Acids, 17: 337-346. (2018 IF: 5.919)
