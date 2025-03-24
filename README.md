Run each script one by one.

Input for script 2.1.2 should contain these mutation types: frameshift insertion, frameshift deletion, non frameshift insertion, non frameshift deletion, and stoploss. Mutation types of startloss and stopgain should not be included in the input as they won't generate no novel epitopes. The script has a mechanism to filter out same epitopes from WT protein. It's just haven't been tested for startloss and stopgain mutations.

The NetMHC output ".xls" file is not a binary XLS file but a tab-separated text file. So before running script 4, simply rename the NetMHC output ".xls" file as ".txt": e.g. `mv Pt19R2_missenceSNVs_NetMHCpan_out.xls Pt19R2_missenceSNVs_NetMHCpan_out.txt`.

Then run script 4 by: `python 4.Format_NetMHC_output.py [NetMHCpan_output].txt [epitope_generation_output].fasta [formatted_output].xlsx`
