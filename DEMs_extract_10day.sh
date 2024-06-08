#!/bin/sh
cat 10day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >10d_XXXYD_orth_XYD_high_degs_list1.csv
cat 10day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >10d_XXXYD_orth_XX_high_degs_list1.csv
cat 10day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >10d_XXXY_orth_XY_high_degs_list1.csv
cat 10day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >10d_XXXY_orth_XX_high_degs_list1.csv
cat 10day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >10d_XYDXY_orth_XYD_high_degs_list1.csv
cat 10day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >10d_XYDXY_orth_XY_high_degs_list1.csv
