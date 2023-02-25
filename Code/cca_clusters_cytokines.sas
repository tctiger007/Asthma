title "Canonical Correlation Analysis";

proc import datafile="/home/wwang750/caa_morning/group3_rmoutlier.csv"
        out=group3
        dbms=csv
        replace;
     getnames=yes;
run;


/* cluster 1  */
proc cancorr out=group3_cca_c1 vprefix=cytokine vname="cytokine"
wprefix=group3_c1_ wname="Group3 Cluster 1 miRNA Variables";
		var PEFPercentPred FVCPercentPred FEV1PercentPred;
	with hsa_let_7f_5p
hsa_miR_26a_5p
hsa_let_7a_5p
hsa_miR_92a_3p
hsa_miR_21_5p
hsa_miR_126_3p
hsa_miR_30c_5p
hsa_miR_22_3p
hsa_miR_199a_3p
hsa_miR_221_5p
hsa_miR_103a_3p
hsa_miR_1_3p;
run;

/* cluster 2  */
proc cancorr out=group3_cca_c2 vprefix=cytokine vname="cytokine"
wprefix=group3_c2_ wname="Group3 Cluster 2 miRNA Variables";
		var PEFPercentPred FVCPercentPred FEV1PercentPred;
	with hsa_miR_122_5p
hsa_miR_148a_3p
hsa_miR_101_3p
hsa_miR_151a_3p
hsa_miR_375_3p
hsa_miR_10b_5p
hsa_miR_423_5p
hsa_miR_671_3p
hsa_miR_127_3p
hsa_miR_4433b_3p
hsa_miR_99a_5p;
run;

