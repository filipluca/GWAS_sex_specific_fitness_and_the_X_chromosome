#######################
## LD pruning ###
#######################

#Prune all SNPs so that they are in approximate linkage equilibrium
~/Downloads/plink_mac/plink --bfile f3c.lhm.snp --clump female_maf0.05.assoc2 --clump-p1 1 --clump-p2 1 --clump-kb 10 --clump-r2 0.4 --clump-field Wald_P_randomised --out female_maf0.05_pruned_r0.4_10kb
#Prune coding SNPs so that they are in approximate linkage equilibrium
~/Downloads/plink_mac_20200428/plink --bfile f3c.lhm.snp --clump female_maf0.05.assoc3 --clump-p1 1 --clump-p2 1 --clump-kb 10 --clump-r2 0.4 --clump-field Wald_P_randomised --out female_maf0.05_coding_pruned_r0.4_10kb
#Prune candidates
~/Downloads/plink_mac_20190304/plink --noweb --bfile f3c.lhm.snp --clump male_maf0.05.assoc2 --clump-p1 0.00009665 --clump-p2 0.1 --clump-kb 10 --clump-r2 0.4 --clump-field Wald_P --out clumped_highly_associated_r0.4_10kb
