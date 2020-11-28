#######################
## Prepare kinships ###
#######################

#Cut genome into sections (default settings)
./ldak5.mac --cut-weights sect_maf0.05 --bfile f3c.lhm.snp
#Calculate weightings
./ldak5.mac --calc-weights-all sect_maf0.05 --bfile f3c.lhm.snp
#Calculate kinships
./ldak5.mac --calc-kins-direct kinsm_maf0.05 --bfile f3c.lhm.snp --weights sect_maf0.05/weights.all --power -.25 --kinship-raw YES

#######################
## GWAS ###
#######################

## GWAS of male and female fitness
./ldak5.mac --linear male_maf0.05 --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_maf0.05 --mpheno 1
./ldak5.mac --linear female_maf0.05 --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_maf0.05 --mpheno 2

#######################
## Permuted phenotypes ###
#######################

#Residual phenotypic values once kinships have been regressed out
./ldak5.mac --reml male_maf0.05_reml --grm kinsm_maf0.05 --pheno pheno.txt --mpheno 1
./ldak5.mac --reml female_maf0.05_reml --grm kinsm_maf0.05 --pheno pheno.txt --mpheno 2

#See R code for how to prepare permuted phenos

#Linear association using plink (no kinship matrix, as it has been regressed out already) for 1000 permuted phenotypes (based on '(fe)male_maf0.05_reml.indi.res')
for i in {1..1000}; do ~/Downloads/plink_mac/plink --assoc --pheno 1000_permuted_phenos/pheno_${i}_.txt --bfile f3c.lhm.snp --mpheno 1 --out 1000_permuted_phenos/male_maf0.05_permuted_${i} --allow-no-sex; awk '{ print $5 }' 1000_permuted_phenos/male_maf0.05_permuted_${i}.qassoc > 1000_permuted_phenos/male_maf0.05_permuted_${i}.effect; rm 1000_permuted_phenos/male_maf0.05_permuted_${i}.qassoc; rm 1000_permuted_phenos/male_maf0.05_permuted_${i}.nosex; rm 1000_permuted_phenos/male_maf0.05_permuted_${i}.log; done
for i in {1..1000}; do ~/Downloads/plink_mac/plink --assoc --pheno 1000_permuted_phenos/pheno_${i}_.txt --bfile f3c.lhm.snp --mpheno 2 --out 1000_permuted_phenos/female_maf0.05_permuted_${i} --allow-no-sex; awk '{ print $5 }' 1000_permuted_phenos/female_maf0.05_permuted_${i}.qassoc > 1000_permuted_phenos/female_maf0.05_permuted_${i}.effect; rm 1000_permuted_phenos/female_maf0.05_permuted_${i}.qassoc; rm 1000_permuted_phenos/female_maf0.05_permuted_${i}.nosex; rm 1000_permuted_phenos/female_maf0.05_permuted_${i}.log; done

#######################
## Repeat all using different values of alpha ###
#######################

## Alpha = -1
#Calculate kinships
./ldak5.mac --calc-kins-direct kinsm_maf0.05_weight_minus_1 --bfile f3c.lhm.snp --weights sect_maf0.05/weights.all --power -1 --kinship-raw YES
#GWAS
./ldak5.mac --linear male_maf0.05_weight_minus_1 --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_maf0.05_weight_minus_1 --mpheno 1;
./ldak5.mac --linear female_maf0.05_weight_minus_1 --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_maf0.05_weight_minus_1 --mpheno 2;
#Residual phenotypic values once kinships have been regressed out
./ldak5.mac --reml male_maf0.05_weight_minus_1_reml --grm kinsm_maf0.05_weight_minus_1 --pheno pheno.txt --mpheno 1;
./ldak5.mac --reml female_maf0.05_weight_minus_1_reml --grm kinsm_maf0.05_weight_minus_1 --pheno pheno.txt --mpheno 2;
#GWAS on permuted phenotypic values
for i in {1..1000}; do ~/Downloads/plink_mac_20200428/plink --assoc --pheno 1000_permuted_phenos_weight_minus_1/pheno_${i}_.txt --bfile f3c.lhm.snp --mpheno 1 --out 1000_permuted_phenos_weight_minus_1/male_maf0.05_permuted_${i} --allow-no-sex; awk '{ print $5 }' 1000_permuted_phenos_weight_minus_1/male_maf0.05_permuted_${i}.qassoc > 1000_permuted_phenos_weight_minus_1/male_maf0.05_permuted_${i}.effect; rm 1000_permuted_phenos_weight_minus_1/male_maf0.05_permuted_${i}.qassoc; rm 1000_permuted_phenos_weight_minus_1/male_maf0.05_permuted_${i}.nosex; rm 1000_permuted_phenos_weight_minus_1/male_maf0.05_permuted_${i}.log; done
for i in {1..1000}; do ~/Downloads/plink_mac_20200428/plink --assoc --pheno 1000_permuted_phenos_weight_minus_1/pheno_${i}_.txt --bfile f3c.lhm.snp --mpheno 2 --out 1000_permuted_phenos_weight_minus_1/female_maf0.05_permuted_${i} --allow-no-sex; awk '{ print $5 }' 1000_permuted_phenos_weight_minus_1/female_maf0.05_permuted_${i}.qassoc > 1000_permuted_phenos_weight_minus_1/female_maf0.05_permuted_${i}.effect; rm 1000_permuted_phenos_weight_minus_1/female_maf0.05_permuted_${i}.qassoc; rm 1000_permuted_phenos_weight_minus_1/female_maf0.05_permuted_${i}.nosex; rm 1000_permuted_phenos_weight_minus_1/female_maf0.05_permuted_${i}.log; done


## Alpha = 0
#Calculate kinships
./ldak5.mac --calc-kins-direct kinsm_maf0.05_weight_0 --bfile f3c.lhm.snp --weights sect_maf0.05/weights.all --power 0 --kinship-raw YES
#GWAS
./ldak5.mac --linear male_maf0.05_weight_0 --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_maf0.05_weight_0 --mpheno 1;
./ldak5.mac --linear female_maf0.05_weight_0 --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_maf0.05_weight_0 --mpheno 2;
#Residual phenotypic values once kinships have been regressed out
./ldak5.mac --reml male_maf0.05_weight_0_reml --grm kinsm_maf0.05_weight_0 --pheno pheno.txt --mpheno 1;
./ldak5.mac --reml female_maf0.05_weight_0_reml --grm kinsm_maf0.05_weight_0 --pheno pheno.txt --mpheno 2;
#Permuted phenotypic values
for i in {1..1000}; do ~/Downloads/plink_mac_20200428/plink --assoc --pheno 1000_permuted_phenos_weight_0/pheno_${i}_.txt --bfile f3c.lhm.snp --mpheno 1 --out 1000_permuted_phenos_weight_0/male_maf0.05_permuted_${i} --allow-no-sex; awk '{ print $5 }' 1000_permuted_phenos_weight_0/male_maf0.05_permuted_${i}.qassoc > 1000_permuted_phenos_weight_0/male_maf0.05_permuted_${i}.effect; rm 1000_permuted_phenos_weight_0/male_maf0.05_permuted_${i}.qassoc; rm 1000_permuted_phenos_weight_0/male_maf0.05_permuted_${i}.nosex; rm 1000_permuted_phenos_weight_0/male_maf0.05_permuted_${i}.log; done
for i in {1..1000}; do ~/Downloads/plink_mac_20200428/plink --assoc --pheno 1000_permuted_phenos_weight_0/pheno_${i}_.txt --bfile f3c.lhm.snp --mpheno 2 --out 1000_permuted_phenos_weight_0/female_maf0.05_permuted_${i} --allow-no-sex; awk '{ print $5 }' 1000_permuted_phenos_weight_0/female_maf0.05_permuted_${i}.qassoc > 1000_permuted_phenos_weight_0/female_maf0.05_permuted_${i}.effect; rm 1000_permuted_phenos_weight_0/female_maf0.05_permuted_${i}.qassoc; rm 1000_permuted_phenos_weight_0/female_maf0.05_permuted_${i}.nosex; rm 1000_permuted_phenos_weight_0/female_maf0.05_permuted_${i}.log; done


#######################
## Gene-based assosiation analysis ###
#######################

## Male
./ldak5.mac --cut-genes gbat_male_maf0.05 --bfile f3c.lhm.snp --genefile ensgenes.txt --ignore-weights YES --gene-buffer 5000 --overlap NO; 
./ldak5.mac --calc-genes-reml gbat_male_maf0.05 --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 1 --covar kinsm_maf0.05.vect --num-perms 1000;
## Female
./ldak5.mac --cut-genes gbat_female_maf0.05 --bfile f3c.lhm.snp --genefile ensgenes.txt --ignore-weights YES --gene-buffer 5000 --overlap NO;
./ldak5.mac --calc-genes-reml gbat_female_maf0.05 --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 2 --covar kinsm_maf0.05.vect --num-perms 1000;


#######################
## Band-based assosiation analysis ###
#######################

## Male
./ldak5.mac --cut-genes band_bat_male_maf0.05 --bfile f3c.lhm.snp --genefile Dmel_chr_bands.txt --power -.25 --ignore-weights YES; 
./ldak5.mac --calc-genes-reml band_bat_male_maf0.05 --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 1 --covar kinsm_maf0.05.vect --num-perms 1000;
## Female
./ldak5.mac --cut-genes band_bat_female_maf0.05 --bfile f3c.lhm.snp --genefile Dmel_chr_bands.txt --power -.25 --ignore-weights YES;
./ldak5.mac --calc-genes-reml band_bat_female_maf0.05 --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 2 --covar kinsm_maf0.05.vect --num-perms 1000;


#######################3
## Random chunk-based assosiation analysis ###
#######################

## Male
#for i in {1..1000}; do ./ldak5.mac --cut-genes chunk_bat_male_maf0.05 --bfile f3c.lhm.snp --genefile Dmel_chr_chunks/Dmel_chr_chunks${i}.txt --power -.25 --ignore-weights YES; ./ldak5.mac --calc-genes-reml chunk_bat_male_maf0.05 --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 1 --covar kinsm_maf0.05.vect; cut -d " " -f2,3,4 chunk_bat_male_maf0.05/remls.1 > chunk_bat_male_maf0.05/remls${i}; done
mkdir chunk_bat_male_maf0.05
parallel --jobs 3 './ldak5.mac --cut-genes chunk_bat_male_maf0.05{} --bfile f3c.lhm.snp --genefile Dmel_chr_chunks/Dmel_chr_chunks{}.txt --power -.25 --ignore-weights YES; ./ldak5.mac --calc-genes-reml chunk_bat_male_maf0.05{} --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 1 --covar kinsm_maf0.05.vect; cut -d " " -f2,3,4 chunk_bat_male_maf0.05{}/remls.1 > chunk_bat_male_maf0.05/remls{}; rm -r chunk_bat_male_maf0.05{}' ::: {1..1000}

## Female
mkdir chunk_bat_female_maf0.05
#for i in {1..1000}; do ./ldak5.mac --cut-genes chunk_bat_female_maf0.05 --bfile f3c.lhm.snp --genefile Dmel_chr_chunks/Dmel_chr_chunks${i}.txt --power -.25 --ignore-weights YES; ./ldak5.mac --calc-genes-reml chunk_bat_female_maf0.05 --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 2 --covar kinsm_maf0.05.vect; cut -d " " -f2,3,4 chunk_bat_female_maf0.05/remls.1 > chunk_bat_female_maf0.05/remls${i}; done
parallel --jobs 3 './ldak5.mac --cut-genes chunk_bat_female_maf0.05{} --bfile f3c.lhm.snp --genefile Dmel_chr_chunks/Dmel_chr_chunks{}.txt --power -.25 --ignore-weights YES; ./ldak5.mac --calc-genes-reml chunk_bat_female_maf0.05{} --bfile f3c.lhm.snp --ignore-weights YES --power -.25 --pheno pheno.txt --mpheno 2 --covar kinsm_maf0.05.vect; cut -d " " -f2,3,4 chunk_bat_female_maf0.05{}/remls.1 > chunk_bat_female_maf0.05/remls{}; rm -r chunk_bat_female_maf0.05{}' ::: {1..1000}
