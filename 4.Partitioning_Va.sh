###############################
## Va comparisons, delta=-0.25
###############################

#List of autosomal/X-linked loci (see R code)
cat list1 list2 > listALL
./ldak5.mac --cut-kins reml_XvsA_maf0.05 --bfile f3c.lhm.snp --partition-number 2 --partition-prefix autosome_X_partition/list
#Calculate kinships for each section
for i in {1..2}; do ./ldak5.mac --calc-kins reml_XvsA_maf0.05 --bfile f3c.lhm.snp --partition $i --weights sect_maf0.05/weights.short --power -0.25; done
#Va estimates using REML
./ldak5.mac --reml reml_XvsA_maf0.05/XvsA.male --mgrm reml_XvsA_maf0.05/partition.list --pheno pheno.txt --mpheno 1  --constrain REGIONS;
./ldak5.mac --reml reml_XvsA_maf0.05/XvsA.female --mgrm reml_XvsA_maf0.05/partition.list --pheno pheno.txt --mpheno 2 --constrain REGIONS;

#Permuted Va estimates using REML
cd reml_XvsA_perm_circ_maf0.05/
for j in {1..1000}; do ./ldak5.mac --cut-kins reml_XvsA_perm --bfile ../f3c.lhm.snp --partition-number 2 --partition-prefix list_${j}_; for i in {1..2}; do ../ldak5.mac --calc-kins reml_XvsA_perm --bfile ../f3c.lhm.snp --partition $i --weights ../sect_maf0.05/weights.short --power -0.25; ../ldak5.mac --reml reml_XvsA_perm/XvsA_male_${j} --mgrm reml_XvsA_perm/partition.list --pheno ../pheno.txt --mpheno 1 --constrain YES; ../ldak5.mac --reml reml_XvsA_perm/XvsA_female_${j} --mgrm reml_XvsA_perm/partition.list --pheno ../pheno.txt --mpheno 2 --constrain YES; done; done;

###############################
## Va comparisons, alternative alphas
###############################

#Delta = -1
for i in {1..2}; do ./ldak5.mac --calc-kins reml_XvsA_maf0.05 --bfile f3c.lhm.snp --partition $i --weights sect_maf0.05/weights.short --power -1; done
./ldak5.mac --reml reml_XvsA_maf0.05/XvsA.weight_minus_1.male --mgrm reml_XvsA_maf0.05/partition.list --pheno pheno.txt --mpheno 1 --constrain REGIONS;
./ldak5.mac --reml reml_XvsA_maf0.05/XvsA.weight_minus_1.female --mgrm reml_XvsA_maf0.05/partition.list --pheno pheno.txt --mpheno 2 --constrain REGIONS;
#Permuted
cd reml_XvsA_perm_circ_maf0.05/
for j in {1..1000}; do ../ldak5.mac --cut-kins reml_XvsA_perm_weight_minus_1 --bfile ../f3c.lhm.snp --partition-number 2 --partition-prefix list_${j}_; for i in {1..2}; do ../ldak5.mac --calc-kins reml_XvsA_perm_weight_minus_1 --bfile ../f3c.lhm.snp --partition $i --weights ../sect_maf0.05/weights.short --power -1; ../ldak5.mac --reml reml_XvsA_perm_weight_minus_1/XvsA_male_${j} --mgrm reml_XvsA_perm_weight_minus_1/partition.list --pheno ../pheno.txt --mpheno 1 --constrain YES; ../ldak5.mac --reml reml_XvsA_perm_weight_minus_1/XvsA_female_${j} --mgrm reml_XvsA_perm_weight_minus_1/partition.list --pheno ../pheno.txt --mpheno 2 --constrain YES; done; done;


#Delta = 0
for i in {1..2}; do ./ldak5.mac --calc-kins reml_XvsA_maf0.05 --bfile f3c.lhm.snp --partition $i --weights sect_maf0.05/weights.short --power 0; done
./ldak5.mac --reml reml_XvsA_maf0.05/XvsA.weight_0.male --mgrm reml_XvsA_maf0.05/partition.list --pheno pheno.txt --mpheno 1  --constrain REGIONS;
./ldak5.mac --reml reml_XvsA_maf0.05/XvsA.weight_0.female --mgrm reml_XvsA_maf0.05/partition.list --pheno pheno.txt --mpheno 2 --constrain REGIONS;
#Permuted
cd reml_XvsA_perm_circ_maf0.05/
for j in {1..1000}; do ../ldak5.mac --cut-kins reml_XvsA_perm_weight_0 --bfile ../f3c.lhm.snp --partition-number 2 --partition-prefix list_${j}_; for i in {1..2}; do ../ldak5.mac --calc-kins reml_XvsA_perm_weight_0 --bfile ../f3c.lhm.snp --partition $i --weights ../sect_maf0.05/weights.short --power 0; ../ldak5.mac --reml reml_XvsA_perm_weight_0/XvsA_male_${j} --mgrm reml_XvsA_perm_weight_0/partition.list --pheno ../pheno.txt --mpheno 1 --constrain YES; ../ldak5.mac --reml reml_XvsA_perm_weight_0/XvsA_female_${j} --mgrm reml_XvsA_perm_weight_0/partition.list --pheno ../pheno.txt --mpheno 2 --constrain YES; done; done;
