## CN=1
awk -F"," '{if($7==1){OFS="\t";print "chr"$2,$3,$4}}' colo829.copynumber.caveman.csv|grep -Ev "chrX|chrY" > colo829.cn1.tsv
awk -F"," '{if($7==1){OFS="\t";print "chr"$2,$3,$4}}' hcc1937.copynumber.caveman.csv|grep -Ev "chrX|chrY" > hcc1937.cn1.tsv

## CN=2
awk -F"," '{if($7==2){OFS="\t";print "chr"$2,$3,$4}}' colo829.copynumber.caveman.csv|grep -Ev "chrX|chrY" > colo829.cn2.tsv
awk -F"," '{if($7==2){OFS="\t";print "chr"$2,$3,$4}}' hcc1937.copynumber.caveman.csv|grep -Ev "chrX|chrY" > hcc1937.cn2.tsv

## CN=2 and het
awk -F"," '{if($7==2 && $8==1){OFS="\t";print "chr"$2,$3,$4}}' colo829.copynumber.caveman.csv|grep -Ev "chrX|chrY" > colo829.cn2het.tsv
awk -F"," '{if($7==2 && $8==1){OFS="\t";print "chr"$2,$3,$4}}' hcc1937.copynumber.caveman.csv|grep -Ev "chrX|chrY" > hcc1937.cn2het.tsv
