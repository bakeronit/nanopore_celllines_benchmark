module load python/3.9.13

for i in `ls  analysis/svs/nanomonsv/COLO829*/*.nanomonsv.result.filt.pass.svtype.txt`; do \
    echo -en "${i}\t";python ~/working/share/sv_comp/svcomp_cli.py dcc \
    -q base/COLO829.somatic.cat12.dcc \
    -s $i --notype|grep Recall
done

for i in `ls  analysis/svs/nanomonsv/HCC1937*/*.nanomonsv.result.filt.pass.svtype.txt`; do \
    echo -en "${i}\t";python ~/working/share/sv_comp/svcomp_cli.py dcc \
    -q base/HCC1937.somatic.cat12.dcc \
    -s $i --notype|grep Recall
done
