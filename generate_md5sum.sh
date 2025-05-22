outfile=md5sum_table.txt
echo "          MD5 Checksum           | Filename                                        " > $outfile
echo "---------------------------------|-------------------------------------------------" >> $outfile

for file in `cat data_essential.list data_big.list`; do
    openssl md5 -r $file >> $outfile
done
