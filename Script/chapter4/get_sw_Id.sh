
#grab the sw list
# ls ${inDir1}  > Id.txt
#
# #copy the list-files from raw data
# rsync -a ${inDir2} --files-from=Id.txt ${ouDir}
# # rsync -a ${inDir2} --files-from=${ouFile} ${ouDir}
#
# #remove .fa
# cat Id.txt | cut -f1 -d '.' > ${ouFile}
# rm Id.txt

#############################################
inD=("../test_human_mouse_rat/Data_3/Mafft/mapped_cds/"
    "../test_human_mouse_rat/Data_6/Mafft/mapped_cds/"
    "../test_human_mouse_rat/Data_9/Mafft/mapped_cds/"
    "../test_human_mouse_rat/Data_12/Mafft/mapped_cds/")

index=('3' '6' '9' '12')

for i in {0..3}
do
   ls ${inD[$i]} > hmr/Id.${i}.txt
   rsync -a hmr/Raw/cds/ --files-from=hmr/Id.${i}.txt hmr/Data/${index[$i]}/
 done

 rm hmr/Id.*.txt
