#grab the raw cds filterd from align_max alignments

inD="../test_90_species/Raw_data/align_max"
ouD="90/species"

for subD in ${inD}/*
do
   seed=$(basename ${subD})
   ls ${subD} > ${ouD}/${seed}.txt
done
