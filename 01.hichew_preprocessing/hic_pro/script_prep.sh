proj_id="snHiChew_20240318"
lib_id_raw=($(ls -d /research/xieyeming1/proj_2023/${proj_id}/raw_data/*V350256554*))
cell_num="2000"   # number of cell to demultiplex, greater than the expected cell number

lib_id=(${lib_id_raw[@]##*/})

for i in ${lib_id[@]}; do
  lane_id_tmp=(${i[@]#*M027_})
  lane_id=(${lane_id_tmp[@]%_CWH*})

  mkdir -p ${i}

  cp template_script/*txt ${i}

  echo "lane_id=${lane_id}" > ${i}/hic-pro-batch-generator.sh
  echo "lib_id=${i}" >> ${i}/hic-pro-batch-generator.sh
  echo "proj_id=${proj_id}" >> ${i}/hic-pro-batch-generator.sh
  echo "cell_num=${cell_num}" >> ${i}/hic-pro-batch-generator.sh
  cat template_script/hic-pro-batch-generator.sh >> ${i}/hic-pro-batch-generator.sh

  echo "lane_id=${lane_id}" > ${i}/hic-pro-single.sh
  echo "lib_id=${i}" >> ${i}/hic-pro-single.sh
  echo "proj_id=${proj_id}" >> ${i}/hic-pro-single.sh
  echo "cell_num=${cell_num}" >> ${i}/hic-pro-single.sh
  cat template_script/hic-pro-single.sh >> ${i}/hic-pro-single.sh

done
