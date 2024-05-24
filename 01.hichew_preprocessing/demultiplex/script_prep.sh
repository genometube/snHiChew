proj_id="snHiChew_20240318"
lib_id_raw=($(ls -d /research/xieyeming1/proj_2023/${proj_id}/raw_data/*V350256554*))
cell_num="2000"   # number of cell to demultiplex, greater than the expected cell number

lib_id=(${lib_id_raw[@]##*/})

for i in ${lib_id[@]}; do
  lane_id_tmp=(${i[@]#*M027_})
  lane_id=(${lane_id_tmp[@]%_CWH*})

  mkdir -p ${i}

  cp template_script/*py template_script/sample* ${i}

  echo "lane_id=${lane_id}" > ${i}/demultiplex_i5_i7_test_run.sh
  echo "lib_id=${i}" >> ${i}/demultiplex_i5_i7_test_run.sh
  echo "proj_id=${proj_id}" >> ${i}/demultiplex_i5_i7_test_run.sh
  echo "cell_num=${cell_num}" >> ${i}/demultiplex_i5_i7_test_run.sh
  cat template_script/demultiplex_i5_i7_test_run.sh >> ${i}/demultiplex_i5_i7_test_run.sh

  echo "lane_id=${lane_id}" > ${i}/demultiplex_i5_i7.sh
  echo "lib_id=${i}" >> ${i}/demultiplex_i5_i7.sh
  echo "proj_id=${proj_id}" >> ${i}/demultiplex_i5_i7.sh
  echo "cell_num=${cell_num}" >> ${i}/demultiplex_i5_i7.sh
  cat template_script/demultiplex_i5_i7.sh >> ${i}/demultiplex_i5_i7.sh

done
