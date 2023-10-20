$MY_TEST_AREA/LArRecoND/bin/PandoraInterface -i $MY_TEST_AREA/LArRecoND/settings/development/PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml \
-j LArTPC -N -n 10000 -d ArgonCube -r AllHitsNu -e LArEDepSim_numu_all_1.root >& log.txt &

mv DUNEND_Accel_Pass1_CaloHitListU.csv DUNEND_Accel_Pass1_CaloHitListU_1.csv
mv DUNEND_Accel_Pass1_CaloHitListV.csv DUNEND_Accel_Pass1_CaloHitListV_1.csv
mv DUNEND_Accel_Pass1_CaloHitListW.csv DUNEND_Accel_Pass1_CaloHitListW_1.csv

$MY_TEST_AREA/LArRecoND/bin/PandoraInterface -i $MY_TEST_AREA/LArRecoND/settings/development/PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml \
-j LArTPC -N -n 10000 -d ArgonCube -r AllHitsNu -e LArEDepSim_numu_all_2.root >& log.txt &

mv DUNEND_Accel_Pass1_CaloHitListU.csv DUNEND_Accel_Pass1_CaloHitListU_2.csv
mv DUNEND_Accel_Pass1_CaloHitListV.csv DUNEND_Accel_Pass1_CaloHitListV_2.csv
mv DUNEND_Accel_Pass1_CaloHitListW.csv DUNEND_Accel_Pass1_CaloHitListW_2.csv

$MY_TEST_AREA/LArRecoND/bin/PandoraInterface -i $MY_TEST_AREA/LArRecoND/settings/development/PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml \
-j LArTPC -N -n 10000 -d ArgonCube -r AllHitsNu -e LArEDepSim_numu_all_3.root >& log.txt &

mv DUNEND_Accel_Pass1_CaloHitListU.csv DUNEND_Accel_Pass1_CaloHitListU_3.csv
mv DUNEND_Accel_Pass1_CaloHitListV.csv DUNEND_Accel_Pass1_CaloHitListV_3.csv
mv DUNEND_Accel_Pass1_CaloHitListW.csv DUNEND_Accel_Pass1_CaloHitListW_3.csv

$MY_TEST_AREA/LArRecoND/bin/PandoraInterface -i $MY_TEST_AREA/LArRecoND/settings/development/PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml \
-j LArTPC -N -n 10000 -d ArgonCube -r AllHitsNu -e LArEDepSim_numu_all_4.root >& log.txt &

mv DUNEND_Accel_Pass1_CaloHitListU.csv DUNEND_Accel_Pass1_CaloHitListU_4.csv
mv DUNEND_Accel_Pass1_CaloHitListV.csv DUNEND_Accel_Pass1_CaloHitListV_4.csv
mv DUNEND_Accel_Pass1_CaloHitListW.csv DUNEND_Accel_Pass1_CaloHitListW_4.csv

$MY_TEST_AREA/LArRecoND/bin/PandoraInterface -i $MY_TEST_AREA/LArRecoND/settings/development/PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml \
-j LArTPC -N -n 10000 -d ArgonCube -r AllHitsNu -e LArEDepSim_numu_all_5.root >& logt.txt &

mv DUNEND_Accel_Pass1_CaloHitListU.csv DUNEND_Accel_Pass1_CaloHitListU_5.csv
mv DUNEND_Accel_Pass1_CaloHitListV.csv DUNEND_Accel_Pass1_CaloHitListV_5.csv
mv DUNEND_Accel_Pass1_CaloHitListW.csv DUNEND_Accel_Pass1_CaloHitListW_5.csv

cat DUNEND_Accel_Pass1_CaloHitListU_1.csv DUNEND_Accel_Pass1_CaloHitListU_2.csv DUNEND_Accel_Pass1_CaloHitListU_3.csv \
DUNEND_Accel_Pass1_CaloHitListU_4.csv DUNEND_Accel_Pass1_CaloHitListU_5.csv > DUNEND_Accel_Pass1_CaloHitListU.csv

cat DUNEND_Accel_Pass1_CaloHitListV_1.csv DUNEND_Accel_Pass1_CaloHitListV_2.csv DUNEND_Accel_Pass1_CaloHitListV_3.csv \
DUNEND_Accel_Pass1_CaloHitListV_4.csv DUNEND_Accel_Pass1_CaloHitListV_5.csv > DUNEND_Accel_Pass1_CaloHitListV.csv

cat DUNEND_Accel_Pass1_CaloHitListW_1.csv DUNEND_Accel_Pass1_CaloHitListW_2.csv DUNEND_Accel_Pass1_CaloHitListW_3.csv \
DUNEND_Accel_Pass1_CaloHitListW_4.csv DUNEND_Accel_Pass1_CaloHitListW_5.csv > DUNEND_Accel_Pass1_CaloHitListW.csv

mkdir -p csv
mv *.csv csv/.
