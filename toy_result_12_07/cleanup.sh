hadd -k -f /home/hep/palvare1/B2Kmumu_Kee/Sulley/toy_result_12_07/toystudyHistBremCatTsallisBkg_resultsMode7L0ETOSOnly_d.root /home/hep/palvare1/B2Kmumu_Kee/Sulley/toy_result_12_07/output_1/toystudyHistBremCatTsallisBkg_resultsMode7L0ETOSOnly_d.root
mkdir iodir
mv runOnBatch.sh.* iodir/.
/home/hep/palvare1/B2Kmumu_Kee/Sulley//bin/maintables /home/hep/palvare1/B2Kmumu_Kee/Sulley/toy_result_12_07/toystudyHistBremCatTsallisBkg_resultsMode7L0ETOSOnly_d.root params_L0ETOSOnly_d 199 92 440 16 /home/hep/palvare1/B2Kmumu_Kee/Sulley/toy_result_12_07//toyresults.dat 0
/home/hep/palvare1/B2Kmumu_Kee/Sulley//bin/maintables /home/hep/palvare1/B2Kmumu_Kee/Sulley/toy_result_12_07/toystudyHistBremCatTsallisBkg_resultsMode7L0ETOSOnly_d.root params_L0ETOSOnly_d 199 92 440 16 /home/hep/palvare1/B2Kmumu_Kee/Sulley/toy_result_12_07//toyresults.tex 1
pdflatex /home/hep/palvare1/B2Kmumu_Kee/Sulley/toy_result_12_07//toyresults.tex
