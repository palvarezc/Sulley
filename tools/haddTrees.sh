#!/bin/bash

# You need SetupProject Urania v3r0 first!!

#path="/vols/lhcbdisk03/alvarezc/gangadir/workspace/palvare1/LocalXML/"
end=".root"

iter=""
subjobs=34

output="tree2Dl2.root"

dir="toy_result_ghi/toy_result2Dl2_"
nametree="toystudyHistBremCatTsallisBkg_results2DL0ETOSOnly_d.root"

for i in `seq 0 $subjobs`; do
   iter+=$dir$i"/"$nametree" "
done

string="hadd -k -f "$output" "$iter

echo "Merging..."
echo $string
$string
echo 
