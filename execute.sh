#!/usr/bin/env bash


cd src

start_time=$(date +%s)

echo "Computing edge weights..."
python compute_edge_weights.py -t 0.7 -ts 0.25

echo "Extract features..."
python feature_extraction.py -d 48 -p 4 -q 1 -wl 80 -nw 400 -k 89

echo "Processing clusters..."
python Module_processing.py -dn 1 -k 3 -ns 2500 -ne 100 -step 100
echo "DONE!!!"

end_time=$(date +%s)
Total_time=$((end_time - start_time))
echo "Total time: $Total_time seconds"



