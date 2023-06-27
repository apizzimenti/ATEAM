
fp=./output/profiles/metadata/$1.json
vprof -c cph test-profile.py --output-file $fp
vprof --input-file $fp
