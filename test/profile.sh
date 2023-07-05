
fp=./output/profiles/metadata/profile.json
vprof -c cph sample.py --output-file $fp

cfp=./output/profiles/metadata/profile-$1.json
vprof -c cph sample-$1.py --output-file $cfp
