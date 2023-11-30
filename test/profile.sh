
fp=./output/profiles/metadata/profile.json
vprof -c mh sample.py --output-file $fp
vprof --input-file $fp

# cfp=./output/profiles/metadata/profile-$1.json
# vprof -c cph sample-$1.py --output-file $cfp
