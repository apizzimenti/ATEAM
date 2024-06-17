
build: clean
	python setup.py build_ext --inplace


clean:
	@rm -f potts/arithmetic/*.c
	@rm -f potts/arithmetic/*.*o
	@rm -rf ./build


fp = ./test/output/profiles/metadata/profile.json

test: FORCE
	# screen -dm vprof -r
	python test/matrices.py
	# vprof -c mh test/matrices.py --output-file $fp
	# vprof --input-file $fp


FORCE: