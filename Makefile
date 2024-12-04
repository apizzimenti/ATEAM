
build: clean
	python setup.py build_ext --inplace
	sh docs.sh


clean:
	@rm -f potts/*.c
	@rm -f potts/*.o
	@rm -f potts/*.so
	@rm -f potts/arithmetic/*.c
	@rm -f potts/arithmetic/*.o
	@rm -f potts/arithmetic/*.so
	@rm -rf ./build


fp = ./test/output/profiles/metadata/.profiles/profile.json

test: FORCE
	# screen -dm vprof -r
	python test/matrices.py
	# vprof -c mh test/matrices.py --output-file $fp
	# vprof --input-file $fp

FORCE: