
build: clean
	python setup.py build_ext --inplace

docs:
	sh docs.sh


clean:
	@rm -f ateam/*.c
	@rm -f ateam/*.o
	@rm -f ateam/*.so
	@rm -f ateam/arithmetic/*.c
	@rm -f ateam/arithmetic/*.o
	@rm -f ateam/arithmetic/*.so
	@rm -rf ./build


fp = ./test/output/profiles/metadata/.profiles/profile.json

test: FORCE
	# screen -dm vprof -r
	python test/matrices.py
	# vprof -c mh test/matrices.py --output-file $fp
	# vprof --input-file $fp

FORCE: