
# python -m cProfile -o test-profile.pstats test-profile.py
# gprof2dot --depth=2 -f pstats test-profile.pstats | dot -Tpng -o output/figures/test-profile.png
python test-profile.py > ./output/test-profile.txt
