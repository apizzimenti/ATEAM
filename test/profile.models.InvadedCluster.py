
import pstats, cProfile
import pyximport
import ic

pyximport.install()

cProfile.runctx("ic.chain()", globals(), locals(), "Profile.prof")

with open('ic.txt', 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

