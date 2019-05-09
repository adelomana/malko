import os,time,numpy
import multiprocessing,multiprocessing.pool

def launcher(iteration):

    print(iteration)
    time.sleep(3)

    return None

###
### MAIN
###

iterations=numpy.arange(32)
numberOfThreads=8

hydra=multiprocessing.pool.Pool(numberOfThreads)

results=hydra.map(launcher,iterations)
