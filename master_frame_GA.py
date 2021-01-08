'''
==================================================================
Master file for GA brace
==================================================================
'''

'''
======================================================
Import Part
======================================================
'''
from frame_GA import *
'''
======================================================
Parameters Part
======================================================
'''

# Structural Model
#======================================================
nx=3
nz=3
span=1
diameter= 0.1
lx=0
ly=-1000
lz=0
young=10000000000
c1=0
c2=0
c3=0
c4=0
c5=0
c6=0
c7=0
brace_type = None # 0 for brace to be truss element , 1 for brace to be beam element
forgame = 1002

# Population
#======================================================
num_gene = 100 # how many gene in the population

# GA
#======================================================
pop_multate_rate = 0.2 # how many gene will mutate
gene_multate_rate = 0.1 # how likely the gene will mutate
cross_rate = 0.2 # how many gene'codes to be swap
max_pop = num_gene # maximum population
survive_rate = 0.50 # how many gene will survive in each iteration
iteration = 200


'''
======================================================
Run Part
======================================================
'''

GA_pop = Population(objective,num_gene,nx,nz,span,diameter,lx,ly,lz,young,c1,c2,c3,c4,c5,c6,c7,forgame)
GA_test = GA(GA_pop,pop_multate_rate,gene_multate_rate,cross_rate,max_pop,survive_rate,iteration)
GA_test.run()
