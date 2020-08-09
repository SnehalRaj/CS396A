import numpy as np
from mip import *

poly=[[0,1,2],[1,2,3]]
coff=[1,2]
sz=0
for term in poly:
    t=max(term)
    sz=max(sz,t)
sz+=1
adj=[set() for _ in range(sz)]

for term in poly:
    for i in term:
        for j in term:
            adj[i].add(j)
        
for i in range(sz):
    adj[i].remove(i)

m=Model()
y=[m.add_var(name='y'+str(i),ub=1) for i in range(sz)]
c=[[m.add_var(var_type=BINARY, name='C'+str(i)+':'+str(j)) for i in range(sz)] for j in range(sz)]

for v in range(sz):
    m+=xsum(c[v][j] for j in range(sz))==1

for j in range(sz):
    for v in range(sz):
        for vp in adj[v]:
            m+=c[v][j]+c[vp][j]<=y[j]


m.objective=xsum(y[i] for i in range(sz))
status = m.optimize(max_seconds=300)
if status == OptimizationStatus.OPTIMAL:
    print('optimal solution cost {} found'.format(m.objective_value))
elif status == OptimizationStatus.FEASIBLE:
    print('sol.cost {} found, best possible: {}'.format(m.objective_value, m.objective_bound))
elif status == OptimizationStatus.NO_SOLUTION_FOUND:
    print('no feasible solution found, lower bound is: {}'.format(m.objective_bound))
if status == OptimizationStatus.OPTIMAL or status == OptimizationStatus.FEASIBLE:
    print('solution:')
    for v in m.vars:
        print('{} : {}'.format(v.name, v.x))