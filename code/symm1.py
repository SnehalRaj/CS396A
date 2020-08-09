# The following program checks if a functions symmetric block multilinear representation can be not bounded by 1.
# We check this over all functions on n variables
# A function is a counterexample if its block multilinear rep is not bounded

from sys import argv
import math
import itertools
from mip import *

#from random import seed
#from random import randint
f = open("output.txt", "w")
print("Give the number of variables")
n=int(input()) # number of variables
size_tt = 2**n  # size of the truth table of f
t=0 # degree of the polynomials
n_bm = t * (n + 1)
count = 0 # global variable to count number of counterexamples


def  countSetBits(n): # Count Number of bits which are 1 in n
    count = 0
    while (n): 
        count += n & 1
        n >>= 1
    return count 

def binary(l,dig): #binary representation of l (0 replaced with -1) with dig number of digits   
    rep = []
    m=int(l)
    for i in range(0,dig):
        k = m%2
        #f.write(m,k)
        if (k == 0):
            rep.append(-1)
        else:
            rep.append(1)
        m = int(m / 2)
    return rep

#f.write(binary(63))
def symmetricSplit(cnt,n): # Split 1 into t parts
    i=cnt%2
    if cnt//2==0:
        return [[0]*cnt]*(n+1)
    else:
        split= [ [j]*((cnt//2)*2) for j in range(n+1) ]
        if i:
            for j in range(len(split)):
                split[j]=[0]+split[j]
        return split


def is_func(func,n,t): # is this function a counterexample, needs input in block multilinear representation
    m=Model()
    f.write("Number of blocks: {}".format(str(t)))
    permutations_t=list(itertools.permutations(list(range(t))))
    x_b=[]
    x_ext=[]
    x_sum=[]
    t_fact=math.factorial(t)
    for x in range(n+1):
        x_b.append( [m.add_var(name="x_b_"+str(x)+"_"+str(i),var_type=BINARY) for i in range(t)] )
    for x in range(size_tt):
        x_ext.append([[m.add_var(name="x_ext_"+str(x)+"_"+str(j)+"_"+str(i),var_type=BINARY) for i in range(t_fact)] for j in range(n+1)] ) # This will store the value of monomial's
        x_sum.append([[m.add_var(name="x_sum_"+str(x)+"_"+str(j)+"_"+str(i),var_type=INTEGER) for i in range(t_fact)] for j in range(n+1)])
        vars=[]
        temp=x
        cnt=0
        for i in range(n):
            if 1&temp==1:
                vars.append(i+1)
                cnt+=1
            temp//=2
        if cnt>t:
            continue
        split=symmetricSplit(t-cnt)
        
        for i in range(n+1):
            temp=vars.copy()+split[i]
            for j in range(t_fact):
                m+=(x_ext[x][i][j]+2*x_sum[x][i][j]==xsum(x_b[temp[z]][permutations_t[j][z]] for z in range(t))) # Adding Constraints of the ILP

    m.objective=xsum( (1/(n+1))*func[x][z]*(1-2*x_ext[x][y][z]) for x in range(size_tt) for y in range(n+1) for z in range(t_fact) )
    status = m.optimize(max_seconds=500)
    if status == OptimizationStatus.OPTIMAL:
        f.write('Optimal value is {}'.format(str(m.objective_value)))
        if m.objective_value<-1-1e-5:
            f.write('Counter example found with {} as its minimum value'.format(m.objective_value))
            for v in m.vars:
                f.write('{} : {}'.format(v.name, v.x))
            exit()

def chi(i,j):  # Calculates \chi_i(j) , where i should be viewed as a set (-1 means index not present) and j should be viewed as an assignment.
    value = 1
    bin_j = binary(j,n)
    bin_i = binary(i,n)
    for k in range(0,n):
        if (bin_i[k] == 1):
            value = value * bin_j[k]
    #f.write(value, bin_i, bin_j)
    return value


def poly(func_tt):   # give the block multilinear rep of f from its truth table representation, the first entry of truth table is -1,-1,..., -1
    global t
    t=0
    func_poly = [0] * size_tt
    for i in range(0,size_tt):
        func_poly[i] = 0
        for j in range(0,size_tt):
            func_poly[i] = func_poly[i] + func_tt[j] * chi(i,j)
        func_poly[i] = (1.0 / size_tt) * func_poly[i]   # ith Fourier coefficient
        if abs(func_poly[i])>1e-6:
            t=max(countSetBits(i),t)    
    t_fact=math.factorial(t)
    h = 1/(t_fact)
    func_bm=[]
   
    for i in range(size_tt):
        func_bm.append([0]*t_fact)
        for j in range(t_fact):
            func_bm[i][j]=func_poly[i]*h
    return func_bm


for i in range(0,2**(size_tt)):
    f_tt = binary(i, size_tt)
    f.write("The truth table of the function is {}".format(str(f_tt)))
    f_bm = poly(f_tt)
    f.write("Its blockmultilinear poly representation is {}".format( str(f_bm)))
    is_func(f_bm,n,t)

    f.write("count is {}".format( count))
f.close()