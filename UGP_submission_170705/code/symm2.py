# The following program checks if a functions symmetric block multilinear representation can be not bounded by 1.
# We check this over all functions on n variables
# A function is a counterexample if its block multilinear rep is not bounded

from sys import argv
import math
import itertools
from mip import *

#from random import seed
#from random import randint

print("Give the number of variables")
n=int(input()) # number of variables
size_tt = 2**n  # size of the truth table of f
t=0 # degree of the polynomials
n_bm = t * (n + 1)
count = 0 # global variable to count number of counterexamples


def  countSetBits(n): 
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
        #print(m,k)
        if (k == 0):
            rep.append(-1)
        else:
            rep.append(1)
        m = int(m / 2)
    return rep

#print(binary(63))


def is_func(func):        # is this function a counterexample, needs input in block multilinear representation
    m=Model()
    global t
    print("t: "+str(t))
    permutations_t=list(itertools.permutations(list(range(t))))
    print("permutation")
    print(permutations_t)
    x_b=[]
    x_ext=[]
    x_sum=[]
    t_fact=math.factorial(t)
    for x in range(n+1):
        x_b.append( [m.add_var(name="x_b_"+str(x)+"_"+str(i),var_type=BINARY) for i in range(t)] )
    for x in range(size_tt):
        x_ext.append([m.add_var(name="x_ext_"+str(x)+"_"+str(i),var_type=BINARY) for i in range(t_fact)])
        x_sum.append([m.add_var(name="x_sum_"+str(x)+"_"+str(i),var_type=INTEGER) for i in range(t_fact)])
        vars=[]
        temp=x
        cnt=0
        for i in range(n):
            if 1&temp==1:
                vars.append(i+1)
                cnt+=1
            temp//=2
        vars+=[0]*(t-cnt)
        for y in range(t_fact):
            m+=(x_ext[x][y]+2*x_sum[x][y]==xsum(x_b[vars[z]][permutations_t[y][z]] for z in range(t)))

    m.objective=xsum( func[x][y]*(1-2*x_ext[x][y]) for x in range(size_tt) for y in range(t_fact) )
    status = m.optimize(max_seconds=500)
    if status == OptimizationStatus.OPTIMAL:
        print('Optimal value is'+str(m.objective_value))
        if m.objective_value<-1-1e-5:
            print('Counter eample found with {} as its minimum value'.format(m.objective_value))
            for k in x_ext:
                # if abs(v.x) > 1e-6: # only printing non-zeros
                for v in k:
                    print('{} : {}'.format(v.name, v.x))
            print(func)
            exit()
            # exit()

def chi(i,j):  # Calculates \chi_i(j) , where i should be viewed as a set (-1 means index not present) and j should be viewed as an assignment.
    value = 1
    bin_j = binary(j,n)
    bin_i = binary(i,n)
    for k in range(0,n):
        if (bin_i[k] == 1):
            value = value * bin_j[k]
    #print(value, bin_i, bin_j)
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
    #print("the poly rep is ", func_poly)
    
    t_fact=math.factorial(t)
    h = 1/(t_fact)
    func_bm=[]
   
    for i in range(size_tt):
        func_bm.append([0]*t_fact)
        for j in range(t_fact):
            func_bm[i][j]=func_poly[i]*h
    print(func_bm)
    return func_bm


for i in range(0,2**(size_tt)):
    f_tt = binary(i, size_tt)
    print("The truth table of the function is ", f_tt)
    f_bm = poly(f_tt)
    #print("Its bm poly representation is ", f_bm)
    is_func(f_bm)

    print("count is ", count)