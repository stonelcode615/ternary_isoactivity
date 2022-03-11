#! /usr/bin/python

from matplotlib import pyplot as plt
import scipy as sp
from scipy.optimize import differential_evolution
from scipy.optimize import newton
import sys
import optparse
import numpy as np
import warnings
import time

start_time = time.time()

np.set_printoptions(formatter={'float':'{:.4f}'.format})

alpha12=np.array([0.2,0.2])
alpha13=np.array([0.2,0.2])
alpha23=np.array([0.2,0.2])
a12=round(alpha12[0],1); a21=round(alpha12[1],1)
a13=round(alpha13[0],1); a31=round(alpha13[1],1)
a23=round(alpha23[0],1); a32=round(alpha23[1],1)
alpha = np.array([[0.000,a12,a13],[a21,0.000,a23],[a31,a32,0.000]])

Tem=273.15+10.00
RT = 8.3144598/1000*Tem

#============$file_reading==============
def file_reading(filename):
    file_handle = open(filename,'r')
    s1 = 'hello'
    index = 0
    row_item = []
    while s1 != '':
          s1 = file_handle.readline()
          s2 = s1.rstrip('\n').split(' ')
          try:
              float(s2[0])
          except ValueError:
                 if s2[0]=='':
                     continue
                 else:
                     topname = s2
                     continue
          content = [float(item) for item in s2]
          row_item.append(content)
    return topname, row_item

#============$gamma==============
def cal_gamma(tao, tie_x):
    warnings.simplefilter("error",RuntimeWarning)
    tao=np.insert(tao,0,0.0)
    tao=np.insert(tao,4,0.0)
    tao=np.insert(tao,8,0.0)
    Tao = tao.reshape((3,3))
    G = np.exp(-alpha*Tao)

    G1   = np.ones((3,3))
    G2   = np.ones((3,3))
    G3   = np.ones((3,3))
    DACT = np.ones((3,3))
    DDG  = np.ones((2,2))
    
    gamma0 = []
    for nd in range(len(tie_x)):
        gamma1 = []
        for phase in range(2):
            x = tie_x[nd][phase]
            for i in range(3):
                if x[i] < 1.0e-10:
                    x[i] = 1.0e-10
            r = np.array([np.dot(x,G[:,i])        for i in range(3)])
            s = np.array([np.dot(x,(Tao*G)[:,i]) for i in range(3)])
            sor = s/r
            ln_gamma = sor
            for i in range(3):
                G1[:,i] = G[:,i]/r[i]
                G2[:,i] = G1[:,i]*(Tao[:,i]-sor[i])
                G3[:,i] = x[i]*G2[:,i]
            for i in range(3):
                ln_gamma[i] = ln_gamma[i]+sum(G3[i,:])
            valid='success'
            try:
                GAM = np.exp(ln_gamma)
            except(RuntimeError,RuntimeWarning):
                valid = 'failed'
                break
            GAMlist = GAM.tolist()
            ACT = x*GAM
            for i in range(3):
                for j in range(3):
                    sumG2 = G2[i,j]+G2[j,i]
                    for k in range(3):
                        sumG2 = sumG2-G1[i,k]*G3[j,k]-G1[j,k]*G3[i,k]
                    DACT[i,j] = sumG2
                    DACT[j,i] = sumG2
            for i in range(3):
                for j in range(3):
                    DACT[i,j] = DACT[i,j]*ACT[i]
                    if j==i:
                        DACT[i,i] = DACT[i,i]+GAM[i]

            for i in range(3):
                for j in range(3):
                    try:
                        DACT[i,j] = DACT[i,j]/ACT[i]
                    except:
                        valid = 'failed'
                        break
                if valid is 'failed':
                    break

            for i in range(1,3):
                ii = i-1
                for j in range(1,3):
                    jj = j-1
                    DDG[ii,jj] = DACT[i,j]-DACT[0,j]-DACT[i,0]+DACT[0,0]
            DET=DDG[0,0]*DDG[1,1]-DDG[1,0]-DDG[1,0]
            if DET <= 0.0 or DDG[0,0] <= 0.0 or DDG[1,1] <= 0.0:
                valid = 'failed'
                #print 'phase %d unstable'%phase
                break   # break phase loop
            gamma1.append(GAM)
        if valid is 'failed':
            break       # break nd tie-line number loop
        gamma0.append(gamma1)
        gamma=np.array(gamma0)
    if valid is 'failed':
	return None
    else:
        return gamma

#============$tie==============
def cal_tie(tao,tie_x):
    gamma = cal_gamma(tao, tie_x)
    if gamma is None:
	    return None
    tie_line=[]
    valid='success'
    warnings.simplefilter("error",RuntimeWarning)
    for nd in range(len(tie_x)):
        try:
	    K = [gamma1/gamma2 for gamma1,gamma2 in zip(gamma[nd][0],gamma[nd][1])]
        except (RuntimeError,RuntimeWarning):
            valid='failed'
#	    print 'K=gamma1/gamma2 ', valid
            break
        l = 0.5
        Z = [l*phase2_x+(1.0-l)*phase1_x for phase1_x,phase2_x in zip(tie_x[nd][0],tie_x[nd][1])]
        def equ_x(x,K,Z):
            return sum([(k-1)*z/(1+(k-1)*x) for k,z in zip(K,Z)])
        try:
            alpha_c = newton(equ_x,0.5,args=(K,Z))
        except (RuntimeError,RuntimeWarning):
            valid='failed'
#	    print 'newtown ', valid
            break
        try:
	    x1 = [  z/(k*alpha_c+(1-alpha_c)) for k,z in zip(K,Z)]
	    x2 = [k*z/(k*alpha_c+(1-alpha_c)) for k,z in zip(K,Z)]
	except (RuntimeError,RuntimeWarning):
            valid = 'failed'
#	    print 'x1 or x2 ', valid
            break
        tie_line.append([x1, x2])
    if valid is 'failed':
        return None
    else:
        return tie_line


#===========$f2===============================
def obj_f2(tao, exp_x):
    sum1 = 0.0
    sum2 = 0.0
    cal_x = cal_tie(tao,exp_x)
    if cal_x is None:
        return 1000.0
    gamma = cal_gamma(tao,cal_x)
    if gamma is None:
        return 1000.0
    else:
        for nd in range(len(exp_x)):
            sum1 += sum([np.power((i-j),2) for i,j in zip(exp_x[nd][0],cal_x[nd][0])])
            sum1 += sum([np.power((i-j),2) for i,j in zip(exp_x[nd][1],cal_x[nd][1])])
        return sum1

#################################

component, row = file_reading('exp-tie.dat')

exp_tielines = []
for item in row:
    x1_ph1 = round(1.0-item[0]-item[1],4)
    x1_ph2 = round(1.0-item[2]-item[3],4)
    phase1=[]
    phase2=[]
    tielines=[]
    phase1 = [x1_ph1, item[0], item[1]]
    phase2 = [x1_ph2, item[2], item[3]]
    exp_tielines.append([phase1, phase2])
exp_x = exp_tielines[:]
print exp_x

f_tie = open("fit-tie.dat","w")
f_tie.write(component[0]+' '+component[1]+' '+component[2]+'\n')
bounds = [(-4.0,6.0)]*6
result = differential_evolution(obj_f2, bounds,(exp_x,),popsize=50)
tao=result.x
print 'DE tao is ', tao, 'with sum square ',result.fun

sum1 = 0.0
sum2 = 0.0
cal_x = cal_tie(tao, exp_x)
for nd in range(len(cal_x)):
    sum1 += sum([np.power((i-j),2) for i,j in zip(exp_x[nd][0],cal_x[nd][0])])
    sum1 += sum([np.power((i-j),2) for i,j in zip(exp_x[nd][1],cal_x[nd][1])])
    f_tie.write('%.4f %.4f %.4f %.4f\n'%(cal_x[nd][0][1],cal_x[nd][0][2],\
            cal_x[nd][1][1],cal_x[nd][1][2]))
f_tie.close()

print 'x     fitting sum2 is %.5f' % sum1
print 'x rmsd is %.5f' %(sum1/(6*len(cal_x)))**(0.5)


print '%.3f  %.3f  %.1f' %(tao[0]*Tem,tao[2]*Tem,a12)
print '%.3f  %.3f  %.1f' %(tao[1]*Tem,tao[4]*Tem,a13)
print '%.3f  %.3f  %.1f' %(tao[3]*Tem,tao[5]*Tem,a23)


print "%f seconds" %(time.time()-start_time)


