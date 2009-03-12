def saveconfig(p):
    diffusion = p[0]
    specfile = open("species.txt", "w")
    print >>specfile, """e-      7e10
He      6e14
H3p     2e10
H3o     5e10"""
    specfile.close()

    reactfile = open("reactions.txt", "w")
    print >>reactfile, """#recombination
7e-7    H3p e- => H2 H
5.7e-8  H3o e- => H2 H

#ternary recombination
2.5e-25   H3p e- He => H2 H He
0.5e-25   H3o e- He => H2 H He

#state change
2e-12   H3p He => H3o He
2e-12   H3o He => H3p He

#diffusion
""", diffusion, """     H3p => H3p_wall
""", diffusion, """     H3o => H3o_wall
""", diffusion, """     e- => e-_wall
#1e3     H3p_wall e- => H2 H"""
    reactfile.close()

    config = open("config.txt", "w")
    print >>config, """outfile = result.dat
specfile = species.txt
reactfile = reactions.txt
plotfile = plotcmd.gnuplot"""
    config.close()

saveconfig([1e3])

import os
os.system("./kyslik")

from pylab import load
result = load("result.dat")

x = result[:,0]
y = result[:,1:]
from scipy.interpolate import interp1d
resfunc = interp1d(x,y,axis=0)
#print result
x = x[:-1]
y = y[:-1,:]


#generate noisy sample
from scipy import random, hstack, c_
y1 = random.normal(y[:,0], y[:,0]*5e-1)
y2 = random.normal(y[:,2], y[:,2]*5e-1)
sample = hstack((c_[x], c_[y1], c_[y2]))

#construct the fitfunc
def fitfuncgen(p):
    saveconfig(p)
    import os
    os.system("./kyslik")

    from pylab import load
    result = load("result.dat")

    x = result[:,0]
    y = result[:,1:]
    from scipy.interpolate import interp1d
    return interp1d(x,y,axis=0, bounds_error=False)

def errfunc(p, sample):
    fitfunc = fitfuncgen(p)

    from scipy import r_
    xsample = sample[:,0]
    fitvalues = fitfunc(sample[:,0])
    err = r_[fitvalues[:,0]-sample[:,1]]
    err1 = r_[(fitvalues[:,0]-sample[:,1])*1.0/(1e7+sample[:,1])]
    err2 = r_[(fitvalues[:,2]-sample[:,2])*1.0/(1e7+sample[:,1])]
    err = r_[err1, err2]
    #print fitvalues[100:105,:]
    print p, sum(abs(err))
    return err

from scipy.optimize import leastsq
p1, tmp, tmp, msg, success = leastsq(errfunc, [9e2], args = (sample,), full_output=1 )
#p1, tmp, tmp, msg, success = leastsq(errfunc, [9e2], args = (sample,), full_output=1 , epsfcn)
print p1
print success, msg

import pylab as p
p.semilogy(x, y1, 'o')
p.semilogy(x, y2, 'o')

myfunc = fitfuncgen([6e3])
p.semilogy(x, myfunc(x)[:,0])
p.semilogy(x, myfunc(x)[:,2])
p.semilogy(x, y[:,0], '.')
p.semilogy(x, y[:,2], '.')
#p.semilogy(x, abs(errfunc([1e3], sample)))
p.xlim(xmax=.5e-2)
p.ylim(ymin=1e7, ymax=2e11)
#plot(x, y)
p.show()
