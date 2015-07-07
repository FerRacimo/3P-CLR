import dadi, numpy, scipy
from optparse import OptionParser

def TwoPop(params,ns,pts):
    driftC,driftA = params    
    Nhum = 1
    Narch = (driftC * 2 * Nhum) / (2 * driftA)
    earliest = max([driftC,driftA])*10
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi,xx,(earliest-driftC),nu=Nhum)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
    phi = dadi.Integration.two_pops(phi,xx,driftC,Nhum,Narch,m12=0,m21=0)
    fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
    fs = numpy.array(fs)
#    fs = fs[1:ns[0]]
    return(fs)


parser = OptionParser("$prog [options]")
parser.add_option("-c", "--driftC", dest="driftC", help="Drift C", default=None, type="float")
parser.add_option("-a", "--driftA", dest="driftA", help="Drift A", default=None, type="float")
parser.add_option("-m", "--nC", dest="nC", help="Number of samples from C", default=None, type="int")
parser.add_option("-b", "--nB", dest="nB", help="Number of samples from B", default=None, type="int")
(options,args) = parser.parse_args()


driftC = options.driftC
driftA = options.driftA
nC = options.nC
nB = options.nB

daditable = TwoPop([driftC,driftA],(nC,nB),20)

#print daditable

daditable.tolist()

for i in daditable:
    for j in i:
        print j
