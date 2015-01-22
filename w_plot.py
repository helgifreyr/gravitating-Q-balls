from w_data import get_data
import sys
from pylab import *

dir = sys.argv[1]

data = get_data(dir)
w = data[:,1]
M = data[:,5]*100
plot(w,M,'.',ms=2.0)
# xlim(0.31,0.40)
yscale('log')
xlabel(r'$\omega$')
ylabel(r'$M$')
savefig(dir.replace('/','-')+'w_plot.png')
