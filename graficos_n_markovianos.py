from numpy import cos, sin, sqrt, pi, exp
import sys 
sys.path.append('src')
from src.theoric_channels import *
a = TheoricMaps()

lamb = 0.01
x1 = np.linspace(0,1,21)
#x = [i/max(x) for i in x]
th = np.pi/2
ph = np.pi/2
#y1 = a.bpf(x1, th, ph)
#ya = a.bpf(xa, th, ph)
#yb = a.bpf(xb, th, ph)

a.plot_storaged('adg',True)
a.plot_theoric(x1,'adg',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()

a.plot_storaged('ad',True)
a.plot_theoric(x1,'ad',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()

a.plot_storaged('pf',True)
a.plot_theoric(x1,'pf',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()

#a.plot_storaged('pd',True)
a.plot_theoric(x1,'pd',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()

a.plot_storaged('bf',True)
a.plot_theoric(x1,'bf',theta=pi/2,phi=pi/2,descript='isometria')
plt.legend(loc=1)
plt.show()

a.plot_storaged('bpf',True)
a.plot_theoric(x1,'bpf',theta=pi/2,phi=0.0,descript='isometria')
plt.legend(loc=1)
plt.show()
#s.exit()
a.plot_storaged('d',True)
a.plot_theoric(x1,'d',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()
a.plot_storaged('l',True)
a.plot_theoric(x1,'l',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()
a.plot_storaged('adg',True)
a.plot_theoric(x1,'adg',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()
a.plot_storaged('hw',True)
a.plot_theoric(x1,'hw',theta=pi/2,phi=0,descript='isometria')
plt.legend(loc=1)
plt.show()
