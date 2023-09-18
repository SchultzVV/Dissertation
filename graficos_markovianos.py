from numpy import cos, sin, sqrt, pi, exp
import sys 
sys.path.append('src')
from src.theoric_channels import *
a = TheoricMaps()

lamb = 0.01
list_1 = np.linspace(0.01,1,21)
list_2 = np.linspace(0.01,1000,21)


list_of_maps = ['bf','pf','l']
list_of_maps = ['l']
list_of_maps = ['ad','pd','adg','bf','bpf','d','l','hw']
list_of_maps = ['ad','pd','adg','bf','pf','bpf','d']
list_of_maps = ['pd','adg','bf','pf','bpf','d','l','hw']
list_of_maps = ['ad']
th = pi/2
ph = 0
for map in list_of_maps:
    if map == 'bf':
        ph = pi/2
    else:
        ph = 0
    a.plot_theoric(list_1,map,theta=th,phi=ph,descript='Teórico Markoviano')
    a.plot_storaged(map,True)
    a.plot_theoric_n_Markov(list_2,map,theta=th,phi=ph,descript='Teórico não Markoviano')
    # a.plot_theoric_n_Markov_B(x1,map,theta=th,phi=ph,descript='Teórico não Markoviano')
    a.plot_storaged(map,False)
    if map == 'l':
        plt.xlabel(fr'$\xi$ ; t (n-Markov)')
    else:
        plt.xlabel('p (Markov) ; t (n-Markov)')
    plt.ylabel('coerência')

    plt.xscale('log')
    plt.xlim(0.01)
    plt.legend(loc=0)
    plt.show()
sys.exit()    

a.plot_theoric(x1,'ad',theta=pi/2,phi=0,descript='Teórico Markoviano')
a.plot_storaged('ad',True)
a.plot_theoric_n_Markov(x1,'ad',theta=pi/2,phi=0,descript='Teórico não Markoviano')
a.plot_storaged('ad',False)
plt.xlabel('p (Markov) ; t (n-Markov)')
plt.ylabel('coerência')

plt.xscale('log')

plt.legend(loc=1)
plt.show()
sys.exit()
a.plot_theoric(x1,'pd',theta=pi/2,phi=0,descript='Phase-Damping')
a.plot_storaged('pd',False)
plt.legend(loc=1)
plt.show()

# 
a.plot_theoric(x1,'adg',theta=pi/2,phi=0,descript='Amplitude damping generalizado')
a.plot_storaged('adg',False)
plt.legend(loc=1)
plt.show()

a.plot_theoric(x1,'bf',theta=pi/2,phi=pi/2,descript='Bit-Flip')
a.plot_storaged('bf',False)
plt.legend(loc=1)
plt.show()

a.plot_theoric(x1,'pf',theta=pi/2,phi=0,descript='Phase-Flip')
a.plot_storaged('pf',False)
plt.legend(loc=1)
plt.show()

a.plot_theoric(x1,'bpf',theta=pi/2,phi=0.0,descript='Bit-Phase-Flip')
a.plot_storaged('bpf',False)
plt.legend(loc=1)
plt.show()
#s.exit()

a.plot_theoric(x1,'d',theta=pi/2,phi=0,descript='Depolarizing')
a.plot_storaged('d',False)
plt.legend(loc=1)
plt.show()

a.plot_theoric(x1,'l',theta=pi/2,phi=0,descript='Lorentz')
a.plot_storaged('l',False)
plt.legend(loc=1)
plt.show()

#plt.show()
a.plot_theoric(x1,'hw',theta=pi/2,phi=0,descript='H-W dephasing')
a.plot_storaged('hw',False)
plt.legend(loc=1)

plt.show()

