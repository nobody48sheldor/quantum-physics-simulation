import matplotlib.pyplot as plt
import numpy as np
import cmath
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
#import cython
import os
from matplotlib import style

style.use('dark_background')

L = 20
n = 300
nt = int(200*n/4)
alpha = 4
m = 1
h = 6.62607015 * (10**(-34))
hbar = 1

T = 2*L/(400)
print(T)

x =np.linspace(-L, L, n)
y =np.linspace(-L, L, n)
t =np.linspace(0, 20, nt)

dx = x[1] - x[0]
dy = y[1] - y[0]
dt = t[1] - t[0]

print(dt/(dx*dx))

X, Y = np.meshgrid(x, y)

def psi_0(x, y):
    S = []
    Sr = []
    Si = []
    S2 = []
    for w in y:
        Sx = []
        Sxr = []
        Sxi = []
        Sx2 = []
        for i in x:
            s = (1-1j) * ((np.exp(-((i-L/2)**2 + (w-L/2)**2)/alpha)) / (np.sqrt(2)*np.pi*alpha))
            Sx.append(s)
            Sxr.append(s.real)
            Sxi.append(s.imag)
            Sx2.append(s.real*s.real + s.imag*s.imag)
        S.append(np.array(Sx))
        Sr.append(np.array(Sxr))
        Si.append(np.array(Sxi))
        S2.append(np.array(Sx2))
    return(Sr, Si, S, S2)

psi0 = psi_0(x, y)
psi0r = np.array(psi0[0])
psi0i = np.array(psi0[1])
psi0N = np.array(psi0[2])
psi02 = np.array(psi0[3])

fig = plt.figure()

axa = plt.subplot2grid((1,2), (0,0), rowspan = 1, colspan = 1, projection = '3d')
axa.plot_surface(X, Y, np.array(psi0r), cmap = cm.Reds, linewidth=0, antialiased=True, alpha = 0.6)
axa.plot_surface(X, Y, np.array(psi0i), cmap = cm.Blues, linewidth=0, antialiased=True, alpha = 0.6)
axa.plot_surface(X, Y, np.array(psi02), cmap = cm.plasma, linewidth=0, antialiased=True, alpha = 0.8)
axa.axes.set_xlim3d(left=-L, right=L)
axa.axes.set_ylim3d(bottom=-L, top=L)
axa.axes.set_zlim3d(bottom=-0.4, top=0.4)

def V(x, y):
    S = []
    for w in y:
        Sx = []
        for i in x:
            s = 0.5*(i*i + w*w)
            Sx.append(s)
        S.append(np.array(Sx))
    return(S)

V = V(x, y)

axb = plt.subplot2grid((1,2), (0,1), rowspan = 1, colspan = 1, projection = '3d')
axb.plot_surface(X, Y, np.array(V), cmap = cm.plasma, linewidth=0, antialiased=True)
axb.axes.set_xlim3d(left=-L, right=L)
axb.axes.set_ylim3d(bottom=-L, top=L)
axb.axes.set_zlim3d(bottom=-1, top= 0.5*L**2 + 1)
plt.show()

PSI = []
PSIR = []
PSII = []
PSI2 = []

PSI.append(psi0N)
PSI.append(psi0N)

PSIR.append(psi0r)
PSIR.append(psi0r)

PSII.append(psi0i)
PSII.append(psi0i)

PSI2.append(psi02)
PSI2.append(psi02)

psix = []
psixR = []
psixI = []
psix2 = []

psiy = []
psiyR = []
psiyI = []
psiy2 = []

psi = 0

def plot():
    fig = plt.figure()
    plt.ion()

    for i in range(nt):
        plt.clf()

        ax1 = plt.subplot2grid((5,5), (0,0), rowspan = 3, colspan = 3, projection = '3d')
        ax2 = plt.subplot2grid((5,5), (0,3), rowspan = 2, colspan = 3, projection = '3d')
        ax3 = plt.subplot2grid((5,5), (3,0), rowspan = 3, colspan = 2, projection = '3d')
        ax4 = plt.subplot2grid((5,5), (3,3), rowspan = 2, colspan = 2, projection = '3d')

        ax1 = plt.subplot2grid((4,4), (0,0), rowspan = 3, colspan = 3, projection = '3d')
        ax2 = plt.subplot2grid((4,4), (3,0), rowspan = 1, colspan = 2, projection = '3d')
        ax3 = plt.subplot2grid((4,4), (3,2), rowspan = 1, colspan = 2, projection = '3d')
        ax4 = plt.subplot2grid((4,4), (0,3), rowspan = 3, colspan = 1, projection = '3d')

        ax1.plot_surface(X, Y, PSIR[i], cmap = cm.Reds, linewidth=0, antialiased=True, alpha = 0.6)
        ax1.plot_surface(X, Y, PSII[i], cmap = cm.Blues, linewidth=0, antialiased=True, alpha = 0.6)
        ax1.plot_surface(X, Y, PSI2[i], cmap = cm.plasma, linewidth=0, antialiased=True, alpha = 0.8)
        ax1.axes.set_xlim3d(left=-L, right=L)
        ax1.axes.set_ylim3d(bottom=-L, top=L)
        ax1.axes.set_zlim3d(bottom = -0.1, top = 0.1)
        ax1.set_title("all together")
        ax1.axes.set_axis_off()

        ax2.plot_surface(X, Y, PSIR[i], cmap = cm.Reds, linewidth=0, antialiased=True, alpha = 0.9)
        ax2.axes.set_xlim3d(left=-L, right=L)
        ax2.axes.set_ylim3d(bottom=-L, top=L)
        ax2.axes.set_zlim3d(bottom = -0.1, top = 0.1)
        ax2.set_title("real part")
        #ax2.axes.set_axis_off()

        ax3.plot_surface(X, Y, PSII[i], cmap = cm.Blues, linewidth=0, antialiased=True, alpha = 0.9)
        ax3.axes.set_xlim3d(left=-L, right=L)
        ax3.axes.set_ylim3d(bottom=-L, top=L)
        ax3.axes.set_zlim3d(bottom = -0.1, top = 0.1)
        ax3.set_title("imaginary part")
        #ax3.axes.set_axis_off()

        ax4.plot_surface(X, Y, PSI2[i], cmap = cm.plasma, linewidth=0, antialiased=True, alpha = 0.9)
        ax4.axes.set_xlim3d(left=-L, right=L)
        ax4.axes.set_ylim3d(bottom=-L, top=L)
        ax4.axes.set_zlim3d(bottom = -0.1, top = 0.1)
        ax4.set_title("probability density")
        ax4.axes.set_axis_off()

        plt.title('Time = {} secondes'.format(t[i]))

        #if i%20 == 0:
            #plt.savefig("{}.png".format(i))
        plt.pause(dt)

for i in range(1, nt - 2):
    os.system("clear")
    print(int(100*i/nt), " %")
    print(i, " / ", nt)
    psiy = []
    psiyR = []
    psiyI = []
    psiy2 = []

    psiy.append(PSI[i][0])
    psiyR.append(PSIR[i][0])
    psiyI.append(PSII[i][0])
    psiy2.append(PSI2[i][0])

    for v in range(1, n-1):
        psix = []
        psixR = []
        psixI = []
        psix2 = []

        psix.append(PSI[i][0][0])
        psixR.append(PSIR[i][0][0])
        psixI.append(PSII[i][0][0])
        psix2.append(PSI2[i][0][0])

        for u in range(1, n-1):
            try:
                psi = PSI[i-1][v][u] + ((2*1j/hbar) * ( ((hbar*hbar/(2*m)) * (0.5*(PSI[i][v][u+1] - 2*PSI[i][v][u] + PSI[i][v][u-1] + PSI[i-1][v][u+1] - 2*PSI[i-1][v][u] + PSI[i-1][v][u-1])/(2*dx*dx) + 0.5*(PSI[i][v+1][u] - 2*PSI[i][v][u] + PSI[i][v-1][u] + PSI[i-1][v+1][u] - 2*PSI[i-1][v][u] + PSI[i-1][v-1][u])/(2*dy*dy)) + V[v][u]*PSI[i][v][u]))*dt)
                psix.append(psi)
                psixR.append(psi.real)
                psixI.append(psi.imag)
                psix2.append(psi.real*psi.real + psi.imag*psi.imag)
            except KeyboardInterrupt:
                plot()

        psix.append(psi)
        psixR.append(psi.real)
        psixI.append(psi.imag)
        psix2.append(psi.real*psi.real + psi.imag*psi.imag)

        psiy.append(np.array(psix))
        psiyR.append(np.array(psixR))
        psiyI.append(np.array(psixI))
        psiy2.append(np.array(psix2))

    psiy.append(np.array(psix))
    psiyR.append(np.array(psixR))
    psiyI.append(np.array(psixI))
    psiy2.append(np.array(psix2))

    PSI.append(np.array(psiy))
    PSIR.append(np.array(psiyR))
    PSII.append(np.array(psiyI))
    PSI2.append(np.array(psiy2))

input("press enter")
plot()
