# IMPORTIAMO LE LIBRERIE NECESSARIE E I FILE OTTENUTI DAL CODICE IN C

from numpy import loadtxt
from pylab import xlabel,ylabel
import matplotlib.pyplot as plt

data = loadtxt("Dati_es_1",float)
databessel = loadtxt("Dati_es_1_Bessel",float)
dataneumann = loadtxt("Dati_es_1_Neumann",float)
datadelta = loadtxt("Dati_es_1_phase_shift",float)
datacross6 = loadtxt("Dati_es_1_cross_section6.txt",float)
datacross8 = loadtxt("Dati_es_1_cross_section8",float)
datadeltaquadro = loadtxt("Dati_es_1_Delta_squared",float)
datcrossperognil = loadtxt("Dati_es_1_cross_section_per_ogni_l", float)

#################################################################################################################

# STEP 1: PLOT DELLE SOLUZIONI DELL'EQUAZIONE DI SCHRODINGER

x0 = data[:10000,0]  # sotto-array per separare le soluzioni al variare di l
y0 = data[:10000,1]

x1 = data[10001:20000,0]
y1 = data[10001:20000,1]

x2 = data[20001:30000,0]
y2 = data[20001:30000,1]

x3 = data[30001:40000,0]
y3 = data[30001:40000,1]

x4 = data[40001:50000,0]
y4 = data[40001:50000,1]

x5 = data[50001:60000,0]
y5 = data[50001:60000,1]

x6 = data[60001:70000,0]
y6 = data[60001:70000,1]


fig, ax = plt.subplots(1, figsize=(8, 6))  # initialize the figure and axes

fig.suptitle('Solutions of the Schrödinger equation', fontsize=15) # set the title for the figure

ax.plot(x0,y0,"r",label="l=0") # Draw all the lines in the same plot, assigning a label for each one to be shown in the legend
ax.plot(x1,y1,"g",label="l=1")
ax.plot(x2,y2,"b",label="l=2")
ax.plot(x3,y3,"c",label="l=3")
ax.plot(x4,y4,"m",label="l=4")
ax.plot(x5,y5,"y",label="l=5")
ax.plot(x6,y6,"k",label="l=6")

ax.legend(loc="upper left", title="Legend", frameon=False) # Add a legend with title, position it (loc) with no box framing (frameon)

xlabel ("$r [\mathring{\mathrm{A}}]$") # label for the two axis
ylabel ("$u_{l}(r)$")

plt.show() # show the plot


#################################################################################################################

# STEP 2.1: PLOT DELLE FUNZIONI DI BESSEL

j_x = databessel[:,0]  # punti in cui calcolare le funzioni di Bessel
j_y0 = databessel[:,1]  # funzioni di Bessel al variare di l
j_y1 = databessel[:,2]
j_y2 = databessel[:,3]
j_y3 = databessel[:,4]
j_y4 = databessel[:,5]
j_y5 = databessel[:,6]
j_y6 = databessel[:,7]

figj, jx = plt.subplots(1, figsize=(8, 6))

figj.suptitle('Spherical Bessel functions', fontsize=15)

jx.plot(j_x,j_y0,"r",label="l=0")
jx.plot(j_x,j_y1,"g",label="l=1")
jx.plot(j_x,j_y2,"b",label="l=2")
jx.plot(j_x,j_y3,"c",label="l=3")
jx.plot(j_x,j_y4,"m",label="l=4")
jx.plot(j_x,j_y5,"y",label="l=5")
jx.plot(j_x,j_y6,"k",label="l=6")

jx.legend(loc="upper right", title="Legend", frameon=False)

xlabel ("$x$")
ylabel ("$j_{l}(x)$")
plt.xlim([0.1, 20])
plt.ylim([-0.25, 1.05])

plt.show()


#################################################################################################################

# STEP 2.2: PLOT DELLE FUNZIONI DI NEUMANN

n_x = dataneumann[:,0]  # punti in cui calcolare le funzioni di Neumann
n_y0 = dataneumann[:,1]  # funzioni di Neumann al variare di l
n_y1 = dataneumann[:,2]
n_y2 = dataneumann[:,3]
n_y3 = dataneumann[:,4]
n_y4 = dataneumann[:,5]
n_y5 = dataneumann[:,6]
n_y6 = dataneumann[:,7]

fign, nx = plt.subplots(1, figsize=(8, 6))

fign.suptitle('Spherical Neumann functions', fontsize=15)

nx.plot(n_x,n_y0,"r",label="l=0")
nx.plot(n_x,n_y1,"g",label="l=1")
nx.plot(n_x,n_y2,"b",label="l=2")
nx.plot(n_x,n_y3,"c",label="l=3")
nx.plot(n_x,n_y4,"m",label="l=4")
nx.plot(n_x,n_y5,"y",label="l=5")
nx.plot(n_x,n_y6,"k",label="l=6")
""" nx.plot(n_x,n_y7,"tab:orange",label="l=7")
nx.plot(n_x,n_y8,"tab:brown",label="l=8") """

nx.legend(loc="lower right", title="Legend", frameon=False)

xlabel ("$x$")
ylabel ("$n_{l}(x)$")
plt.xlim([0, 20])
plt.ylim([-1.5, 1])

plt.show()


#################################################################################################################

# STEP 2.3: PLOT DELLE PHASE SHIFT

r2 = datadelta[0:5,0]  # punti r2 della coppia (r1,r2) per calcolare le phase shift
delta0 = datadelta[0:5,1]  # phase shift al variare di l
delta1 = datadelta[0:5,2]
delta2 = datadelta[0:5,3]
delta3 = datadelta[0:5,4]
delta4 = datadelta[0:5,5]
delta5 = datadelta[0:5,6]
delta6 = datadelta[0:5,7]

figph, phx = plt.subplots(1, figsize=(8, 6))

figph.suptitle('Phase shifts', fontsize=15)

phx.plot(r2,delta0,"ro",label="l=0")
phx.plot(r2,delta1,"go",label="l=1")
phx.plot(r2,delta2,"bo",label="l=2")
phx.plot(r2,delta3,"co",label="l=3")
phx.plot(r2,delta4,"mo",label="l=4")
phx.plot(r2,delta5,"yo",label="l=5")
phx.plot(r2,delta6,"ko",label="l=6")

phx.legend(loc="upper left", title="Legend", frameon=False, bbox_to_anchor=(1.0, 1.0))

xlabel ("$r_{2} [\mathring{\mathrm{A}}]$")
ylabel ("$\delta_{l} [rad]$")

plt.show()


#################################################################################################################

# STEP 3: PLOT DELLE CROSS SECTION
energy7 = datcrossperognil[:350,0]

cross_section7_0  = datcrossperognil[:350,1]
cross_section7_1  = datcrossperognil[350:700,1]
cross_section7_2  = datcrossperognil[700:1050,1]
cross_section7_3  = datcrossperognil[1050:1400,1]
cross_section7_4  = datcrossperognil[1400:1750,1]
cross_section7_5  = datcrossperognil[1750:2100,1]
cross_section7_6  = datcrossperognil[2100:2450,1]
cross_section7_7  = datcrossperognil[2450:2800,1]
cross_section7_8  = datcrossperognil[2800:3150,1]

fige, ex = plt.subplots(1, figsize=(8, 6))
fige.suptitle('Contribution to the total cross section', fontsize=15)
ex.plot(energy7, cross_section7_0,"r", label="l=0")
ex.plot(energy7, cross_section7_1,"g", label="l=1")
ex.plot(energy7, cross_section7_2,"b", label="l=2")
ex.plot(energy7, cross_section7_3,"m", label="l=3")
ex.plot(energy7, cross_section7_4,"y", label="l=4")
ex.plot(energy7, cross_section7_5,"k", label="l=5")
ex.plot(energy7, cross_section7_6,"c", label="l=6")
ex.plot(energy7, cross_section7_7,"violet", label="l=7")
ex.plot(energy7, cross_section7_8,"palegreen", label="l=8")
xlabel ("$E [meV]$")
ylabel ("$\sigma_{tot}(E) [\mathring{\mathrm{A}}^{2}]$")
plt.xlim([0, 3.5])
plt.ylim([0, 460])
ex.legend(loc="upper left", title="Legend", frameon=False, bbox_to_anchor=(1.0, 1.0))
plt.show()




energy6 = datacross6[:350,0]  # valori di energia 
cross_section6 = datacross6[:350,1]  # valori di cross section per l_max = 6

energy8 = datacross8[:350,0]
cross_section8 = datacross8[:350,1] 

figc, cx = plt.subplots(1, figsize=(8, 6))

figc.suptitle('Total cross section', fontsize=15)

cx.plot(energy6, cross_section6,"b", label="l_max=6")
cx.plot(energy8, cross_section8, "r--", label="l_max=8")

cx.legend(loc="upper right", title="Legend", frameon=False)

xlabel ("$E [meV]$")
ylabel ("$\sigma_{tot}(E) [\mathring{\mathrm{A}}^{2}]$")

plt.show()




#################################################################################################################

# STEP 4: PLOT DEI DELTA^{2}

delta_x = datadeltaquadro[:,0]  # valori di sigma
delta_y = datadeltaquadro[:,1]  # valori di Delta^{2} corrispondenti

figd, dx = plt.subplots(1, figsize=(8, 6))

figd.suptitle('$\Delta^{2}$ for the parameter $\sigma$ of the Lennard-Jones potential', fontsize=15)

dx.scatter(delta_x,delta_y)

xlabel ("$\sigma [\mathring{\mathrm{A}}]$")
ylabel ("$\Delta^{2}(\sigma) [meV^{2}]$")

plt.show()