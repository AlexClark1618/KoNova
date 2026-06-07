import numpy as np
import matplotlib.pyplot as plt
from statistics import NormalDist 

overlap1 = NormalDist(mu=430, sigma=(75/2.355)).overlap(NormalDist(mu=425, sigma=(50/2.355)))

print(overlap1)

overlap2 = 0.4*NormalDist(mu=475, sigma=(100/2.355)).overlap(NormalDist(mu=450, sigma=(300/2.355)))
print(overlap2)

print("overlap:", overlap1*overlap2)

photons_dep = 20000
fiber_sa = 1/38.1
r = 0.94

tot= 0
tri_tot=0
tot_list = []
tri_tot_list = []
tri_len = []
for path_len in np.arange(0, 20, 0.5): #Varying path length
    tri_path_len = (2*(19.05-(np.sqrt(19.05**2-(path_len/2)**2))))/np.sqrt(3)
    print("path_len", path_len)
    tri_len.append(tri_path_len)
    print("tri_len", tri_path_len)
    tri_initial_photons = photons_dep * tri_path_len/10 * 0.5 *fiber_sa 
    initial_photons = photons_dep * path_len/10 * 0.5 *fiber_sa 

    for i in range(51):
        tot += initial_photons * ((r**i * (1-fiber_sa)**i) + (r**(i+1)*(1-fiber_sa)**i))
        tri_tot += tri_initial_photons * ((r**i * (1-fiber_sa)**i) + (r**(i+1)*(1-fiber_sa)**i))

    print("Tot:", tot*(1)*0.0173)
    print("Tri Tot:", tri_tot*(1)*0.0173)

    tot_list.append(tot*(1)*0.0173)
    tri_tot_list.append(tri_tot*(1)*0.0173)

    tot = 0
    tri_tot = 0

plt.scatter(np.arange(0, 20, 0.5), tot_list)
plt.show()
plt.scatter(tri_len, tri_tot_list)

plt.ylabel("Photons Collected at SiPM")
plt.xlabel("Muon Path Length (mm)")
plt.show()

#5.5% capture efficiency spread across emission spectra
#10,000 photons produced across scint emission spectra
#


#tot1=initial_photons
#tot2=initial_photons * r
'''
for i in range(50):

    print(i)
    tot1 += initial_photons*((r**(i+1)) * ((1-(fiber_sa))**(i+1)))
    print(tot1)
    tot2 += initial_photons*((r**(i+2)) * ((1-(fiber_sa))**(i+1)))
    print(tot2)
print((tot1+tot2)*(0.0275*(2))*0.4)

for i in range(51):
    tot += initial_photons * ((r**i * (1-fiber_sa)**i) + (r**(i+1)*(1-fiber_sa)**i))

print(tot*(0.0275*(2))*0.4)
'''