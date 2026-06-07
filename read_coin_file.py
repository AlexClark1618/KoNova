import time
import numpy as np
import matplotlib.pyplot as plt

file_path = 'KNVA-20260514-01-00079_coinc.dat'

coincidence_list = []

with open(file_path, 'r') as file:

    while True:

        header = file.readline()

        if not header:   # EOF
            break

        n1, n2 = map(int, header.split())
        n = n1 + n2

        block = []

        for _ in range(n):

            ts, qdc, ch = file.readline().split()

            block.append([
                int(ts),
                float(qdc),
                int(ch)
            ])

        coincidence_list.append(block)

print(f"Total number of coincidence events: {len(coincidence_list)}")
print(f'Coincidence rate: {len(coincidence_list) / 100} Hz')

bar_ch_map = 'bar_ch_map.csv'
with open(bar_ch_map, 'r') as f:
    next(f) #Skip header
    lines = f.readlines()

    bar_to_ch = {}

    for line in lines:
        bar, ch = line.strip().split(',')
        bar_to_ch[int(bar)] = int(ch)

offsets = {1: 64, 2: 320, 3: 576, 4: 832} #Layer-ch num offset dict

ch_to_bar = {}

for bar, base_ch in bar_to_ch.items():

    for layer, offset in offsets.items():

        scaled_ch = base_ch + offset
        ch_to_bar[scaled_ch] = (layer, bar)

layer_freq_dict = {
    layer: {bar: 0 for bar in range(1, 65)}
    for layer in range(1, 5)
}

layer_qdc_dict = {
    layer: []
    for layer in range(1, 5)
}

for block in coincidence_list:
    #print(block)

    for event in block:

        #print(event)
        ch = event[2]

        if ch in ch_to_bar:

            layer, bar = ch_to_bar[ch]

            layer_freq_dict[layer][bar] += 1

            layer_qdc_dict[layer].append(event[1])

#print(layer_freq_dict)
'''
        if 64 <= event[2] <=127:
            print("Port 1")
            rescaled_ch = event[2] - 64
            layer_freq_dict[1].append(event)

        elif 256 <= event[2] <=319:
            print("Port 3")

        elif 576 <= event[2] <=639:
            print("Port 5")

        elif 832 <= event[2] <=895:
            print("Port 7")
   '''         
#print(layer_qdc_dict.values())
#print(layer_freq_dict[1].keys())

coincidence_array = np.array(coincidence_list, dtype=object)

new_coincidence_array = []

for event in coincidence_array.flatten():
    block = []
    for hit in event:
        #print(hit[2])
        if hit[2] in ch_to_bar:
            layer, bar = ch_to_bar[hit[2]]
            #print(f"Layer: {layer}, Bar: {bar}, QDC: {hit[1]}, Timestamp: {hit[0]}")
            block.append([layer, bar, hit[1], hit[0]])
    
    new_coincidence_array.append(block)

new_coincidence_array = np.array(new_coincidence_array,  dtype=object)

#print(new_coincidence_array[0])
#for hit in new_coincidence_array[0]:

bar_center_to_center_sep = 16.5 #mm
hits_xy = []
for event in new_coincidence_array:
    #Bottom x-y
    layer_1_hit_list = np.array([hit for hit in event if hit[0]==1],  dtype=object)
    #print(layer_1_hit_list)

    if len(layer_1_hit_list) == 0:
        continue
    #for hit in layer_1_hit_list:
    x2 = np.mean(layer_1_hit_list[:,1])*16.5
    #print(np.mean(layer_1_hit_list[:,1])*16.5)

    layer_2_hit_list = np.array([hit for hit in event if hit[0]==2],  dtype=object)
    #print(layer_2_hit_list)

    if len(layer_2_hit_list) == 0:
        continue
    #for hit in layer_1_hit_list:
    y2 = np.mean(layer_2_hit_list[:,1])*16.5
    #print(np.mean(layer_2_hit_list[:,1])*16.5)

    #Top x-y
    layer_3_hit_list = np.array([hit for hit in event if hit[0]==3],  dtype=object)
    #print(layer_3_hit_list)
    if len(layer_3_hit_list) == 0:
        continue
    #for hit in layer_1_hit_list:
    x1 = np.mean(layer_3_hit_list[:,1])*16.5

    #print(np.mean(layer_3_hit_list[:,1])*16.5)

    layer_4_hit_list = np.array([hit for hit in event if hit[0]==4],  dtype=object)
    #print(layer_4_hit_list)
    if len(layer_4_hit_list) == 0:
        continue
    #for hit in layer_1_hit_list:
    y1 = np.mean(layer_4_hit_list[:,1])*16.5

    #print(np.mean(layer_4_hit_list[:,1])*16.5)

    hits_xy.append([(x1, y1), (x2, y2)])

#print(hits_xy)

import numpy as np

delta_z =  11.75 * 2.54 *10 # in mm from x layer on top to x layer on bottom.  Same for y.  Measured 4/29/2026 see Photos
def zenith_angle(p1, p2):
    """
    Calculate the zenith angle between two x,y points relative to the z-axis.
    
    The zenith angle is the angle between the vector connecting p1→p2
    and the positive z-axis [0, 0, 1].
    
    Parameters
    ----------
    p1, p2 : array-like, shape (2,)
        Points as (x, y) coordinates.
    
    Returns
    -------
    float
        Zenith angle in degrees.
    """
    # Build the 2D displacement vector, then lift it into 3D (z=0 plane)
    v = np.array([p2[0] - p1[0], p2[1] - p1[1], delta_z])
    #print(v)
    z_axis = np.array([0.0, 0.0, 1.0])

    # cos θ = (v · z) / (|v| |z|)  — |z| = 1, so it simplifies
    magnitude = np.linalg.norm(v)
    #print(magnitude)
    if magnitude == 0:
        raise ValueError("Points are identical — zenith angle is undefined.")

    cos_theta = np.dot(v, z_axis) / magnitude   # dot product with z-hat = v[2] / |v|
    zenith_deg = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

    return zenith_deg

import numpy as np

def azimuth_angle(p1, p2):
    """
    Calculate the azimuth angle of the vector p1→p2 with respect to the x-axis.
    
    Parameters
    ----------
    p1, p2 : array-like, shape (2,)
        Points as (x, y) coordinates.
    
    Returns
    -------
    float
        Azimuth angle in degrees, in range [0, 360).
    """
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]

    # atan2 gives angle in (-180, 180] measured CCW from +x axis
    angle_rad = np.arctan2(dy, dx)
    angle_deg = np.degrees(angle_rad)

    # Normalise to [0, 360)
    return angle_deg % 360

#print(hits_xy[0])
zenith_dist = []
azimuth_dist = []
for hit in hits_xy:
    p1, p2 = hit
    #print(p1, p2)
    zenith = zenith_angle(p1, p2)
    azimuth = azimuth_angle(p1, p2)

    #print(f"Zenith angle: {zenith:.2f}°")
    #print(f"Azimuth angle: {azimuth:.2f}°")
    zenith_dist.append(zenith)
    azimuth_dist.append(azimuth)
# → 90.00° (vector lies flat in the x,y plane, perpendicular to z)

plt.hist(zenith_dist, bins=int(np.sqrt(len(zenith_dist))))
plt.title(f'Zenith Angle Distribution \n Mean: {np.mean(zenith_dist):.2f}° ')

plt.show()

plt.hist(azimuth_dist, bins=int(np.sqrt(len(azimuth_dist))))
plt.title('Azimuth Angle Distribution')
plt.show()

for layer in layer_freq_dict:
    plt.bar(layer_freq_dict[layer].keys(), layer_freq_dict[layer].values())
    plt.title(f'Layer {layer} Bar Hit Frequency')
    plt.xlabel('Bar Number')
    plt.ylabel('Frequency')
    #plt.show()

    layer_qdc_mean = np.mean(layer_qdc_dict[layer])
    plt.hist(layer_qdc_dict[layer], bins=100)
    plt.title(f'Layer {layer} QDC Distribution | Mean: {layer_qdc_mean:.2f}')
    plt.xlabel('QDC Value')
    plt.ylabel('Frequency')
    #plt.show()