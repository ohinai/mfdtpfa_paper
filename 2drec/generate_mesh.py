
import random
import subprocess

import numpy as np 

count = 0
scale = [4, 8, 16, 32, 64, 128, 256]

for n in scale:
    points_file = open("mesh"+str(n), 'w')
    count = 0
    current_list = []
    for j in range(0, n):
        for i in range(0, n):
            new_x = 1./n*i + .5*1./float(n)
            new_x += 3./50.*abs(np.sin(4.*np.pi*new_x))
            new_y = 1./n*j + .5*1./float(n)
            new_y += 3./50.*abs(np.sin(4.*np.pi*new_y))
            print >>points_file, count, new_x, new_y, .5
                ##new_point = np.array([new_x, new_y, new_z])
                ##current_list.append(new_point)
            count += 1

    points_file.close()
    subprocess.call("voro++ -c \"%i$%q$%v$%C$%P$%t$%f$%l$%n\" 0 1 0 1 0 1 " + "mesh"+str(n), shell=True)

