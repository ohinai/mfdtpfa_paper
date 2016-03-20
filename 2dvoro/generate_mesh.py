
import random
import subprocess

import numpy as np 

count = 0
scale = [4, 8, 16, 32, 64, 128]

for n in scale:
    points_file = open("mesh"+str(n), 'w')
    count = 0
    current_list = []
    while count < n*n:
        new_x = random.random()
        new_y = random.random()
        if .9999 >new_x > 1.e-5 and .9999 > new_y > 1.e-5:
            print >>points_file, count, new_x, new_y, .5
                ##new_point = np.array([new_x, new_y, new_z])
                ##current_list.append(new_point)
            count += 1
    print count
    points_file.close()
    subprocess.call("voro++ -c \"%i$%q$%v$%C$%P$%t$%f$%l$%n\" 0 1 0 1 0 1 " + "mesh"+str(n), shell=True)

