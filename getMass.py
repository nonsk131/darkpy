import numpy as np
import gizmo_analysis as gizmo
import utilities as ut
import rockstar_analysis as rockstar
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# read in merger tree using the new merger tree
halt = rockstar.io.IO.read_tree(simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100/',
    rockstar_directory='halo/rockstar_dm_new')

i = []
for i in range(590, 601,1):
    hal_ind = np.where(halt['snapshot'] == i)[0]
    i.append(hal_ind)

i = np.array(i)
all_m = halt['mass'][i]

np.savetxt('/mnt/home/npanithanpaisal/darkpy/halos/all_masses.txt', all_m,
            header='masses of halos in snapshot 590 to 600')
