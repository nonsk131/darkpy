import numpy as np
import wutilities as ut
import rockstar_analysis as rockstar
import matplotlib.pyplot as plt

halt = rockstar.io.IO.read_tree(simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
n_tree = np.zeros(601)
n_cata = np.zeros(601)
n_tree[0] = 1
n_cata[0] = 1
for i in range(1,601,1):
    hal_i_tree = rockstar.io.IO.get_catalog_from_tree(halt, i)
    hal_i_cata = rockstar.io.IO.read_catalogs('index', i, simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
    n_tree[i] = len(hal_i_tree['host.distance'])
    n_cata[i] = len(hal_i_cata['host.distance'])

np.savetxt('/mnt/home/npanithanpaisal/darkpy/ratio.txt', np.column_stack((n_tree.astype(int), n_cata.astype(int), n_tree/n_cata)), header='halos from tree, halos from read_catalogs, ratio')
