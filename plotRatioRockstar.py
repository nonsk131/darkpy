import numpy as np
import utilities as ut
import rockstar_analysis as rockstar
import matplotlib.pyplot as plt

halt = rockstar.io.IO.read_tree(simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
n_tree = np.full(601,1)
n_cata = np.full(601,1)
for i in range(250,601,1):
    print('doing index {}'.format(i))
    hal_i_tree = rockstar.io.IO.get_catalog_from_tree(halt, i)
    hal_i_cata = rockstar.io.IO.read_catalogs('index', i, simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
    n_tree[i] = len(hal_i_tree['am.phantom']) - (hal_i_tree['am.phantom'] != 0).sum()
    n_cata[i] = len(hal_i_cata['am.phantom']) - (hal_i_cata['am.phantom'] != 0).sum()


np.savetxt('/mnt/home/npanithanpaisal/darkpy/ratio.txt', np.column_stack((n_tree.astype(int), n_cata.astype(int), n_tree/n_cata)), header='halos from tree, halos from read_catalogs, ratio')
