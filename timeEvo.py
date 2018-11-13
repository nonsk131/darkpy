import numpy as np
import gizmo_analysis as gizmo
import wutilities as ut
import rockstar_analysis
import matplotlib.pyplot as plt
import os

file_path = "/mnt/home/npanithanpaisal/darkpy/figs/"
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
    os.makedirs(directory)

def make_fig(pos_pa, i):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1,1,1)
    ax.scatter(pos_pa[:,0], pos_pa[:,2], s=4, c='black', alpha=0.3)
    ax.set_xlim((-80, 80))
    ax.set_ylim((-80, 80))
    ax.set_xlabel('x [kpc]')
    ax.set_ylabel('z [kpc]')
    fig.savefig('/mnt/home/npanithanpaisal/darkpy/figs/xz_snap{}.png'.format(i), dpi=300)
    print('xz_snap{}.png has been created'.format(i))
    close()

# read in stars (no dark matter) at z = 0, snapshot 600
part_600 = gizmo.io.Read.read_snapshots(['star'], 'index', 600, assign_principal_axes=True,
                                 assign_orbit=True,
                                 simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
# read in indices of stars in the stream
st = np.loadtxt('one-stream-ids.txt', dtype=int)
pos_pa600 = ut.coordinate.get_coordinates_rotated(part_600['star']['host.distance'][st],part_600.principal_axes_vectors)
make_fig(pos_pa600, 600)

for i in range(250, 600, 1):
    try:
        # read in stars at snapshot i
        part_i = gizmo.io.Read.read_snapshots(['star'], 'index', i, assign_principal_axes=True,
                                         assign_orbit=True,
                                         simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')

        # read star index pointer
        gizmo.track.ParticleIndexPointer.io_pointers(part_i, directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100/track/')
        st_i = part_i.index_pointers[st]

        # check that all the indices are not null
        if np.isnan(st_i).sum() > 0:
            print(i)

        pos_pa_i = ut.coordinate.get_coordinates_rotated(part_i['star']['host.distance'][st_i],part_i.principal_axes_vectors)
        make_fig(pos_pa_i, i)
    except:
        continue
