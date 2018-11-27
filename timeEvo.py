import numpy as np
import gizmo_analysis as gizmo
import wutilities as ut
import rockstar_analysis as rockstar
import matplotlib.pyplot as plt
import os

file_path = "/mnt/home/npanithanpaisal/darkpy/figs/"
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
    os.makedirs(directory)

def make_fig(pos_pa, i, pos_pa_hal, count):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1,1,1)
    ax.scatter(pos_pa[:,0], pos_pa[:,2], s=4, c='black', alpha=0.3)
    ax.scatter(pos_pa_hal[:,0], pos_pa_hal[:,2], s=50, c='red')
    ax.set_xlim((-80, 80))
    ax.set_ylim((-80, 80))
    ax.set_xlabel('x [kpc]')
    ax.set_ylabel('z [kpc]')
    ax.set_title('snapshot: {}   count: {}'.format(i, count))
    fig.savefig('/mnt/home/npanithanpaisal/darkpy/figs/xz_snap{}.png'.format(i), dpi=300)
    print('xz_snap{}.png has been created'.format(i))
    plt.close(fig)

# pos = part['star']['host.distance'][st]
# pos_halo = hal['host.distance']
def compute_dist(pos, pos_hal, hal, threshold = 2):

    print('printing first 10 halo id-------------------------------')
    print(hal['id'][:10])
    count = 0
    halo_ind = []
    for i in range(len(pos_hal)):
        # throw away those that are clearly too far away
        if np.absolute(pos_hal[i,0]) > 80 or np.absolute(pos_hal[i,1]) > 80 or np.absolute(pos_hal[i,2]) > 80:
            continue

        else:
            d_hal = np.tile(pos_hal[i], len(pos)).reshape((len(pos),3))
            d_star = pos
            d_star_hal = np.sqrt(((d_hal - d_star)**2).sum(axis=1))
            if np.min(d_star_hal) < threshold:
                count += 1
                halo_ind.append(i)

    return count, np.array(halo_ind)

def compute_dist_old(pos_pa_hal, pos_pa, threshold = 2):
    encounter = []
    dist_list = []
    for i in range(len(pos_pa_hal)):
        if np.absolute(pos_pa_hal[i,0]) > 80 or np.absolute(pos_pa_hal[i,1]) > 80 \
        or np.absolute(pos_pa_hal[i,2]) > 80:
            continue
        else:
            mindist = np.infty
            for j in range(len(pos_pa)):
                d_sq = (pos_pa_hal[i,0]-pos_pa[j,0])**2 + (pos_pa_hal[i,1]-pos_pa[j,1])**2 + \
                (pos_pa_hal[i,2]-pos_pa[j,2])**2
                if d_sq < mindist:
                    mindist = d_sq

            if mindist < threshold**2:
                encounter.append(i)
                dist_list.append(mindist)
    return np.array(encounter), np.sqrt(np.array(dist_list))


#-------------------------------------------------------------------------------

# read in merger tree
#halt = rockstar.io.IO.read_tree(simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
#hal_600 = rockstar.io.IO.get_catalog_from_tree(halt, 600)
hal_600 = rockstar.io.IO.read_catalogs('index', 600, simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
# read in stars (no dark matter) at z = 0, snapshot 600
part_600 = gizmo.io.Read.read_snapshots(['star'], 'index', 600, assign_principal_axes=True,
                                 assign_orbit=True,
                                 simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
# read in indices of stars in the stream
st = np.loadtxt('one-stream-ids.txt', dtype=int)
pos_pa600 = ut.coordinate.get_coordinates_rotated(part_600['star']['host.distance'][st],part_600.principal_axes_vectors)
pos_pa600_hal = ut.coordinate.get_coordinates_rotated(hal_600['host.distance'],part_600.principal_axes_vectors)
# compute interacting halos using principal axis coordinates
count, interacting_hal_id = compute_dist(pos_pa600, pos_pa600_hal, hal_600)
make_fig(pos_pa600, 600, pos_pa600_hal[interacting_hal_id], count)

# c is an array that keep track of how many interacting halo are there as a function of index
c = np.zeros(601)
c[600] = count
for i in range(250, 600, 1):
    try:
        # read in stars at snapshot i
        part_i = gizmo.io.Read.read_snapshots(['star'], 'index', i, assign_principal_axes=True,
                                         assign_orbit=True,
                                         simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
        hal_i = rockstar.io.IO.read_catalogs('index', i, simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')

        # read star index pointer
        gizmo.track.ParticleIndexPointer.io_pointers(part_i, directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100/track/')
        st_i = part_i.index_pointers[st]

        pos_pa_i = ut.coordinate.get_coordinates_rotated(part_i['star']['host.distance'][st_i],part_600.principal_axes_vectors)
        pos_pa_hal_i = ut.coordinate.get_coordinates_rotated(hal_i['host.distance'],part_600.principal_axes_vectors)
        count, interacting_hal_id = compute_dist(pos_pa_i, pos_pa_hal_i, hal_i)
        make_fig(pos_pa_i, i, pos_pa_hal_i[interacting_hal_id], count)

        c[i] = count
    except:
        continue
np.savetxt('/mnt/home/npanithanpaisal/darkpy/count.txt', c, header='number of interacting halos as a function of index, which is represented by row')
