import numpy as np
import gizmo_analysis as gizmo
import utilities as ut
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


#-------------------------------------------------------------------------------
# store indices of halos that we have to draw in each snapshot
indices = {k: [] for k in range(1,601,1)}
masses = []
all_i = []

# read in merger tree using the new merger tree
halt = rockstar.io.IO.read_tree(simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100/',
    rockstar_directory='halo/rockstar_dm_new')

# read in stars (no dark matter) at z = 0, snapshot 600 just to get the principal axis
part_600 = gizmo.io.Read.read_snapshots(['star'], 'index', 600, assign_host_principal_axes=True,
                                 assign_host_orbits=True,
                                 simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
# read in indices of stars in the stream
st = np.loadtxt('one-stream-ids.txt', dtype=int)

for i in range(600, 601, 1):
    #try:
    # read in stars at snapshot i
    part_i = gizmo.io.Read.read_snapshots(['star'], 'index', i, assign_host_principal_axes=True,
                                     assign_host_orbits=True,
                                     simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
    hal_i_ind = np.where(halt['snapshot'] == i)[0]

    if i < 600:
        # read star index pointer
        gizmo.track.ParticleIndexPointer.io_pointers(part_i, directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100/track/')
        st_i = part_i.index_pointers[st]
    else:
        st_i = st
    pos_pa_i = ut.coordinate.get_coordinates_rotated(part_i['star']['host.distance'][st_i],part_600.principal_axes_vectors)
    pos_pa_hal_i = ut.coordinate.get_coordinates_rotated(halt['host.distance'][hal_i_ind],part_600.principal_axes_vectors)
    count, id = compute_dist(pos_pa_i, pos_pa_hal_i)
    jj = hal_i_ind[id]
    for j in jj:
        if j not in indices[i]:
            indices[i].append(j)
        if j not in all_i:
            mass.append(halt['mass'][j])
            all_i.append(j)

    c = np.zeros(2)
    c[1] = count
    np.savetxt('/mnt/home/npanithanpaisal/darkpy/t.txt', c)

    # track the halos in N previous and N after snapshots
    N = 3
    n = 1
    prog_index = jj
    while n <= N:
        # going backward
        prog_index = halt['progenitor.main.index'][prog_index]
        prog_index = prog_index[np.where(prog_index > 0)]
        for k in prog_index:
            if k not in indices[i-n]:
                indices[i-n].append(k)
        n += 1

    if i + N > 600:
        N = 600 - i
    n = 1
    des_index = jj
    while n <= N:
        # going forward
        des_index = halt['descendant.index'][des_index]
        des_index = des_index[np.where(des_index > 0)]
        for k in des_index:
            if k not in indices[i+n]:
                indices[i+n].append(k)
        n += 1


        #make_fig(pos_pa_i, i, pos_pa_hal_i[interacting_hal_id], count)

    #except:
    #    continue
#np.savetxt('/mnt/home/npanithanpaisal/darkpy/count.txt', c, header='number of interacting halos as a function of index, which is represented by row')
