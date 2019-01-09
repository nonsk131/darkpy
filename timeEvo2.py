import sys
sys.path.insert(0,'/mnt/home/npanithanpaisal/packages')
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import gizmo_analysis as gizmo
import utilities as ut
import rockstar_analysis as rockstar
import os

file_path = "/mnt/home/npanithanpaisal/darkpy/figs2/"
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
    os.makedirs(directory)
file_path = "/mnt/home/npanithanpaisal/darkpy/halos/"
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
    os.makedirs(directory)

def make_fig(pos_pa, i, pos_pa_hal, pos_pa_hal2, count):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1,1,1)
    ax.scatter(pos_pa[:,0], pos_pa[:,2], s=4, c='black', alpha=0.3)
    ax.scatter(pos_pa_hal[:,0], pos_pa_hal2[:,2], s=50, c='red', alpha=0.2)
    ax.scatter(pos_pa_hal[:,0], pos_pa_hal[:,2], s=50, c='red')
    ax.set_xlim((-80, 80))
    ax.set_ylim((-80, 80))
    ax.set_xlabel('x [kpc]')
    ax.set_ylabel('z [kpc]')
    ax.set_title('snapshot: {}   count: {}'.format(i, count))
    fig.savefig('/mnt/home/npanithanpaisal/darkpy/figs2/xz_snap{}.png'.format(i), dpi=300)
    plt.close(fig)

# pos = part['star']['host.distance'][st]
# pos_halo = hal['host.distance']
def compute_dist(pos, pos_hal, threshold = 2):

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

def compute_dist_vel(pos, vel, pos_hal, vel_hal, threshold = 2, N = 10):
        count = 0
        halo_ind = []
        mindist = []
        relative = []
        for i in range(len(pos_hal)):
            # throw away those that are clearly too far away
            if np.absolute(pos_hal[i,0]) > 80 or np.absolute(pos_hal[i,1]) > 80 or np.absolute(pos_hal[i,2]) > 80:
                continue

            else:
                # compute distance for all stars in the stream
                d_hal = np.tile(pos_hal[i], len(pos)).reshape((len(pos),3))
                d_star = pos
                d_star_hal = np.sqrt(((d_hal - d_star)**2).sum(axis=1))
                if np.min(d_star_hal) < threshold:
                    count += 1
                    halo_ind.append(i)
                    mindist.append(np.min(d_star_hal))

                    d_and_v = np.column_stack((d_star_hal, vel))
                    # sort by first column
                    d_and_v = d_and_v[d_and_v[:,0].argsort()]
                    v_star_mean = np.mean(d_and_v[:N,1:], axis=0)

                    v_hal = vel_hal[i]
                    v_rel = np.sqrt(((v_hal - v_star_mean)**2).sum())
                    relative.append(v_rel)

        halo_ind = np.array(halo_ind)
        mindist = np.array(mindist)
        relative = np.array(relative)
        return count, halo_ind, mindist, relative

def compute_dist_vperp(pos, vel, pos_hal, vel_hal, threshold = 2, N = 10):
        count = 0
        halo_ind = []
        mindist = []
        relative = []
        for i in range(len(pos_hal)):
            # throw away those that are clearly too far away
            if np.absolute(pos_hal[i,0]) > 80 or np.absolute(pos_hal[i,1]) > 80 or np.absolute(pos_hal[i,2]) > 80:
                continue

            else:
                # compute distance for all stars in the stream
                d_hal = np.tile(pos_hal[i], len(pos)).reshape((len(pos),3))
                d_star = pos
                d_star_hal = np.sqrt(((d_hal - d_star)**2).sum(axis=1))
                if np.min(d_star_hal) < threshold:
                    count += 1
                    halo_ind.append(i)
                    mindist.append(np.min(d_star_hal))

                    d_and_v = np.column_stack((d_star_hal, vel))
                    # sort by first column
                    d_and_v = d_and_v[d_and_v[:,0].argsort()]
                    v_star_mean = np.mean(d_and_v[:N,1:], axis=0)
                    size = np.sqrt((v_star_mean**2).sum())
                    unit_vec = v_star_mean/size

                    v_hal = vel_hal[i]
                    v_rel = np.sqrt(((v_hal - np.dot(v_hal, unit_vec)*unit_vec)**2).sum())
                    relative.append(v_rel)

        halo_ind = np.array(halo_ind)
        mindist = np.array(mindist)
        relative = np.array(relative)
        return count, halo_ind, mindist, relative



#-------------------------------------------------------------------------------
# store indices of interacting halos in each snapshot
interacting = {k: [] for k in range(1,601,1)}
# interpolate
interpolate = {k: [] for k in range(1,601,1)}
masses = []
all_i = []
all_d = []
all_v = []
s = []
c = np.zeros(601)

# read in merger tree using the new merger tree
halt = rockstar.io.IO.read_tree(simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100/',
    rockstar_directory='halo/rockstar_dm_new')

# read in stars (no dark matter) at z = 0, snapshot 600 just to get the principal axis
part_600 = gizmo.io.Read.read_snapshots(['star'], 'index', 600, assign_host_principal_axes=True,
                                 assign_host_orbits=True,
                                 simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
# read in indices of stars in the stream
st = np.loadtxt('one-stream-ids.txt', dtype=int)

start = 375
for i in range(start, 601, 1):
    print('analyzing snapshot {}'.format(i))

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
    pos_pa_i = ut.coordinate.get_coordinates_rotated(part_i['star']['host.distance'][st_i],part_600.host_rotation_tensors[0])
    pos_pa_hal_i = ut.coordinate.get_coordinates_rotated(halt['host.distance'][hal_i_ind],part_600.host_rotation_tensors[0])
    vel_i = part_i['star']['host.velocity'][st_i]
    vel_hal_i = halt['host.velocity'][hal_i_ind]
    count, id, d, v = compute_dist_vperp(pos_pa_i, vel_i, pos_pa_hal_i, vel_hal_i)
    all_d.extend(d)
    all_v.extend(v)
    c[i] = count
    # where in the big merger tree
    jj = hal_i_ind[id]
    for j in jj:
        if j not in interacting[i]:
            interacting[i].append(j)
        if j not in all_i:
            masses.append(halt['mass'][j])
            all_i.append(j)
            s.append(i)

    # track the halos in N previous and N after snapshots
    N = 3
    n = 1
    prog_index = jj
    while n <= N:
        # going backward
        prog_index = halt['progenitor.main.index'][prog_index]
        prog_index = prog_index[np.where(prog_index > 0)]
        for k in prog_index:
            if k not in interpolate[i-n]:
                interpolate[i-n].append(k)
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
            if k not in interpolate[i+n]:
                interpolate[i+n].append(k)
        n += 1


        #make_fig(pos_pa_i, i, pos_pa_hal_i[interacting_hal_id], count)

# save snapshot, id, masses, min distance, relative velocity
id_mass = np.column_stack((np.array(s), np.array(all_i), np.array(masses), np.array(all_d), np.array(all_v)))
np.savetxt('/mnt/home/npanithanpaisal/darkpy/halos/id_mass_perp.txt', id_mass, header='save snapshot, id, masses, min distance, relative velocity')

cc = np.column_stack((np.array(range(0,601,1)), c))
np.savetxt('/mnt/home/npanithanpaisal/darkpy/halos/count.txt', cc)

for i in range(start,601,1):
    name = 'interacting_snap{}.txt'.format(i)
    path = os.path.join('/mnt/home/npanithanpaisal/darkpy/halos/', name)
    np.savetxt(path, np.array(interacting[i]))

    name = 'interpolate_snap{}.txt'.format(i)
    path = os.path.join('/mnt/home/npanithanpaisal/darkpy/halos/', name)
    np.savetxt(path, np.array(interpolate[i]))

    part_i = gizmo.io.Read.read_snapshots(['star'], 'index', i, assign_host_principal_axes=True,
                                     assign_host_orbits=True,
                                     simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')
    if i < 600:
        # read star index pointer
        gizmo.track.ParticleIndexPointer.io_pointers(part_i, directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100/track/')
        st_i = part_i.index_pointers[st]
    else:
        st_i = st
    pos_pa_i = ut.coordinate.get_coordinates_rotated(part_i['star']['host.distance'][st_i],part_600.host_rotation_tensors[0])
    pos_pa_hal_i = ut.coordinate.get_coordinates_rotated(halt['host.distance'][interacting[i]],part_600.host_rotation_tensors[0])
    pos_pa_hal_i2 = ut.coordinate.get_coordinates_rotated(halt['host.distance'][interpolate[i]],part_600.host_rotation_tensors[0])
    make_fig(pos_pa_i, i, pos_pa_hal_i, pos_pa_hal_i2, c[i])
