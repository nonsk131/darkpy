import numpy as np
import gizmo_analysis as gizmo
import wutilities as ut
import rockstar_analysis

# read in stars (no dark matter) at z = 0
part600 = gizmo.io.Read.read_snapshots(['star'], 'redshift', 0, assign_principal_axes=True,
                                 assign_orbit=True,
                                 simulation_directory='/mnt/ceph/users/firesims/fire2/metaldiff/m12i_res7100')

print('end of script')
