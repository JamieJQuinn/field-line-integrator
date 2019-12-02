import numpy as np
from mayavi import mlab
import sdf
import sys
import argparse
import math

def field_twist(bx, by, bz, dx, dy, dz):
    gradbx = np.gradient(bx, dx, dy, dz)
    gradby = np.gradient(by, dx, dy, dz)
    gradbz = np.gradient(bz, dx, dy, dz)
    
    twist = (gradbz[1] - gradby[2])*bx - (gradbz[0] - gradbx[2])*by + (gradby[0] - gradbx[1])*bz
    twist = twist / (bx*bx + by*by + bz*bz)
    
    return twist

def parallel_electric_field(bx, by, bz, dx, dy, dz):
    gradbx = np.gradient(bx, dx, dy, dz)
    gradby = np.gradient(by, dx, dy, dz)
    gradbz = np.gradient(bz, dx, dy, dz)
    
    twist = (gradbz[1] - gradby[2])*bx - (gradbz[0] - gradbx[2])*by + (gradby[0] - gradbx[1])*bz
    twist = twist / (np.sqrt(bx*bx + by*by + bz*bz))
    
    return twist

def calculate_field_lines(sdfFile, n_lines, direction, integration_variable, starting_z=0):
    """Uses Mayavi to calculate streamlines. Returns Mayavi flow structure."""
    bx = sdfFile.Magnetic_Field_bx_centred.data
    by = sdfFile.Magnetic_Field_by_centred.data
    bz = sdfFile.Magnetic_Field_bz_centred.data

    dx = (sdfFile.Grid_Grid.extents[3] - sdfFile.Grid_Grid.extents[0])/sdfFile.Magnetic_Field_bx_centred.dims[0]
    dy = (sdfFile.Grid_Grid.extents[4] - sdfFile.Grid_Grid.extents[1])/sdfFile.Magnetic_Field_bx_centred.dims[1]
    dz = (sdfFile.Grid_Grid.extents[5] - sdfFile.Grid_Grid.extents[2])/sdfFile.Magnetic_Field_bx_centred.dims[2]

    print("Plotting " + integration_variable)
    if integration_variable == "twist":
        integrand = field_twist(bx, by, bz, dx, dy, dz)
    elif integration_variable == "parallel_electric_field":
        integrand = parallel_electric_field(bx, by, bz, dx, dy, dz)
    elif integration_variable == "connectivity":
        integrand = None
    else:
        integrand = None
        print(integration_variable, "not found! Purely calculating field lines.")

    # field = mlab.pipeline.vector_field(bx, by, bz)
    # magnitude = mlab.pipeline.extract_vector_norm(field)

    if integrand is not None:
        flow = mlab.flow(bx, by, bz,
                         scalars=integrand,
                         seedtype='plane', 
                         seed_resolution=n_lines)
    else:
        flow = mlab.flow(bx, by, bz,
                         seedtype='plane', 
                         seed_resolution=n_lines)

    flow.update_mode = 'non-interactive'
    flow.seed.update_mode = 'non-interactive'

    nx = sdfFile.Magnetic_Field_bx_centred.dims[0]
    nz = sdfFile.Magnetic_Field_bz_centred.dims[2]

    z_plane = float(starting_z) * nz

    # print("Resizing grid to ", nx)

    flow.seed.widget.normal = np.array([0.0, 0.0, 1.0])
    flow.seed.widget.origin = np.array([0.0, 0.0, float(z_plane) + 2.0])
    flow.seed.widget.point1 = np.array([float(nx), 0.0, float(z_plane) + 2.0])
    flow.seed.widget.point2 = np.array([0.0, float(nx), float(z_plane) + 2.0])

    flow.seed.widget.update_placement()

    flow.stream_tracer.maximum_propagation = 1000
    flow.stream_tracer.integration_direction = direction

    flow.update_streamlines = 0

    del bx, by, bz, integrand

    print("Calculated field lines")
    
    return flow

def extract_twist_and_initial_points(flow, map_point):
    """Extracts twist and starting points from streamline structure. Returns arrays of this data."""
    
    line_indices = flow.stream_tracer.get_output().lines.to_array()
    line_points = flow.stream_tracer.get_output().points.to_array()
    line_scalar = flow.stream_tracer.get_output()._get_point_data().get_array(1).to_array()

    integrated_quantities = []
    grid_points = []

    current_line_index = 0
    while current_line_index < len(line_indices):
        n_points = line_indices[current_line_index]
        mask = line_indices[current_line_index + 1 : current_line_index + 1 + n_points]
        line_length = np.sum(np.linalg.norm(line_points[mask][1:] - line_points[mask][:-1], axis=1))

        ds = np.linalg.norm(line_points[mask][1:] - line_points[mask][:-1], axis=1)
        averaged_integrand = (line_scalar[mask][1:] + line_scalar[mask][:-1])/2

        if map_point > 0:
            iz_middle = (np.abs(line_points[mask].transpose()[2,:] - map_point)).argmin()
            grid_points.append(line_points[mask][iz_middle])
        else:
            grid_points.append(line_points[mask[0]])

        integration_result = np.sum(averaged_integrand*ds)

        integrated_quantities.append(integration_result)
        current_line_index += n_points + 1
    return integrated_quantities, grid_points

def line_connection(start, end, nx, ny):
    x_min = -2
    y_min = -2
    x_max = 2
    y_max = 2

    x_width = x_max - x_min
    y_width = y_max - y_min

    start[0] = x_width*start[0]/nx + x_min
    end[0] = x_width*end[0]/nx + x_min
    start[1] = y_width*start[1]/ny + y_min
    end[1] = y_width*end[1]/ny + y_min

    start_r = math.sqrt(pow(start[0], 2) + pow(start[1], 2))
    end_r = math.sqrt(pow(end[0], 2) + pow(end[1], 2))

    boundary1 = 0.5 # boundary at r = 0.5
    boundary2 = 1.0 # boundary at r = 1.0

    if start_r <= boundary1 and end_r <=boundary1:
        # Field connects to inner region
        return 1

    if (start_r > boundary1 and start_r <= boundary2) and (end_r > boundary1 and end_r <= boundary2):
        # Field connects to outer region
        return 2

    if start_r > boundary2 and end_r > boundary2:
        # Field is connected to outer straight field
        return 4

    if (start_r <= boundary1 and (  end_r > boundary1)) or \
       (  end_r <= boundary1 and (start_r > boundary1)):
        # field is connected to inner and outer
        return 1+2

    if ((start_r > boundary1 and start_r < boundary2) and (end_r > boundary2)) or \
       ((end_r > boundary1 and end_r < boundary2) and (start_r > boundary2)):
        # field is connected to outer and straight
        return 2+4

    if (start_r <= boundary1 and (  end_r > boundary2)) or \
       (  end_r <= boundary1 and (start_r > boundary2)):
        # field is connected to inner and straight
        return 1+4

    return 0

def extract_connectivity(sdfFile, flow):
    line_indices = flow.stream_tracer.get_output().lines.to_array()
    line_points = flow.stream_tracer.get_output().points.to_array()

    nx = sdfFile.Magnetic_Field_bx_centred.dims[0]
    ny = sdfFile.Magnetic_Field_by_centred.dims[1]

    connectivities = []
    grid_points = []

    current_line_index = 0
    while current_line_index < len(line_indices):
        n_points = line_indices[current_line_index]
        mask = line_indices[current_line_index + 1 : current_line_index + 1 + n_points]

        start_pt = line_points[mask][0]
        end_pt = line_points[mask][-1]

        line_connectivity = line_connection(start_pt, end_pt, nx, ny)

        connectivities.append(line_connectivity)
        grid_points.append(line_points[mask[0]])
        current_line_index += n_points + 1
    return connectivities, grid_points

def calculate_twist_on_grid(sdfFile, n_lines, direction, map_point, integration_variable, starting_z=0):
    """Returns 2D grid of twist integrated along field lines from points on that grid"""

    if integration_variable == "connectivity":
        direction = 'forward'

    print("Calulating field lines")
    flow = calculate_field_lines(sdfFile, n_lines, direction, integration_variable, starting_z)

    print("Integrating along field lines")
    if integration_variable == "connectivity":
        twist, pts = extract_connectivity(sdfFile, flow)
    else:
        twist, pts = extract_twist_and_initial_points(flow, map_point)

    n_lines = flow.stream_tracer.get_output()._get_number_of_lines()

    del flow

    min_x = int(pts[ 0][0])
    max_x = int(pts[-1][0])
    min_y = int(pts[ 0][1])
    max_y = int(pts[-1][1])
    
    # Normalise grid back to indices
    if direction == 'both':
        n_points = int(np.sqrt(n_lines/2))
    else:
        n_points = int(np.sqrt(n_lines))

    grid_spacing = (max_x - min_x)/(n_points-1)
    
    print("Aligning to grid")
    twist_points = ((np.array(pts)[:,:2] - min_x)/grid_spacing).astype(int)
    twist_2d = np.zeros((n_points, n_points))

    # Copy twist values into 2d grid
    for i, twist_val in enumerate(twist):
        ix, iy = twist_points[i]
        twist_2d[ix, iy] += twist_val

    return twist_2d

def main():
    parser = argparse.ArgumentParser(description='Integrates a given quantity along field lines')
    parser.add_argument('--input', required=True, help='Input SDF file')
    parser.add_argument('--output', required=True, help='Output data file')
    parser.add_argument('--grid_width', required=True, help='Width of grid')
    parser.add_argument('--map_point', default=0.0, help='Option to remap twist values to centre of z-domain')
    parser.add_argument('--integration_variable', default='twist', help='Variable to integrate along field lines')
    parser.add_argument('--starting_z', default=0.0, help='Plane of initial positions of field lines')
    args = parser.parse_args()

    # Disable window popping up
    # mlab.options.offscreen = True

    print("Opening SDF file: " + args.input)
    sdfFile = sdf.read(args.input)
    twist = calculate_twist_on_grid(sdfFile, args.grid_width, 'both', float(args.map_point), args.integration_variable, args.starting_z)
    print("Saving: " + args.output)
    np.save(args.output, twist)

if __name__ == "__main__":
    main()
