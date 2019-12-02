# Field Line Integrator

This is a small tool designed to take SDF files of magnetic fields produced with lare3d and integrate certain quanities along the magnetic field lines.

The lines are calculated using the tool [Mayavi](https://docs.enthought.com/mayavi/mayavi/) and the given quantity integrated along with a simple midpoint rule integration. This produces a grid of integrated values.

Currently the quantities available are:

- Parallel electric field - $(\vec{j} \cdot \vec{B})/|\vec{B}|$
- Twist - $(\vec{j} \cdot \vec{B})/|\vec{B}|^2$
- "Connectivity" - a field line is labeled with a certain connectivity if it starts and ends in specifically defined boundary zones. This is very situation dependent.

This tool is purely provided for interest/inspiration and is unlikely to be of direct use in other research. If you would like to improve it, please do submit a pull request!

## Usage

Running `python twist_calc.py` will display a brief help. Descriptions of some options are below.

`grid_width` sets the width of the square grid of seed points.

`map_point` attempts to remap the grid of integrated quanities back to the given point. This often doesn't work very well.

`integration_variable` sets which variable should be integrated

`starting_z` sets the z-location of the initial plane of seed points.

## Example

The following example will integrate the parallel electric field along field lines derived from the field in `input.sdf` producing a grid of 700x700 outputs in `output.npy`:

```
python twist_calc.py --integration_variable parallel_electric_field --input input.sdf --output output.npy --grid_width 700
```

## Limitations

- Will *not* work with the default staggered grid magnetic field output from lare. The field must be centred.
- Grid of field line seed points is always normal to $z$ and square. This can be changed by editing the script.

## Visualisation of outputs

Also included is a very basic visualisation tool for quick plotting. For production plots, I recommend importing the outputs into python and plotting them properly.
