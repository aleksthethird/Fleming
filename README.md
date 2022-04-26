
# Fleming Survey

## About

The Fleming survey aims to survey short period variable stars ($P<8$ hours) along
the Galactic plane using the James Gregory Telescope at the University of St Andrews, Fife.

This is the Python pipeline used to reduce the data and automatically select variable stars.

- [DR1: Scholz, Warwick & van Aalten (2019)](https://iopscience.iop.org/article/10.3847/2515-5172/ac32bc)


## Credits

- Concept: Aleks Scholz
- Pipeline: Thomas van Aalten, Callum McGregor
- DR1: Ben Warwick, Aleks Scholz
- DR2: Callum McGregor, Aleks Scholz
- Observations: James Gregory Telescope, University of St Andrews

## Usage

See [`docs/MANUAL.md`](docs/MANUAL.md).


## Requirements

Uses Python `>= 3.8.3`

Packages:

- `numpy`
- `matplotlib`
- `ipython`
- `astropy`
- `astroplan`
- `parse`
