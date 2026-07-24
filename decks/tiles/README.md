# Tile benchmarks

## 2d thermal plasma

The input deck os-stdin-2d-therm was used to test the speed of the tiles code compared to OSIRIS without tiles.  The listed variables were varied in the following ways:
- Number of processors in _x_ and _y_, varied between 2 and 16 by factors of 2.
- Number of threads, varied between 64 and 2 by factors of 2.
- Number of tiles in _x_ and _y_, varied between 2 and 128 by factors of 2.
- Max number of particles per tile, set to 1.5 times the expected amount for a single tile.

Optimal tile size was found to be 16x16 cells with 144 particles per cell, and between 32x32 and 64x64 cells for 1 ppc.

## 3d thermal plasma

The input deck os-stdin-3d-therm was used to test the speed of the tiles code compared to OSIRIS without tiles.  The listed variables were varied in the following ways:
- Number of processors in _x_, _y_, and _z_, varied between 2 and 8 by factors of 2.
- Number of threads, varied between 64 and 2 by factors of 2.
- Number of tiles in _x_, _y_, and _z_, varied between 2 and 16 by factors of 2.
- Max number of particles per tile, set to 1.5 times the expected amount for a single tile.

Optimal tile size was found to be 8x8x8 cells with 64 particles per cell, and 16x16x16 cells for 1 ppc.
