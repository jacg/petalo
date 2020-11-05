![Build and Test](https://github.com/jacg/petalo/workflows/Build%20and%20Test/badge.svg)
[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)


# Installation


The following instructions assume that you have
[installed](https://github.com/jacg/IC/tree/manage-with-nix/doc/nix/install-nix)
[Nix](https://nixos.org/) and [enabled direnv](https://github.com/jacg/IC/tree/manage-with-nix/doc/nix/home-manager#automatic-environment-switching-with-direnv).

If you have not enabled direnv, you can still use this software: see the note below, at the point where you `cd petalo`.

If you haven't installed Nix ... well, it's a few orders of magnitude easier to
install Nix, than to attempt to run this software without it.

Assuming you have Nix ...

+ `git clone https://github.com/jacg/petalo.git`

+ `cd petalo`

   - if [`direnv` is
     enabled](https://github.com/jacg/IC/tree/manage-with-nix/doc/nix/home-manager#automatic-environment-switching-with-direnv),
     the environment required to run this software will be enabled
     automatically. The first time this is done, it will have to install all the
     required dependencies, which may take a significant amount of time.
     Thereafter, it should be almost instantaneous.

   - if you have not enabled `direnv` you will have to type `nix-shell` (inside
     this directory) to enable the environment.

# Running the code

**NB: the code currently fails to run on MacOS***

## Compilation

Nix has taken care of all the non-Rust dependencies, `cargo` (Rust's package
manager) will take care of compiling our code, as well as any Rust dependencies.

The first line below will have to compile a lot of Rust code (mine, and the
libraries on which it depends) and download the data file containing 1 million
LORs, which is used in the reconstructions: all this will take a few minutes the
first time around around!

## Some examples to try

+ `cargo run --bin mlem --release -- -i 6`

  Perform 6 MLEM iterations (`-i 6`) without TOF.

+ `cargo run --bin mlem --release -- -i 6 -r 20`

  Perform 6 MLEM iterations with TOF resolution of 20 picoseconds (`-r 20`).

+ `cargo run --bin mlem --release -- -i 6 -r 20 -c`

  Perform 6 MLEM iterations with TOF resolution of 20 picoseconds using the C
  implementation (`-c`) from https://github.com/jerenner/tofpet3d

+ `cargo run --bin mlem --release -- -i 6 -r 20 -n 50,50,50`

  Use a different voxel resolution (`-n 50,50,50`)

By default the resulting images will be written to the `data/out/mlem` and
`data/out/cmlem` directories.

To see all CLI options: `cargo run --bin mlem --release -- --help`

# View results in Jupyter notebook

The notebook `notebooks/view_reconstructed.ipynb` compares the images generated
by the Rust and C versions of the algorithm, over a variety of TOF resolutions
and voxel resolutions.

The notebook will look for images in the `data/out/mlem` and `data/out/cmlem`
directories, and automatically run the code to generate any images that are
missing.

I will only add stripped notebooks to the `master` branch. For now, I will try
to keep the `notebook` branch up to date with unstripped versions of notebooks
that are included in `master`. Thus, you should hopefully be able to view the
results without having to run the code yourself,
[here](https://github.com/jacg/petalo/blob/notebook/notebooks/view_reconstructed.ipynb)

# Discussion of results

Without TOF, the Rust and C result look identical, by eye.

With TOF enabled, the Rust version produces much less clear images. To be
investigated.

# Pretty things

Try the following

+ `cargo run --bin vislor -- ball`
+ `cargo run --bin vislor -- ball -r 200`
+ `cargo run --bin vislor -- ball -r 20`

  These are visualizations of the voxel weights along a LOR in a 180mm / 60
  voxel cube. The first case applies no TOF weighting, the next two apply TOF
  weigting with resolutions of 200 and 20 ps respectively.

+ `cargo run --bin vislor -- ball -f data/in/mlem/run_fastfastmc_1M_events.bin32 -e 10`

  Visualize the 10th event (`-e 10`) taken from the data file specified with
  `-f`.

  Note, currently this reads in the whole file (containing a million
  coincidences), so there is a noticeable delay.

+ `cargo run --bin vislor -- box --vbox-size 160,100,80 --nvoxels 8,5,5 --lor '0 32   -155 126 -33    151 -131 35'`

  Specify the voxel and LOR geometries to be visualized, directly on the command
  line. (When a property-based test of the voxel weights fails, the test output
  includes a CLI incantation of this sort which will visualize the failing
  case.)

## Extra

You might like to try replacing `ball` with `box` in the visualization examples,
and maybe even add a `-t 0.01`.
