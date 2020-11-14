![Build and Test](https://github.com/jacg/petalo/workflows/Build%20and%20Test/badge.svg)
[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)


# Installation


The following instructions assume that you have
[installed](doc/nix/install-nix/README.md)
[Nix](https://nixos.org/) and enabled [`direnv`](https://direnv.net/).

If you haven't installed Nix ... well, it's a few orders of magnitude easier to
install Nix, than to attempt to run this software without it.

If you have not enabled `direnv`, you can still use this software: see the note
below, at the point where you `cd petalo`. If you want to enable direnv, This
[automated procedure](doc/nix/home-manager/README.md) will (among other things!)
do it for you.

Assuming you have Nix ...

+ `git clone https://github.com/jacg/petalo.git`

+ `cd petalo`

   - if `direnv` is enabled, the environment required to run this software will
     be enabled automatically, after you have given `direnv` permission to work
     in this directory direnv `with allow`. The first time this is done, it will
     have to install all the required dependencies, which may take a significant
     amount of time. Thereafter, it should be almost instantaneous.

   - if you have not enabled `direnv` you will have to type `nix-shell` (inside
     this directory) to enable the environment.

# Running the code

## Input data

You will need an input data file containing a table with the following 21
columns of float-64s:

```
event_id, true_energy,
true_r1, true_phi1, true_z1, true_t1,
true_r2, true_phi2, true_z2, true_t2,
phot_like1, phot_like2,
reco_r1, reco_phi1, reco_z1, reco_t1,
reco_r2, reco_phi2, reco_z2, reco_t2,
not_sel
```

Tell the program where this file is located with the `-f` (`--input-file`)
option. Oterwise it will look for `data/in/full_body_phantom_reco_combined.h5`

The program assumes that the table is in a dataset called `reco_info/table`. You
can override this with the `-d` (`--dataset`) CLI option.

By default, the reco coordinates are used, but you can opt for the true ones
with `--use-true`.

## Compilation

Nix has taken care of all the non-Rust dependencies, `cargo` (Rust's package
manager) will take care of compiling our code, as well as any Rust dependencies.

The first line below will have to compile a lot of Rust code (mine, and the
libraries on which it depends): all this will take a few minutes the first time
around around!

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

# Results

There are clearly some differences between the TOF runs on Rust vs C.

# Pretty things

[If your machine is not running NixOS, Nix may have trouble adapting to your
graphics hardware, and these visualizations might not work: TODO!]

Try the following

+ `cargo run --bin vislor -- ball`
+ `cargo run --bin vislor -- ball -r 200`
+ `cargo run --bin vislor -- ball -r 20`

  These are visualizations of the voxel weights along a LOR in a 180mm / 60
  voxel cube. The first case applies no TOF weighting, the next two apply TOF
  weigting with resolutions of 200 and 20 ps respectively.

+ `cargo run --bin vislor -- ball -f data/in/full_body_phantom_reco_combined.h5 -e 10`

  Visualize the 10th event (`-e 10`) taken from the data file specified with `-f`.

+ `cargo run --bin vislor -- box --vbox-size 160,100,80 --nvoxels 8,5,5 --lor '0 32   -155 126 -33    151 -131 35'`

  Specify the voxel and LOR geometries to be visualized, directly on the command
  line. (When a property-based test of the voxel weights fails, the test output
  includes a CLI incantation of this sort which will visualize the failing
  case.)

## Extra

You might like to try replacing `ball` with `box` in the visualization examples,
and maybe even add a `-t 0.01`.
