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

# Running some code

I suggest you try this:

```sh
cargo run --bin mlem --release -- -i 6 -f mlem_output/tof_OFF
cargo run --bin mlem --release -- -i 6 -f mlem_output/tof_200 -r 200
cargo run --bin mlem --release -- -i 6 -f mlem_output/tof_20  -r 20
cargo run --bin mlem --release -- -i 6 -f mlem_output/tof_50  -r 50
jupyter notebook notebooks/view_reconstructed.ipynb
```

# What does this do?

The first line will have to compile a lot of Rust code (mine, and the libraries
on which it depends) and download the data file containing 1 million LORs, which
is used in the reconstructions: all this will take a few minutes the first time
around around!

Each of the first 4 lines performs 6 MLEM iterations (`-i 6`), with different
TOF resolutions (`-r <ps>`: No TOF, 200 ps, 20 ps, 50 ps).

The last line uses a Jupyter notebook to display all the generated images.

(You can view the notebook without running any code [here](https://github.com/jacg/petalo/blob/master/notebooks/view_reconstructed.ipynb))

# Discussion of results

According to the
[documentation](https://tofpet3d.readthedocs.io/en/latest/src/code_reconstruction.html)
the default TOF sensitivity in `MLEMReconstructor` is 200ps, which I am told is
to be interpreted as the sigma of the TOF gaussian. I don't understand how this
can have any significant effect on the resulting image, as the whole image
occupies a 180mm cube, while 200ps is around 70mm. The results that are shown
the notebook seem to agree: by eye I cannot distinguish the no-TOF and 200ps
case, while the 20ps case is *much* cleaner.

How does this compare to your results?

# Pretty things

Try the following

+ `cargo run --bin vislor -- ball`
+ `cargo run --bin vislor -- ball -r 200`
+ `cargo run --bin vislor -- ball -r 20`

This is a visualization of the voxel weights along a LOR in a 180mm / 60 voxel
cube. The first case applies no TOF weighting, the next two apply TOF weigting
with resolutions of 200 and 20 ps respectively.

To me, once again, the no-TOF and 200ps cases are indistinguishable by eye.

## Extra

You might like to try replacing `ball` with `box` in the visualization examples,
and maybe even add a `-t 0.01`.
