![Build and Test](https://github.com/jacg/petalo/workflows/Build%20and%20Test/badge.svg)
[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)

The following instructions assume that you have installed
[Nix](https://nixos.org/) and enabled the [modern nix
interface](https://nixos.wiki/wiki/Nix_command) which includes the `nix` command
and flakes. In short, after installing Nix, you should add the line

``` conf
experimental-features = nix-command flakes
```

to one of

+ `/etc/nix/nix.conf`
+ `~/.config/nix/nix.conf`


# For users

The instructions in this section are useful in cases where you want to run the
software, but don't intend to modify it. If you want to modify the software, see
the [For developers](#for-developers) section.

To install the tools provided in this repository, as a user, you have a number
of options. For casual users, I recommend the first two.

1. Ephemeral installation with `nix shell`

   ```shell
   nix shell github:jacg/petalo
   ```

   This will place you in a shell in which (the most recent version of) all the
   tools are available. In this shell, try `mlem --help`. If the installation
   succeeded, this last command should show you the help for the main image
   reconstruction program.

   See [Available tools](#available-tools) for a list of the provided utilities.

   As soon as you exit the shell, the tools 'disappear'.

   On the first invocation, entering the shell will take a considerable amount
   of time as the software needs to be downloaded and compiled. On subsequent
   invocations, it should be instantaneous, unless

   1. A newer version is available (TODO discuss selecting revs)
   2. You have garbage collected your Nix store

2. Ad-hoc installation with `nix build`

   ```shell
   nix build github:jacg/petalo
   ```

   This will create a link called `result` in the working directory, which
   points to a `bin` directory containing all the available tools.

   After having run this command, you can try `result/bin/mlem --help`.

3. Classical package manager-like installation with `nix profile`

   This looks convenient and 'normal' to those used to non-Nix package managers,
   but swims against the current of the principles which make Nix great: it is
   therefore **NOT RECOMMENDED**!

   ```shell
   nix profile install github:jacg/petalo
   ```

   All the tools should now be available in your `PATH`: try `mlem --help`.

4. Home Manager TODO

   This is the recommended approach for the long term. Utilities to make
   home-manager installation easy, are in the works.

# For developers

## Getting the source

You can obtain the source code by cloning the repository directly via git. But
Nix flakes provide an alternative:

``` shell
nix flake clone github:jacg/petalo --dest petalorust
```

## Development tools

The Nix flake in the repository provides an environment containing the tools
needed to work on the code, which can be activated with `nix develop`:

``` shell
cd petalorust  # or wherever you cloned the repository
nix develop
```
`nix develop` places you in a shell containing all the necessary development
tools for the project. Exit the shell do disable the development environment.

### Automatic environment switching with `direnv`

I recommend using [`direnv`](https://direnv.net/) to automate provision of
development tools.

If `direnv` is enabled, the environment required to work on this project will be
enabled automatically, after you have given `direnv` permission to work in this
directory with `direnv allow`. You can revoke permission with `direnv deny`. The
first time the environment is enable, it will have to install all the required
dependencies, which may take a significant amount of time. Thereafter, it should
be instantaneous.

`direnv` will disable the environment as soon as you leave the directory. If you
want to use the environment in a different directory, you can always ask for it
explicitly with either of these commands

+ `nix develop <path-to-your-clone-of-the-project>`
+ `nix develop github:jacg/petalo

`direnv` also automatically switches the environment to match the currently
checked-out version of the project.

## Without Nix

If you haven't installed Nix ... well, it's a few orders of magnitude easier to
install Nix, than to attempt to run this software without it.

# Testing

Running `just` should compile and test all components, including Rust, Python and Julia.

TODO: Python and Julia not yet enabled in the flake

# Available tools

   + `mlem`: Reconstruction of images via Maximum Likelihood Expectation
     Maximization (MLEM) or Ordered Subsets Expectation Maximization (OSEM).

     This includes the ability to correct for scatter and attenuation.

   + `makelor`: Reconstruction of coincidence lines of response (LORs) from
     Monte Carlo simulations of detectors.

   + `make_sensitivity_image`: Generate sensitivity image (for use in `mlem`
     attenuation correction) from a density map of the field of view (FOV).

   + `viewraw.py`: Interactively view 2D slices through 3D reconstructed images
     produced by `mlem`.

   + `foms`: Calculate Figures of Merit (FOMs) for reconstructed images of
     Jaszczak or NEMA7 phantoms.

   + TODO: This on is currently missing! `foms.py`: Calculate and plot FOMs for
     a sequence of reconstructed images produced by `mlem`.

   + `imageprimaries`: Create a 3D image of the distribution of primary vertices
     (back-to-back gamma production points) in a MC simulation. Can be viewed
     with `viewraw.py`.

   + `vislor`: Visualize the interaction of an individual LOR with the FOV.
