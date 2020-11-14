# Bootstrap your Rust development environment

## TL;DR:

Ensure that

+ [Nix](https://nixos.org/) is [enabled](../install-nix/README.md) for your user
  account
+ you have a network connection which is not absurdly slow

Execute:

```shell
sh petalo/doc/nix/home-manager/install-home-manager-and-maybe-spacemacs.sh
```

This will install `home-manager` and use it to provide you with

+ Automatic direnv caching: this dramatically reduces the delay you experience
  when switching to a nix environment, such as the one used to develop this
  project.

+ An installation of [Spacemacs](https://www.spacemacs.org/) (tagline: The best
  editor is neither Emacs nor Vim, it's Emacs *and* Vim!), pre-configured with
  Rust development tools such as
  [rust-analyzer](https://rust-analyzer.github.io/) and [Rust Language
  Server](https://rls.booyaa.wtf/), as well as more general tools such as Magit.

## Details

+ By default this will create a directory `$HOME/my-home-manager` into which it
  will place your new `home-manager` configuration files. If you would like to
  use a different location, run the installation script thus:

  ```shell
  HM_DIR=<location of your choice> sh petalo/doc/nix/home-manager/install-home-manager-and-maybe-spacemacs.sh
  ```

+ If you already have `home-manager` and only want to install Spacemacs for Rust
  development ... I'm working on it: TODO.

+ By default, Spacemacs will be configured to use Emacs editing style. If you
  prefer vim, then run the installation script thus:

  ```shell
  SPACEMACS=vim sh petalo/doc/nix/home-manager/install-home-manager-and-maybe-spacemacs.sh
  ```
  You can override both the `SPACEMACS` and `HM_DIR` defaults when running the script: `SPACEMACS=vim HM_DIR=wherever sh ...`.

+ If you only want `home-manager` to be installed, without Spacemacs:

  ```shell
  SPACEMACS=no sh petalo/doc/nix/home-manager/install-home-manager-and-maybe-spacemacs.sh
  ```

+ If you are already using Emacs and use `$HOME/.emacs.d` for its configuration,
  `home-manager` will move that configuration to `$HOME/.emacs.d.bck`.

  Eventually I aim to provide an option to install Spacemacs in a way that does
  not interfere with a pre-existing Emacs installation: TODO.

# Working with `home-manager`

## TL;DR:

To install (or remove) software from your environment

1. Edit `$HM_DIR/nixpkgs.home.nix`, adding (or removing) items in the
   `home.packages` list.

2. Instruct home manager to switch to the new specification, with `home-manager
   switch -b bck switch`.

The `-b bck` option ensures that any pre-existing files in locations to which
`home-manager` tries to write, will not be deleted, but moved aside to
`<original name>.bck`.

## What is home-manager?

Home manager uses Nix to manage the installation and configuration of software
that you use in your per-user working environment. Benefits include:

+ Install and uninstall software without admin privileges.

+ Try new packages with guarantees that you can revert to the previous state
  without anything breaking.

+ Configurations are reliable and reproducible.

+ Transferring your personal environment to a different machine is trivially
  easy (as long as Nix is available on that machine).

  This means that **installing your environment on *ANY* nix-enabled machine is as
  simple as launching your script, and letting home-manager download and
  configure all the packages you are used to.**

## Ephemeral package installation with `nix-shell`

Nix makes it very easy to try out software packages without fear of breaking
anything in your existing environment. (With other package managers, installing
and uninstalling some package often leaves the system in a modified state.)

One particularly convenient way of trying a package in Nix is with ephemeral
installations using `nix-shell`.

As an example, let's take `lolcat` for a spin. Just like `cowsay` and `figlet`,
`lolcat` is a package which isn't terribly useful, but makes a good guinea-pig
for installation/uninstallation experiments. Install `lolcat` ephemerally and
try it out, with:

```shell
nix-shell -p lolcat # This will drop you in a shell where lolcat is available
ls -l / | lolcat    # A pretty, colourful listing should appear
```

Within this shell (and *only* this shell), you can play around with `lolcat` to
your heart's content: it adds pretty colours to whatever you pipe into it.

When you have seen enough, simply exit this shell, and `lolcat` disappears.
While you were in the shell, the rest of your system was unaware of the
existence of `lolcat`: it could not have interfered with anything else. In the
case of trivial packages like `lolcat` this isn't terribly important, but in the
case of packages on which other software might depend, this can be crucial.

The first time you run `nix-shell -p lolcat` it is likely to take some time as
it might need to download `lolcat` and its dependencies. On subsequent
invocations, it is likely to run much faster, as everything that is needed has
been kept in the nix store on your machine ... at least until you instruct Nix
to collect garbage. You may want to run `nix-collect-garbage` if the Nix store
grows too large and you want to free up disk space wasted on packages (or older
versions) you are not using.

## Persistent package installation with `home-manager`

In the previous section we installed `lolcat` *ephemerally* with `nix-shell`:
the installed package was only available in a single shell, and disappeared as
soon as the shell was exited.

It is possible to install and remove packages *persistently* in your personal
environment by using `nix-env` to *mutate* your personal profile (you can think
of this as a personal `brew`, `apt-get`, `pacman`, `yum`, etc.). I do **NOT**
recommend this.

I suggest you use `home-manager` to manage your personal configuration
*declaratively*. The declarative approach makes it *much* easier to transfer
your environment to different machines, to understand what is present in your
environment, to share package configuration wisdom with your colleagues, to get
help from others, and much more.

If you have tried a package by installing it ephemerally with `nix-shell` and
found that it is useful to you, then you might want to make it persistently
available. With `home-manager` this is a two-step process:

+ Add the package to the `home.packages` list in your `home.nix`

+ Instruct `home-manager` to create and switch to an environment matching your
  new specification:

  ```shell
  home-manager switch -b bck
  ``**

## Version control your environments

**NOTE**: in all that follows, `$HM_DIR` refers to the location you chose for
your `home-manager` configuration directory. By default, this is
`$HOME/my-home-manager`, but you might have overridden this by prepending
`HM_DIR=<some place of your choice>` to the execution of the installaiton.

I recommend that you place your environment specification in version control:
```shell
git init $HM_DIR
```

The `home.nix` that is installed by the installation script, contains the line

```nix
  home.file.".config/nixpkgs".source = link ../nixpkgs;
```

This instructs `home-manager` to create a link that informs your system where
your `home-manager` configuration can be found. This is the simplest short-term
way of getting home-manager to work for you. However, I would recommend a more
powerful and all-encompassing approach for the long-term:

1. Place *all your personal configuration files/directories*

   - `$HOME/.bash{rc,_profile}`
   - `$HOME/.emacs.d`
   - `$HOME/.config/htop`
   - etc.

   into `$HM_DIR`

2. Instruct home-manager to make these files appear in their standard locations.
   See `$HM_DIR/{eg-ro,eg-rw,nixpkgs/home.nix}` for examples of this.

   (If you allowed the installaiton script to install Spacemacs for you, then it
   has already done this for `.emacs.d` and `.spacemacs`.)

3. Keep track of changes with Git (or any VCS of your choice).

In this scheme, installing your personal environment on *any* Nix-enabled
machine amounts to no more than

1. `git clone <your-home-manager-repo> $HM_DIR`
2. `cd $HM_DIR`
3. `nix-shell bootstrap-home-manager`
4. `mkdir -p $HOME/.config/nixpkgs`
5. `ln -s $HM_DIR/nixpkgs/home.nix $HOME/.config/nixpkgs/home.nix`
6. DONE!

Note that all these steps can be wrapped into a single script, so **you could
install your whole environment with a single command!**

Caveat: when you take this approach, you must be careful not to commit any
*secrets* (e.g. private ssh keys) into a repository which you will ever place in
some publicly accessible locations (such as GitHub).
