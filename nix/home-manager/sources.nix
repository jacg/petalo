# To get a more recent version of nixpkgs, go to https://status.nixos.org/,
# which lists the latest commit that passes all the tests for any release.
# Unless there is an overriding reason, pick the latest stable NixOS release, at
# the time of writing this is nixos-20.09.

{ pkgs }:

let
  random_pkgs = import <nixpkgs> {};
  nixpkgs-commit-id = "896270d629efd47d14972e96f4fbb79fc9f45c80"; # nixos-20.09 on 2020-11-11
  nixpkgs-url = "https://github.com/nixos/nixpkgs/archive/${nixpkgs-commit-id}.tar.gz";
  pkgs = import (fetchTarball nixpkgs-url) {
      overlays = map (uri: import (fetchTarball uri)) [
        https://github.com/mozilla/nixpkgs-mozilla/archive/master.tar.gz
      ];
    };
in
{

  ####### Pinned nixpkgs ##################################################

  pkgs = pkgs;

  ####### home-manager ####################################################

  home-manager = let
    src = builtins.fetchGit {
      name = "home-manager-2020-11-06";
      url = https://github.com/nix-community/home-manager;
      rev = "4cc1b77c3fc4f4b3bc61921dda72663eea962fa3";
    };
  # `path` is required for `home-manager` to find its own sources
  in pkgs.callPackage "${src}/home-manager" { path = "${src}"; };

  ####### SPACEMACS ########################################################

  # Bumping the version gave all sorts of problems with Emacs packages not being able to be installed.
  # completely deleting the .emacs.d/.cache and .emacs.d/elpa directories seems to have fixed the problem.
  # After this, the first emacs installation (which installs a ton of stuff) complains a number of times like this:
  #
  #    <package-name>.<date>.<time>.el is write-protedcet; try to save anyway? (y or n)
  #
  # answering with `y` seems to get around these problems ... until the next version bump
  spacemacs = {
    # Don't make the directory read-only so that packages can be installed and
    # caches written.
    recursive = true;
    source = pkgs.fetchFromGitHub {
      owner = "syl20bnr";
      repo = "spacemacs";
      rev = "c3872f165c3ea0862cdb939c5e7b7494b5ce0e72"; # develop branch on 2020-11-09
      #sha256 = lib.fakeSha256;
      sha256 = "1dp1rffyqrnxd31m1x1ldja24db9wml117915wa9h7ixinvqjdfv";
    };
  };
}
