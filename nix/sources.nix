# To get a more recent version of nixpkgs, go to https://status.nixos.org/,
# which lists the latest commit that passes all the tests for any release.
# Unless there is an overriding reason, pick the latest stable NixOS release, at
# the time of writing this is nixos-21.05.

{ pkgs }:

let
  random_pkgs = import <nixpkgs> {};
  nixpkgs-commit-id = "4f37689c8a219a9d756c5ff38525ad09349f422f"; # nixos-21.05 on 2021-11-28
  nixpkgs-url = "https://github.com/nixos/nixpkgs/archive/${nixpkgs-commit-id}.tar.gz";
  oxalica-commit-id = "d9a664513558376595e838b21348cdac0ba3115e"; # 2021-11-28
  pkgs = import (fetchTarball nixpkgs-url) {
      overlays = map (uri: import (fetchTarball uri)) [
        "https://github.com/oxalica/rust-overlay/archive/${oxalica-commit-id}.tar.gz"
      ];
    };
in
{

  ####### Pinned nixpkgs ##################################################

  pkgs = pkgs;

  ####### home-manager ####################################################

  home-manager = let
    src = builtins.fetchGit {
      name = "home-manager-release-21.05-2021-08-03";
      url = https://github.com/nix-community/home-manager;
      rev = "b39647e52ed3c0b989e9d5c965e598ae4c38d7ef";
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
