# To get a more recent version of nixpkgs, go to https://status.nixos.org/,
# which lists the latest commit that passes all the tests for any release.
# Unless there is an overriding reason, pick the latest stable NixOS release, at
# the time of writing this is nixos-20.09.

{ pkgs }:

let
  random_pkgs = import <nixpkgs> {};
  nixpkgs-commit-id = "95d26c9a9f2a102e25cf318a648de44537f42e09"; # nixos-20.09 on 2020-10-24
  nixpkgs-url = "https://github.com/nixos/nixpkgs/archive/${nixpkgs-commit-id}.tar.gz";
  pkgs = import (fetchTarball nixpkgs-url) {
      overlays = map (uri: import (fetchTarball uri)) [
        https://github.com/mozilla/nixpkgs-mozilla/archive/master.tar.gz
      ];
    };
in
{
  pkgs = pkgs;
}
