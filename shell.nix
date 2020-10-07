# To get more recent versions of the packages, you need to update
# `nixpkgs-commit-id`: go to https://status.nixos.org/, which lists the latest
# commit that passes all the tests for any release. Unless there is an
# overriding reason, pick the latest stable NixOS release, at the time of
# writing this is nixos-20.09.
{
 nixpkgs-commit-id ? "7badbf18c45b7490d893452beb8950d966327831" # nixos-20.09 on 2020-10-07
}:
let
  nixpkgs-url = "https://github.com/nixos/nixpkgs/archive/${nixpkgs-commit-id}.tar.gz";
  pkgs = import (fetchTarball nixpkgs-url) {
      overlays = map (uri: import (fetchTarball uri)) [
        https://github.com/mozilla/nixpkgs-mozilla/archive/master.tar.gz
      ];
    };
    rust-stable  = pkgs.latest.rustChannels.stable .rust;
    rust-nightly = pkgs.latest.rustChannels.nightly.rust;
    rust-specific = (pkgs.rustChannelOf { date = "2020-08-28"; channel = "nightly"; }).rust;
    # to use the project's rust-toolchain file:
    rust-toolchain = (pkgs.rustChannelOf { rustToolchain = ./rust-toolchain; }).rust;
    buildInputs = [
      (rust-stable.override {
        extensions = [
          "rust-analysis" # Rust Language Server (RLS)
          "rust-src"      # Needed by RLS? or only rust-analyzer?
          "rls-preview"   # What do *-preview offer?
          "rustfmt-preview"
          "clippy-preview"
        ];
      })

      pkgs.pkgconfig

      # for criterion plot generation
      pkgs.gnuplot

      # for kiss3d visualization
      pkgs.xorg.libX11
      pkgs.xorg.libXcursor
      pkgs.xorg.libXrandr
      pkgs.xorg.libXi
      pkgs.libGL
    ];
in
  pkgs.stdenv.mkDerivation {
    name = "rust";
    buildInputs = buildInputs;
    LD_LIBRARY_PATH = "${pkgs.stdenv.lib.makeLibraryPath buildInputs}";
    # RA_LOG = "info";
    # RUST_BAfCKTRACE = 1;
  }
