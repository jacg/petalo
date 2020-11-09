# To get more recent versions of the packages, you need to update
# `nixpkgs-commit-id`: go to https://status.nixos.org/, which lists the latest
# commit that passes all the tests for any release. Unless there is an
# overriding reason, pick the latest stable NixOS release, at the time of
# writing this is nixos-20.09.
{ nixpkgs-commit-id ? "95d26c9a9f2a102e25cf318a648de44537f42e09" # nixos-20.09 on 2020-10-24
, py ? "38" # To override the default python version:  nix-shell shell.nix --argstr py 37
}:
let
  nixpkgs-url = "https://github.com/nixos/nixpkgs/archive/${nixpkgs-commit-id}.tar.gz";
  pkgs = import (fetchTarball nixpkgs-url) {
      overlays = map (uri: import (fetchTarball uri)) [
        https://github.com/mozilla/nixpkgs-mozilla/archive/master.tar.gz
      ];
    };

  # ----- Rust --------------------------------------------------------------------
  rust-stable  = pkgs.latest.rustChannels.stable .rust;
  rust-nightly = pkgs.latest.rustChannels.nightly.rust;
  rust-specific = (pkgs.rustChannelOf { date = "2020-10-19"; channel = "nightly"; }).rust;
  # to use the project's rust-toolchain file:
  rust-toolchain = (pkgs.rustChannelOf { rustToolchain = ./rust-toolchain; }).rust;
  # Rust system to be used in buldiInputs. Choose between
  # stable/nightly/specific on the next line
  rust = (rust-stable.override {
    extensions = [
      "rust-analysis" # Rust Language Server (RLS)
      "rust-src"      # Needed by RLS? or only rust-analyzer?
      "rls-preview"   # What do *-preview offer?
      "rustfmt-preview"
      "clippy-preview"
    ];
  });

  # ----- Python -------------------------------------------------------------------
  python = builtins.getAttr ("python" + py) pkgs;
  pypkgs = python.pkgs;

  # ----- tofpet3d original C version ----------------------------------------------
  tofpet3d = pkgs.stdenv.mkDerivation {
    name = "tofpet3d";
    version = "2020-05-06";
    src = pkgs.fetchFromGitHub {
      owner = "jerenner";
      repo = "tofpet3d";
      rev= "a8c3fc8293fd05547f9ca752abc259173bac57af";
      sha256= "0zimqadzva4p69ipanpxkzd9q12nm6q3qq9f5qiraxzy38px6rln";
    };

    buildPhase = ''
      mkdir -p $out/lib
    	${pkgs.clang_9}/bin/clang++ -fpic -c -O2 -pg cc/mlem.cc        -o mlem.o
	    ${pkgs.clang_9}/bin/clang++ -shared -o $out/lib/libmlemORIG.so    mlem.o
  '';

    installPhase = ''
      echo "Nothing to be done here ... already done in buildPhase"
    '';
  };

  # ----- Conditional inclusion ----------------------------------------------------
  nothing = pkgs.coreutils;
  linux      = drvn: if pkgs.stdenv.isLinux  then drvn else nothing;
  darwin     = drvn: if pkgs.stdenv.isDarwin then drvn else nothing;
  linuxList  = list: if pkgs.stdenv.isLinux  then list else [];
  darwinList = list: if pkgs.stdenv.isDarwin then list else [];

  # ----- Darwin-specific ----------------------------------------------------------
  darwin-frameworks = pkgs.darwin.apple_sdk.frameworks;

  # --------------------------------------------------------------------------------
  buildInputs = [
    rust

    pkgs.pkgconfig

    # for criterion plot generation
    pkgs.gnuplot

    # for kiss3d visualization
    pkgs.xorg.libX11
    pkgs.xorg.libXcursor
    pkgs.xorg.libXrandr
    pkgs.xorg.libXi
    pkgs.libGL

    # for plotters
    pkgs.cmake
    pkgs.freetype
    pkgs.expat

    # HDF5
    pkgs.hdf5

    # profiling
    (linux pkgs.linuxPackages.perf)
    pkgs.cargo-flamegraph

    # python
    pypkgs.cffi
    pypkgs.jupyter
    pypkgs.matplotlib

    # Needed for compilation to succeed on Macs
    (darwin darwin-frameworks.AppKit)
    (darwin darwin-frameworks.CoreText)

    # Downoad LOR data file
    pkgs.wget

    # C version of mlem
    tofpet3d

    # Used by the Jupyter notebook stripping git filte
    pkgs.jq

    # Hacking around bindgen/dyld/MacOS: a compiler that can be used in
    # buildPhase on both Mac/Linux
    pkgs.clang_9

  ];

in

pkgs.stdenv.mkDerivation {
  name = "petalorust";
  buildInputs = buildInputs;
  LD_LIBRARY_PATH = "${pkgs.stdenv.lib.makeLibraryPath buildInputs}";
  HDF5_DIR = pkgs.symlinkJoin { name = "hdf5"; paths = [ pkgs.hdf5 pkgs.hdf5.dev ]; };

  # Needed if using bindgen to wrap C libraries in Rust
  LIBCLANG_PATH = "${pkgs.llvmPackages_9.libclang}/lib";

  # RA_LOG = "info";
  # RUST_BAfCKTRACE = 1;
}
