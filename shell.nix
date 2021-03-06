# To get more recent versions of the packages, you need to update
# `nixpkgs-commit-id` in `nix/sources.nix`.
{
  py ? "38" # To override the default python version:  nix-shell shell.nix --argstr py 37
}:
let
  random_pkgs = import <nixpkgs> {};
  sources = random_pkgs.callPackage nix/sources.nix {};
  pkgs = sources.pkgs;

  # ----- Rust --------------------------------------------------------------------
  rust-stable  = pkgs.latest.rustChannels.stable .rust;
  rust-nightly = pkgs.latest.rustChannels.nightly.rust;
  rust-specific = (pkgs.rustChannelOf { date = "2020-10-19"; channel = "nightly"; }).rust;
  rust-pinned-stable = (pkgs.rustChannelOf { channel = "1.47.0"; }).rust;
  # to use the project's rust-toolchain file:
  rust-toolchain = (pkgs.rustChannelOf { rustToolchain = ./rust-toolchain; }).rust;
  # Rust system to be used in buldiInputs. Choose between
  # stable/nightly/specific on the next line
  rust = (rust-pinned-stable.override {
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
      substituteInPlace makefile --replace "CXX = g++" "CXX=$CXX"
      substituteInPlace cc/mlem.cc --replace "float ToFFunction"  "extern \"C\" float ToFFunction"
      source tofpet3d_setup.sh
      make
    '';

    installPhase = ''
      mkdir -p $out/lib
      mv lib/libmlem.so $out/lib/libmlem.so
    '';

    fixupPhase = darwinStr ''
      install_name_tool -id $out/lib/libmlem.so $out/lib/libmlem.so
    '';
  };

  # ----- Conditional inclusion ----------------------------------------------------
  nothing = pkgs.coreutils;
  linux      = drvn: if pkgs.stdenv.isLinux  then drvn else nothing;
  darwin     = drvn: if pkgs.stdenv.isDarwin then drvn else nothing;
  linuxList  = list: if pkgs.stdenv.isLinux  then list else [];
  darwinList = list: if pkgs.stdenv.isDarwin then list else [];
  darwinStr  = str:  if pkgs.stdenv.isDarwin then str  else "";

  # ----- Darwin-specific ----------------------------------------------------------
  darwin-frameworks = pkgs.darwin.apple_sdk.frameworks;

  # ----- Getting OpenGL to work on non-NixOS --------------------------------------
  nixGL = pkgs.callPackage "${builtins.fetchTarball {
    url = https://github.com/guibou/nixGL/archive/7d6bc1b21316bab6cf4a6520c2639a11c25a220e.tar.gz;
    #sha256 = pkgs.lib.fakeSha256;
    sha256 = "02y38zmdplk7a9ihsxvnrzhhv7324mmf5g8hmxqizaid5k5ydpr3";
  }}/nixGL.nix" {};
  # nixGL = pkgs.callPackage "${builtins.fetchTarball {
  #       url = https://github.com/guibou/nixGL/archive/17c1ec63b969472555514533569004e5f31a921f.tar.gz;
  #       sha256 = "0yh8zq746djazjvlspgyy1hvppaynbqrdqpgk447iygkpkp3f5qr";
  #     }}/nixGL.nix" {};
  # --------------------------------------------------------------------------------
  buildInputs = [
    rust

    pkgs.pkgconfig

    # for kiss3d visualization
    pkgs.xorg.libX11
    pkgs.xorg.libXcursor
    pkgs.xorg.libXrandr
    pkgs.xorg.libXi
    pkgs.libGL
    # For graphics hardware matching on non-NixOS
    (linux nixGL.nixGLDefault)

    # HDF5
    pkgs.hdf5

    # python
    pypkgs.cffi
    pypkgs.jupyter
    pypkgs.matplotlib

    # Needed for compilation to succeed on Macs
    (darwin darwin-frameworks.AppKit)
    (darwin darwin-frameworks.CoreText)

    # C version of mlem
    tofpet3d

    # Used by the Jupyter notebook stripping git filte
    pkgs.jq

    # Hacking around bindgen/dyld/MacOS: a compiler that can be used in
    # buildPhase on both Mac/Linux
    pkgs.clang_9

    ######## Tools for profiling, debugging, benchmarking #######

    # Profiling
    (linux pkgs.linuxPackages.perf)
    (linux pkgs.oprofile)
    (linux pkgs.kcachegrind)
    (linux pkgs.graphviz) # used by kcachegrind
    pkgs.cargo-flamegraph
    pkgs.cargo-asm
    pkgs.valgrind

    # Benchmarking
    pkgs.hyperfine
    pkgs.gnuplot # for criterion plot generation

    # Debugging
    pkgs.gdbgui
    pkgs.gdb
    pkgs.lldb

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
