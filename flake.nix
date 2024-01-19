{
  description = "Image Reconstruction for PET";

  inputs = {
    nixpkgs         .url = "github:nixos/nixpkgs/nixos-23.11";
    oldpkgs         .url = "github:nixos/nixpkgs/nixos-22.11";
    utils           .url = "github:numtide/flake-utils";
    rust-overlay = { url = "github:oxalica/rust-overlay"; inputs.nixpkgs    .follows = "nixpkgs";
                                                          inputs.flake-utils.follows = "utils"; };
    crate2nix    = { url = "github:kolloch/crate2nix";     flake = false; };
    flake-compat = { url = "github:edolstra/flake-compat"; flake = false; };
  };

  outputs = { self, nixpkgs, oldpkgs, utils, rust-overlay, crate2nix, ... }:
    let
      # This name must match the name in Cargo.toml
      name = "petalo";
    in
    utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [
              # ===== Specification of the rust toolchain to be used ====================
              rust-overlay.overlays.default (final: prev:
                let
                  # If you have a rust-toolchain file for rustup, choose `rustup =
                  # rust-tcfile` further down to get the customized toolchain
                  # derivation.
                  rust-tcfile  = final.rust-bin.fromRustupToolchainFile ./rust-toolchain;
                  rust-latest  = final.rust-bin.stable .latest      ;
                  rust-beta    = final.rust-bin.beta   .latest      ;
                  rust-nightly = final.rust-bin.nightly."2024-01-16";
                  rust-stable  = final.rust-bin.stable ."1.75.0"    ; # nix flake lock --update-input rust-overlay
                  rust-analyzer-preview-on = date:
                    final.rust-bin.nightly.${date}.default.override
                      { extensions = [ "rust-analyzer-preview" ]; };
                in
                  rec {
                    # The version of the Rust system to be used in buildInputs. Choose between
                    # tcfile/latest/beta/nightly/stable (see above) on the next line
                    rustup = rust-stable;

                    rustc = rustup.default;
                    #cargo = rustup.default; # overriding cargo causes problems on 23.11, but we don't needed it?
                    rust-analyzer-preview = rust-analyzer-preview-on "2024-01-16";
                  })
            ];
          };

          # The Rust HDF5 crate doesn't support HDF 1.14.0 yet, which is what comes with nixpkgs 23.05.
          # nixpkgs 23.11 has HDF5 1.14.3; use in the rust crate is blocked by https://github.com/aldanor/hdf5-rust/pull/243
          old = import oldpkgs { inherit system; };

          # X11 support
          libPath = pkgs.lib.makeLibraryPath [
            pkgs.libGL
            pkgs.libxkbcommon
            (linux pkgs.wayland)
            pkgs.xorg.libX11
            pkgs.xorg.libXcursor
            pkgs.xorg.libXi
            pkgs.xorg.libXrandr
          ];

          # ----- Conditional inclusion ----------------------------------------------------
          nothing = pkgs.coreutils;
          linux      = drvn: if pkgs.stdenv.isLinux  then drvn else nothing;
          darwin     = drvn: if pkgs.stdenv.isDarwin then drvn else nothing;
          linuxList  = list: if pkgs.stdenv.isLinux  then list else [];
          darwinList = list: if pkgs.stdenv.isDarwin then list else [];
          darwinStr  = str:  if pkgs.stdenv.isDarwin then str  else "";

          # ----- Darwin-specific ----------------------------------------------------------
          darwin-frameworks = pkgs.darwin.apple_sdk.frameworks;

          # ----- Python -------------------------------------------------------------------
          python-version = "python311";

          my-python-packages = pypkgs: [
            pypkgs.ipython
            pypkgs.numpy
            pypkgs.matplotlib
            pypkgs.pytest
            pypkgs.h5py
            pypkgs.scikit-learn
            pypkgs.docopt
          ];


          inherit (import "${crate2nix}/tools.nix" { inherit pkgs; }) generatedCargoNix;

          # Create the crate2nix project
          project = import
            (generatedCargoNix {
              inherit name;
              src = ./.;
            })
            {
              inherit pkgs;
              # See "Handling external dependencies" on
              # https://ryantm.github.io/nixpkgs/languages-frameworks/rust/#rust
              defaultCrateOverrides = pkgs.defaultCrateOverrides // {

                petalo = old-attributes: {
                  buildInputs = [ (darwin darwin-frameworks.AppKit) ];
                };

                hdf5-sys = old-attributes: {
                  HDF5_DIR = pkgs.symlinkJoin { name = "hdf5"; paths = [ old.hdf5 old.hdf5.dev ]; };
                };
              };
            };

          # non-Rust dependencies
          buildInputs = [ old.hdf5 ];
          nativeBuildInputs = [ pkgs.rustc pkgs.cargo ];
        in
        rec {
          packages."${name}-rust" = project.workspaceMembers.petalo.build;

          packages."${name}-python" = pkgs.${python-version}.pkgs.buildPythonApplication {
            pname = "${name}-python";
            version = "TODO-version";
            src = ./.;
            propagatedBuildInputs = [ (my-python-packages pkgs.${python-version}.pkgs) ];
          };

          packages."${name}-all" = pkgs.buildEnv {
            name = "${name}-all";
            paths = [
              packages."${name}-rust"
              packages."${name}-python"
            ];
          };

          # ========== nix build =========================================================
          defaultPackage = packages."${name}-all";

          # ========== nix run ============================================================
          defaultApp = apps.mlem;

          apps.${name} = utils.lib.mkApp {
            inherit name;
            drv = packages.${name};
          };

          # TODO: adding and removing items from the following list seems to
          # have no effect on what ends up in the build's `bin` directory.
          # Regardless of what is added or removed, we always end up with the
          # following (although with different hashes):
          #
          # result
          # ├── bin
          # │  ├── foms -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/foms
          # │  ├── imageprimaries -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/imageprimaries
          # │  ├── joinlorhdf -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/joinlorhdf
          # │  ├── make_sensitivity_image -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/make_sensitivity_image
          # │  ├── makelor -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/makelor
          # │  ├── mlem -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/mlem
          # │  ├── show_lorogram -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/show_lorogram
          # │  ├── viewraw.py -> /nix/store/q2zbc0vrpikc1i79isv4ai4a6cffwcmh-petalo-python-TODO-version/bin/viewraw.py
          # │  ├── vislor -> /nix/store/9vay8qps8xwslqymmlrv7657nwxf4km8-rust_petalo-0.1.0/bin/vislor
          # │  └── xenon_thickness_from_h5.py -> /nix/store/q2zbc0vrpikc1i79isv4ai4a6cffwcmh-petalo-python-TODO-version/bin/xenon_thickness_from_h5.py
          # └── lib -> /nix/store/q2zbc0vrpikc1i79isv4ai4a6cffwcmh-petalo-python-TODO-version/lib

          apps.mlem                   = utils.lib.mkApp { drv = packages."${name}-all"; name = "mlem"                  ; };
          apps.makelor                = utils.lib.mkApp { drv = packages."${name}-all"; name = "makelor"               ; };
          apps.foms                   = utils.lib.mkApp { drv = packages."${name}-all"; name = "foms"                  ; };
          apps.foms-all               = utils.lib.mkApp { drv = packages."${name}-all"; name = "foms.py"               ; };
          apps.imageprimaries         = utils.lib.mkApp { drv = packages."${name}-all"; name = "imageprimaries"        ; };
          apps.joinlorhdf             = utils.lib.mkApp { drv = packages."${name}-all"; name = "joinlorhdf"            ; };
          apps.vislor                 = utils.lib.mkApp { drv = packages."${name}-all"; name = "vislor"                ; }; # TODO X11 missing at runtime
          apps.make_sensitivity_image = utils.lib.mkApp { drv = packages."${name}-all"; name = "make_sensitivity_image"; };
          apps.viewraw                = utils.lib.mkApp { drv = packages."${name}-all"; name = "viewraw.py"            ; };

          # ========== nix develop ========================================================
          devShell = pkgs.mkShell {
            inputsFrom = builtins.attrValues self.packages.${system};
            buildInputs = buildInputs ++ [
              # Tools you need for development go here.
              pkgs.just
              pkgs.rust-analyzer-preview
              pkgs.cargo-nextest
              #pkgs.rustup.rls pkgs.rustup.rust-analysis
              pkgs.bacon
            ];
            PIP_DISABLE_PIP_VERSION_CHECK = 1;
            RUST_SRC_PATH = "${pkgs.rustup.rust-src}/lib/rustlib/src/rust/library";
            HDF5_DIR = pkgs.symlinkJoin { name = "hdf5"; paths = [ old.hdf5 old.hdf5.dev ]; };
            LD_LIBRARY_PATH = libPath;
          };
        }
      );
}
