{
  description = "Image Reconstruction for PET";

  inputs = {
    nixpkgs         .url = "github:nixos/nixpkgs/nixos-unstable"; # cargo2nix broken in 21.11
    utils           .url = "github:numtide/flake-utils";
    rust-overlay = { url = "github:oxalica/rust-overlay"; inputs.nixpkgs    .follows = "nixpkgs";
                                                          inputs.flake-utils.follows = "utils"; };
    crate2nix    = { url = "github:kolloch/crate2nix";     flake = false; };
    flake-compat = { url = "github:edolstra/flake-compat"; flake = false; };
  };

  outputs = { self, nixpkgs, utils, rust-overlay, crate2nix, ... }:
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
              rust-overlay.overlay (final: prev:
                let
                  # If you have a rust-toolchain file for rustup, choose `rustup =
                  # rust-tcfile` further down to get the customized toolchain
                  # derivation.
                  rust-tcfile  = final.rust-bin.fromRustupToolchainFile ./rust-toolchain;
                  rust-latest  = final.rust-bin.stable .latest      ;
                  rust-beta    = final.rust-bin.beta   .latest      ;
                  rust-nightly = final.rust-bin.nightly."2022-04-05";
                  rust-stable  = final.rust-bin.stable ."1.59.0"    ; # nix flake lock --update-input rust-overlay
                  rust-analyzer-preview-on = date:
                    final.rust-bin.nightly.${date}.default.override
                      { extensions = [ "rust-analyzer-preview" ]; };
                in
                  rec {
                    # The version of the Rust system to be used in buildInputs. Choose between
                    # tcfile/latest/beta/nightly/stable (see above) on the next line
                    rustup = rust-stable;

                    rustc = rustup.default;
                    cargo = rustup.default;
                    rust-analyzer-preview = rust-analyzer-preview-on "2022-04-05";
                  })
              # ==== Cargo nextest ========================================================
              (final: prev: {
                cargo-nextest = final.callPackage ./overlays/cargo-nextest.nix {};
              })
            ];
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

          # ----- Python -------------------------------------------------------------------
          python-version = "python39";

          my-python-packages = pypkgs: [
            pypkgs.numpy
            pypkgs.matplotlib
            pypkgs.pytest
            pypkgs.h5py
            pypkgs.scikit-learn
            pypkgs.docopt
          ];


          inherit (import "${crate2nix}/tools.nix" { inherit pkgs; }) generatedCargoNix;

          # Create the cargo2nix project
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
                  HDF5_DIR = pkgs.symlinkJoin { name = "hdf5"; paths = [ pkgs.hdf5 pkgs.hdf5.dev ]; };
                };
              };
            };

          # non-Rust dependencies
          buildInputs = [ pkgs.hdf5 ];
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

          apps.mlem                   = utils.lib.mkApp { drv = packages."${name}-all"; name = "mlem"                  ; };
          apps.makelor                = utils.lib.mkApp { drv = packages."${name}-all"; name = "makelor"               ; };
          apps.foms                   = utils.lib.mkApp { drv = packages."${name}-all"; name = "foms"                  ; };
          apps.imageprimaries         = utils.lib.mkApp { drv = packages."${name}-all"; name = "imageprimaries"        ; };
          apps.joinlorhdf             = utils.lib.mkApp { drv = packages."${name}-all"; name = "joinlorhdf"            ; };
          apps.vislor                 = utils.lib.mkApp { drv = packages."${name}-all"; name = "vislor"                ; }; # TODO X11 missing at runtime
          apps.make_sensitivity_image = utils.lib.mkApp { drv = packages."${name}-all"; name = "make_sensitivity_image"; };
          apps.viewraw                = utils.lib.mkApp { drv = packages."${name}-all"; name = "viewraw.py"; };

          # ========== nix develop ========================================================
          devShell = pkgs.mkShell {
            inputsFrom = builtins.attrValues self.packages.${system};
            buildInputs = buildInputs ++ [
              # Tools you need for development go here.
              pkgs.just
              pkgs.rust-analyzer-preview
              pkgs.cargo-nextest
              #pkgs.rustup.rls pkgs.rustup.rust-analysis
            ];
            RUST_SRC_PATH = "${pkgs.rustup.rust-src}/lib/rustlib/src/rust/library";
            HDF5_DIR = pkgs.symlinkJoin { name = "hdf5"; paths = [ pkgs.hdf5 pkgs.hdf5.dev ]; };
          };
        }
      );
}
