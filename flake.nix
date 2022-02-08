{
  description = "Nixified Rust project";

  inputs = {
    nixpkgs         .url = "github:nixos/nixpkgs/nixos-unstable";
    utils           .url = "github:numtide/flake-utils";
    rust-overlay = { url = "github:oxalica/rust-overlay"; inputs.nixpkgs    .follows = "nixpkgs";
                                                          inputs.flake-utils.follows = "utils"; };
    crate2nix    = { url = "github:kolloch/crate2nix";     flake = false; };
    flake-compat = { url = "github:edolstra/flake-compat"; flake = false; };
  };

  outputs = { self, nixpkgs, utils, rust-overlay, crate2nix, ... }:
    let
      # This name must match the name in Cargo.toml
      name = "nixified-rust-project";
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
                  rust-nightly = final.rust-bin.nightly."2022-02-07";
                  rust-stable  = final.rust-bin.stable ."1.58.1"    ; # nix flake lock --update-input rust-overlay
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
                    rust-analyzer-preview = rust-analyzer-preview-on "2022-02-07";
                  })
            ];
          };
          inherit (import "${crate2nix}/tools.nix" { inherit pkgs; }) generatedCargoNix;

          # Create the cargo2nix project
          project = pkgs.callPackage
            (generatedCargoNix {
              inherit name;
              src = ./.;
            })
            {
              # See "Handling external dependencies" on
              # https://ryantm.github.io/nixpkgs/languages-frameworks/rust/#rust
              defaultCrateOverrides = pkgs.defaultCrateOverrides // {
                ${name} = old-attributes: {
                  inherit buildInputs nativeBuildInputs;
                };
              };
            };

          # non-Rust dependencies
          buildInputs = [ pkgs.openssl.dev ];
          nativeBuildInputs = [ pkgs.rustc pkgs.cargo ];
        in
        rec {
          packages.${name} = project.rootCrate.build;

          # ========== nix build =========================================================
          defaultPackage = packages.${name};

          # ========== nix run ============================================================
          defaultApp = apps.${name};

          apps.${name} = utils.lib.mkApp {
            inherit name;
            drv = packages.${name};
          };

          apps.one = utils.lib.mkApp { drv = packages.${name}; name = "one"; };
          apps.two = utils.lib.mkApp { drv = packages.${name}; name = "two"; };

          # ========== nix develop ========================================================
          devShell = pkgs.mkShell {
            inputsFrom = builtins.attrValues self.packages.${system};
            buildInputs = buildInputs ++ [
              # Tools you need for development go here.
              pkgs.rust-analyzer-preview
              #pkgs.rustup.rls pkgs.rustup.rust-analysis
            ];
            RUST_SRC_PATH = "${pkgs.rustup.rust-src}/lib/rustlib/src/rust/library";
          };
        }
      );
}
