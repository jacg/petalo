use std::env;
use std::path::PathBuf;

fn main() {

    // Tell cargo to tell rustc to link the tofpet3d libmlem that was installed
    // by shell.nix
    println!("cargo:rustc-link-lib=mlem");

    // Tell cargo to invalidate the built crate whenever the headers change
    println!("cargo:rerun-if-changed=wrapper.h");

    let bindings = bindgen::Builder::default()
        // The header for which we want to generate bindings
        .header("wrapper.h")
        // Invalidate the built crate if headers change
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate bindings
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings");
}
