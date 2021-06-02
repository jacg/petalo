# Quick and dirty instructions

(Assuming you are on Linux: MacOS (as always) is more complicated.)

+ cd into this directory
+ `cargo build --release`
+ `cd ..`
+ `ln -s target/release/libfulano.so libfulano.so`
+ `python`
+ `import fulano`
