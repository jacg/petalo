use petalo::weights::{Point, VoxelBox};
use petalo::visualize;

fn main() {

    let p1 = Point::new(-265.1371093069, 76.0, 0.0);
    let p2 = Point::new( 276.002,         0.0, 0.0);
    let vbox = VoxelBox::new((120.0, 100.0), (5, 5));

    visualize::lor_weights(p1, p2, vbox);
}
