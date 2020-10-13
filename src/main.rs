use petalo::weights::{Point, VoxelBox};
use petalo::visualize;

fn main() {

    let p1 = Point::new(-265.1371093069, 76.0, -200.0);
    let p2 = Point::new( 276.002,         0.0,  200.0);
    let vbox = VoxelBox::new((120.0, 100.0, 80.0), (50, 40, 30));

    visualize::lor_weights(10.00, 10.03, p1, p2, vbox);
}
