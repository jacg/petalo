use petalo::weights::{Point, VoxelBox};
use petalo::visualize;

fn main() {

    let p1 = Point::new(242.73334, 0.0, 95.23425);
    let p2 = Point::new(-236.53258, -54.5143, -94.40064);
    let vbox = VoxelBox::new((115.8324, 100.0, 100.0), (5, 5, 5));

    visualize::lor_weights(10.00, 10.03, p1, p2, vbox);
}
