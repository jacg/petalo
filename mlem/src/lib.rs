#[no_mangle]
pub extern "C" fn sum_floats(x: *const f32, n: usize) -> f32 {
    let slice = unsafe { std::slice::from_raw_parts(x, n) };
    let vec = slice.to_vec();
    vec.iter().sum()
}
