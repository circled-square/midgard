use crate::world_generator::isize_index_matrix::IsizeIndexMatrix;
use num_traits::Zero;

pub fn vec_mul_by_scalar(v: (f64, f64), s: f64) -> (f64, f64) {
    (v.0 * s, v.1 * s)
}
pub fn vec_sum(v: (f64, f64), u: (f64, f64)) -> (f64, f64) {
    (v.0 + u.0, v.1 + u.1)
}
pub fn vec_subtract(v: (f64, f64), u: (f64, f64)) -> (f64, f64){
    (v.0 - u.0, v.1 - u.1)
}
pub fn vec_dot(v: (f64, f64), u: (f64, f64)) -> f64 {
    v.0 * u.0 + v.1 + u.1
}
pub fn vec_module(v: (f64, f64)) -> f64 {
    (v.0 * v.0 + v.1 * v.1).sqrt()
}
pub fn vec_clamp(v: (f64, f64), max: f64) -> (f64, f64) {
    let module = vec_module(v);
    if module > max {
        vec_mul_by_scalar(v, max / module)
    } else {
        v
    }
}
pub fn vec_normalize(v: (f64, f64)) -> (f64, f64) {
    if v.0.is_zero() && v.1.is_zero() {
        v
    } else {
        vec_mul_by_scalar(v, 1.0 / vec_module(v))
    }
}
pub fn get_gradient(field: &Vec<Vec<f64>>, coords: (isize, isize)) -> Option<(f64, f64)> {
    Some((
        field.at_checked((coords.0 + 1, coords.1))? - field.at_checked((coords.0 - 1, coords.1))?,
        field.at_checked((coords.0, coords.1 + 1))? - field.at_checked((coords.0, coords.1 - 1))?,
    ))
}
