use noise::NoiseFn;

pub struct FallOffNoiseFunction<T>{
    source: T,
    center : f64,
    radius : f64,
}

impl<T> FallOffNoiseFunction<T> {
    pub fn new(source : T, center : f64, radius : f64) -> Self {
        Self { source, center, radius}
    }
    pub fn falloff_offset_function(point: [f64;2], center : f64, radius : f64) -> f64 {
        let dist_x = point[0] - center;
        let dist_y = point[1] - center;
        let dist = (dist_x*dist_x + dist_y*dist_y).sqrt();
        let dist_coeff = dist / radius;

        if dist_coeff < 0.9 {
            0.
        } else {
            (dist_coeff-0.9) * 10.0
        }
    }
}

impl<T : NoiseFn<f64, 2>> NoiseFn<f64, 2> for FallOffNoiseFunction<T> {
    fn get(&self, point: [f64; 2]) -> f64 {
        let val = self.source.get(point);

        let offset = Self::falloff_offset_function(point, self.center, self.radius);

        val - offset
    }
}
