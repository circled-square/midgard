use noise::NoiseFn;
pub struct Multi<F: NoiseFn<f64, 2>> {
    f: F,
    octaves: u8,
    freq: f64,
    ampl_decay: f64,
}
impl<F: NoiseFn<f64, 2>> Multi<F> {
    pub fn new(f: F, octaves: u8, freq: f64) -> Self {
        Self {
            f,
            octaves,
            freq,
            ampl_decay: 0.5,
        }
    }
    #[allow(dead_code)]
    pub fn set_ampl_decay(&mut self, ampl_decay: f64) {
        self.ampl_decay = ampl_decay;
    }
}
impl<F: NoiseFn<f64, 2>> NoiseFn<f64, 2> for Multi<F> {
    fn get(&self, point: [f64; 2]) -> f64 {
        let mut ret = 0f64;

        let mut freq = self.freq;
        let mut ampl = 1.0;
        for _ in 0..self.octaves {
            let point = [point[0] * freq, point[1] * freq];
            ret += self.f.get(point) * ampl;

            freq *= 2.0;
            ampl *= self.ampl_decay;
        }

        let divisor = 2.0 - 0.5f64.powi(self.octaves as i32 - 1);

        ret / divisor
    }
}
