use std::time::SystemTime;

pub struct PerformanceProfiler {
    start_time: SystemTime,
    last_time: SystemTime,
}

impl PerformanceProfiler {
    pub fn new(start_time: SystemTime) -> Self {
        Self {
            start_time,
            last_time: start_time,
        }
    }

    #[allow(dead_code)]
    pub fn get_current_time_in_ms(&self) -> SystemTime {
        let current_time = SystemTime::now();
        return current_time;
    }

    pub fn get_elapsed_time_in_ms(&mut self, from: SystemTime) -> u128 {
        let current_time = SystemTime::now();
        let difference = current_time.duration_since(from).expect("Clock error");
        let time_in_ms = difference.as_millis();
        self.last_time = current_time;
        return time_in_ms;
    }

    pub fn print_elapsed_time_in_ms(&mut self, text: &str) {
        println!(
            "Profiler: {text} {} ms",
            self.get_elapsed_time_in_ms(self.last_time)
        );
    }

    pub fn print_total_elapsed_time_in_ms(&mut self, text: &str) {
        println!(
            "Profiler: {text} {} ms",
            self.get_elapsed_time_in_ms(self.start_time)
        );
    }
}
