#![allow(dead_code)]
use std::time::{Duration, Instant};

pub struct Timer {
    start: Instant,
    active: bool,
}

impl Timer {
    pub fn new(active: bool) -> Self {
        Self {
            start: Instant::now(),
            active,
        }
    }

    pub fn time_elapsed(&self) -> Duration {
        Instant::now().duration_since(self.start)
    }

    fn print_time_message(&self, duration: Duration, before: &str, after: &str) {
        if self.active {
            let time = duration.as_secs_f64();
            let seconds = time % 60.0;
            let minutes = (time as isize) / 60;

            if minutes > 0 {
                println!("{}{}m {}s{}", before, minutes, seconds, after);
            } else {
                println!("{}{}s{}", before, seconds, after);
            }
        }
    }

    pub fn time_stamp(&self, before: &str, after: &str) {
        self.print_time_message(self.time_elapsed(), before, after);
    }

    pub fn time_remaining(&self, portion: f64, before: &str, after: &str) {
        self.print_time_message(
            Duration::from_secs_f64(self.time_elapsed().as_secs_f64() * portion),
            before,
            after,
        );
    }
}
