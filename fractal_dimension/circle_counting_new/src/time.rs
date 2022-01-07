#![allow(dead_code)]
use std::time::{Duration, Instant};

/// Handy struct to keep track of timing things and making predictions
pub struct Timer {
    start: Instant,
    active: bool,
}

impl Timer {
    /// Constructor returning a timer. `active` tells whether or not to enable message printing.
    pub fn new(active: bool) -> Self {
        Self {
            start: Instant::now(),
            active,
        }
    }

    /// Returns the time elapsed since the timer was constructed
    pub fn time_elapsed(&self) -> Duration {
        Instant::now().duration_since(self.start)
    }

    /// Given a duration, prints `before` followed by the duration, in minutes/seconds, followed by
    /// `after`.
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

    /// Prints `before` + time elapsed since timer was activated + `after`
    pub fn time_stamp(&self, before: &str, after: &str) {
        self.print_time_message(self.time_elapsed(), before, after);
    }

    /// Exactly like `time_stamp` by prints `factor * self.time_elapsed()` rather than just
    /// `self.time_elapsed()`
    pub fn time_remaining(&self, factor: f64, before: &str, after: &str) {
        self.print_time_message(
            Duration::from_secs_f64(self.time_elapsed().as_secs_f64() * factor),
            before,
            after,
        );
    }
}
