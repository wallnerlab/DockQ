//! Opt-in aggregate scorer phase timers for benchmark builds.

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

pub const PHASES: [&str; 7] = [
    "native_contacts",
    "sequence_alignment",
    "model_contacts",
    "contact_stats",
    "interface_selection",
    "geometry_rmsd",
    "total_scoring",
];

static NANOS: [AtomicU64; 7] = [
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
];

pub fn reset() {
    for counter in &NANOS {
        counter.store(0, Ordering::Relaxed);
    }
}

pub fn record(phase: usize, start: Instant) {
    let elapsed = start.elapsed().as_nanos().min(u64::MAX as u128) as u64;
    NANOS[phase].fetch_add(elapsed, Ordering::Relaxed);
}

pub fn totals() -> [u64; 7] {
    std::array::from_fn(|index| NANOS[index].load(Ordering::Relaxed))
}
