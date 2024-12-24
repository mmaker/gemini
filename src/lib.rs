//! Gemini: elastic arguments for R1CS.
//!
//! This library provides essentially two arguments:
//! - [`snark::Proof`], for non-preprocessing SNARKs.
//!     It provides a non-interactive succinct argument of knowledge for R1CS
//!     without indexer, and where the verifier complexity is linear in the circuit size.
//! - [`psnark::Proof`] for preprocessing SNARKs.
//!     It provides a non-interactive succinct argument of knowledge for R1CS
//!     where the verifier complexity is logarithmic in the circuit size.
//!
//! Choice of the pairing-friendly elliptic curve,
//! is entirely up to the implementor.
//! All arguments are internally using the [`kzg`](crate::kzg) commitment scheme.
//! Support for generic univariate or multivariate commitments will is scheduled and will
//! happen at some point in the future.
//!
//! Both arguments rely on some sub-protocols, implemented as separate modules in [`subprotocols`]
//! and free of use for other protocols.
//!
//! # Building
//!
//! This package can be compiled with `cargo build`, and requires rust nightly at least
//! until [`Iterator::advance_by`] hits stable. There's a bunch of feature flags that can be turned
//! on:
//!
//! - `asm`, to turn on the assembly backend within [`ark-ff`](https://docs.rs/ark-ff/);
//! - `parallel`, to turn on multi-threading. This requires the additional dependency [`rayon`](https://docs.rs/rayon/latest/rayon/);
//! - `std`, to rely on the Rust Standard library;
//! - `print-trace`, to print additional information concerning the execution time of the sub-protocols. **This feature must be enabled if you want to print the execution time of the examples.**
//!
//! # Benchmarking
//!
//! Micro-benchmarks aare available and can be fired with:
//!
//! ```bash
//! cargo bench
//! ```
//!
//! Execution of (preprocessing-)SNARK for arbitrary instance sizes can be done running
//! the examples with:
//!
//! ```bash
//! cargo run --example snark -- -i <INSTANCE_LOGSIZE>
//! ```
//!
//! # License
//!
//! This package is licensed under MIT license.
//!

#![feature(
    iter_advance_by,
    iter_next_chunk,
    slice_as_chunks,
    array_chunks
)]
#![no_std]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use)]
// #![forbid(unsafe_code)]
#![deny(trivial_numeric_casts)]
#![deny(unused_allocation)]

// Lints disable from other arkworks packages:
// - trivial_casts: this causes errors for "trivial casts" when converting an Iterable `It` to
// `dyn Iterable`.

// #![deny(stable_features, , non_shorthand_field_patterns)]
// #![deny(unused_attributes, unused_imports, unused_mut, unused_import_braces, unused_qualifications, missing_docs)]
// #![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#[macro_use]
extern crate ark_std;

/// The domain separator, used when proving statements on gemini.
pub(crate) const PROTOCOL_NAME: &[u8] = b"GEMINI-v0";
/// The threshold for switching from space to time prover within the sumcheck.
const SPACE_TIME_THRESHOLD: usize = 22;
// const SUMCHECK_BUF_SIZE: usize = 1 << 20;

pub mod errors;
pub mod herring;
pub mod iterable;
pub mod kzg;
pub mod psnark;
pub mod snark;
pub mod subprotocols;

// private modules

#[doc(hidden)]
pub mod circuit;
pub mod misc;
mod transcript;
