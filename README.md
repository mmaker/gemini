<p align="center">
  <img src="doc/logo.svg" width=600px />
</p>


[Gemini](https://ia.cr/2022/420) is elastic proof system, FFT-free, blazingly fast and space-conscious.
**This code is **not** meant for production use and has not been audited.**


## Documentation

Gemini's API can be accessed via:

```bash
cargo rustdoc --open
```
This include a detailed protocol description.

## Features

Gemini can be compiled with the following feature flags:

- `asm`, to turn on the assembly backend within [`ark-ff`](https://docs.rs/ark-ff/);
- `parallel`, to turn on multi-threading. This requires the additional dependency [`rayon`](https://docs.rs/rayon/latest/rayon/);
- `std`, to rely on the Rust Standard library;
- `print-trace`, to print additional information concerning the execution time of the sub-protocols. This feature must be enabled in order to monitor the execution time of the examples.


## Performance

Micro-benchmarks for single components are done using [`criterion.rs`](https://github.com/bheisler/criterion.rs).
Benchmarks for the elastic and time-efficient provers can be run through the `examples/`.
More specifically, to benchmark an instance of *logarithmic* size `INSTANCE_LOGSIZE`, run:
```
cargo run --features asm,print-trace,parallel --example snark -- -i <INSTANCE_LOGSIZE>
```
You monitor the memory used setting the environment variable `RUST_LOG=debug`.
We tested so far ranges from 18 up to 35.
It is also possible to run a purely-linear time prover with the additional option `--time-prover`.

## Elasticity

Space footprint for the prover can be tweaked playing with the following constants:
- `TENSOR_EXPANSION_LOG` (set to 16) which sets the space budget for expanding the tensor products;
- `MAX_MSM_BUFFER_LOG` (set to 20) the size of the buffers over which Gemini performs multi-scalar multiplication;
- `SPACE_TIME_THRESHOLD` (set to 22) the threshold for converting the space-prover into the time-prover.
It will run a separate thread measuring stack+heap consumption.

The so-called `time_prover` modules contain a full-speed implementation of the proving algorithm without constraints on space.
