# C and Assembly Implementations of do255e and do255s

This repository contains several implementations of do255e and
do255s, written in C and/or assembly.

These implementations were written for research and demonstration
purposes; however, they are all considered correct and safe, thus
usable in production. The C language being what it is, it is up
to the caller to perform integration of this code into the calling
application's source code and build system.

## Implementation List

For each of the do255e and do255s curves, the following implementations
are provided:

  - `w64`: C implementation meant for 64-bit platforms. It uses the
    intrinsic functions `_addcarry_u64()` and `_subborrow_u64()`,
    which are available on at least x86 platforms with the
    Clang, GCC and MSVC compilers. They _may_ be available on other
    architectures. This implementation can compile on at least
    Linux (Ubuntu 20.04), macOS (10.14.6 Mojave) and Windows 10
    (with MSVC 2019).

  - `bmi2`: a variant of `w64` with some inline assembly to leverage
    the BMI2 (`mulx`) and ADX (`adcx`, `adox`) opcodes. The inline
    assembly uses a syntax which is supported by GCC and Clang; it
    has been known to compile and run successfully on Linux
    (Ubuntu 20.04) and macOS (10.14.6 Mojave).

  - `w32`: a 32-bit variant of `w64`; it uses `_addcarry_u32()` and
    `_subborrow_u32()`. It is meant for 32-bit systems.

  - `cm0`: an implementation for ARM Cortex M0+ CPUs, written mostly
    in assembly.

## Code Roadmap

This repository contains many source files with an intrincate and
somewhat unusual structure: most C source files are not meant to be
compiled separately; instead, the source files _include_ each other. The
purpose of this construction is to make it so that, from the compiler's
point of view, all internal functions are `static` and visible. This is
supposed to help the compiler with inlining functions when necessary. It
also reduces namespace pollution at the link level.

Some source files are included from the main source files of several
distinct implementations. Notably, do255e and do255s use different base
fields, but still share the implementation; the relevant source file
relies on macros that define the exact field on which operations are to
be performed. We thus get inclusion patterns such as the following:

  - `w64` implementation for do255e: main file is `do255e_w64.c`; it
    includes:

      - `support.c`
      - `gf_do255e_w64.c`, which includes:
          - `gf_w64.c`
      - `sqrt_do255e_w64.c`
      - `padd_do255e_w64.c`, which includes:
          - `pcore_w64.c`
      - `icore_w64.c`
      - `scalar_do255e_w64.c`
      - `pmul_base_w64.c`
      - `pmul_do255e_w64.c`
      - `pvrfy_do255e_w64.c`
      - `pmap_do255e_w64.c`

  - `bmi2` implementation for do255e: main file is `do255e_bmi2.c`; it
    includes:

      - `support.c`
      - `gf_do255e_bmi2.c`, which includes:
          - `gf_bmi2.c`
      - `sqrt_do255e_w64.c`
      - `padd_do255e_w64.c`, which includes:
          - `pcore_w64.c`
      - `icore_w64.c`
      - `scalar_do255e_w64.c`
      - `pmul_base_w64.c`
      - `pmul_do255e_w64.c`
      - `pvrfy_do255e_w64.c`
      - `pmap_do255e_w64.c`

As seen above, the `bmi2` and `w64` implementations have many files in
common; only a couple of files differ (for the low-level implementation
of operations in the base field). Compilation of either implementation
is a single invocation of the C compiler on the main file, and it
produces a single object file.

For the ARM Cortex M0+ implementations, a similar mechanism is used, with
some differences:

  - Most of the included C files use the `arm` name.
  - Many functions are not implemented in C but in assembly.
  - The assembly code is provided in a collection of `.S` files. These
    files follow a similar pattern of inclusion, so that they are
    compiled with a single invocation of the C compiler, on the main
    assembly file, which is `do255e_cm0.S` or `do255s_cm0.S`.

Take care that the assembly file names use an uppercase `.S`, which
instructs the compiler to use the C preprocessor on them; this is how
inclusion works, and a few helper macros are used internally as well.

On top of that, the high-level cryptographic functionalities (key
exchange, signature...) are compiled separately, with `alg_do255e.c`
(for do255e) or `alg_do255s.c` (for do255s). Again, file inclusion is
used: both files include `alg.c`, in which the actual code is written.
The high-level functionalities are entirely implemented using the
low-level API (which is documented in `do255.h`) and thus are
architecture-agnostic.

A SHA-3 / SHAKE implementation is provided in `sha3.c` (with header
`sha3.h`). It is self-contained and compiled separately. It is used
by the high-level functionalities implementation, and by the tests
(`test_do255.c`). There is no dependency to SHA-3 or SHAKE in the
low-level API.

## API

The low-level API is documented in `do255.h`. This header defines all
types and functions that allow using the do255e and do255s curves,
as well as scalars for the prime order group (i.e. integers modulo
the group order). This API is deemed sufficient to implement arbitrary
cryptographic protocols over these groups without relying on any
detail of how the low-level operations are implemented.

High-level functionalities such as key exchange and signatures are
provided through the API documented in `do255_alg.h`. In this API,
points and scalars are opaque sequences of 32 bytes.

## Compilation

Type `make`. This should produce test binaries under the names
`test_do255_bmi2`, `test_do255_w64` and `test_do255_w32`, that run
internal tests and benchmarks. Benchmarks return values in clock cycles
(median over 1000 runs, as well as 10%-90% range over these 1000 runs).
For benchmarks, there is always some noise, hence measurement
variations, but they should be independent of any secret data. The CPU
cycle counter is used; if the test CPU has TurboBoost enabled, these are
likely not the "real" cycles. For reliable benchmarking, TurboBoost must
be disabled.

The default `Makefile` assumes that the compiler is Clang, and that the
current system is an Intel Skylake or newer. Adjust as needed.

With MSVC, use a Visual Studio command-line prompt, then type `nmake -f
Makefile.win32`. Only the `w32` and `w64` implementations will be built.
If targeting 32-bit mode, then the `w64` code will not compile; in that
case, compile `test_do255_w32.exe` explciitly.

For the ARM Cortex M0+ code, use `make -f Makefile.cm0`. This will
invoke a cross-compiler under the name `arm-linux-gcc`. To obtain an
appropriate cross-compiler, consider using
[Buildroot](https://buildroot.org/); this will allow the production of a
test binary `test_do255_cm0` that can then be executed with `qemu-arm`
(from [QEMU](https://www.qemu.org/)).

QEMU is very convenient for development and tests, with two caveats:

  - Under QEMU emulation, unaligned memory accesses work fine, but they
    trigger faults on the real hardware.

  - QEMU cannot provide speed benchmarks. These must be done on an
    actual Cortex M0+ CPU.

The ARM Cortex M0+ code has also been compiled with `arm-none-eabi-gcc`
(from the `gcc-arm-none-eabi` package on Ubuntu 20.04) with a custom
firmware skeleton, then executed on a SAM D20 Xplained Pro test board,
on which cycle-accurate benchmarks were performed.
