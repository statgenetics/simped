# SimPed

SimPed simulates haplotype and genotype data for pedigrees of arbitrary size
and structure. Founders are drawn from user-supplied haplotype frequencies
(haplotype mode) or allele frequencies under Hardy-Weinberg equilibrium
(genotype mode); non-founder haplotypes are passed down meiotically with
recombination governed by user-supplied inter-marker recombination fractions
or map distances. The two modes can be interleaved within a single
simulation. SimPed handles diallelic and multiallelic markers and has been
exercised on simulations with more than 20,000 loci.

This repository hosts a modified version of the original 2005 SimPed source
that compiles cleanly with current C compilers (clang or gcc, `-std=c99` or
later) on macOS, Linux, and Windows. The algorithm and output format are
unchanged; output is byte-identical to the original precompiled binaries on
all bundled examples.

**Reference.** Leal SM, Yan K, Müller-Myhsok B. *SimPed: a simulation program
to generate haplotype and genotype data for pedigree structures.* Human
Heredity 60(2):119–122 (2005). PMID: [16224187](https://pubmed.ncbi.nlm.nih.gov/16224187/).

---

## Contents

1. [Installation](#installation)
2. [Running SimPed](#running-simped)
3. [Input files](#input-files)
   - [Pedigree file](#pedigree-file)
   - [Parameter file](#parameter-file)
4. [Worked examples](#worked-examples)
5. [Output file](#output-file)
6. [Algorithm](#algorithm)
7. [Compile-time limits](#compile-time-limits)
8. [Changes from the upstream source](#changes-from-the-upstream-source)
9. [License and contact](#license-and-contact)

---

## Installation

### Build from source (recommended)

A C99 (or later) compiler and `make` are the only requirements. On macOS the
Xcode Command Line Tools provide both (`xcode-select --install`); on Linux,
any distribution package containing `gcc` or `clang` plus `make` will do.

```sh
make                              # builds ./simped
make debug                        # debug build (-g -O0)
make check                        # builds and runs two example inputs
make install                      # installs to /usr/local/bin
make install PREFIX=$HOME/.local  # user-local install, no sudo
```

### Use a precompiled binary (not maintained; use as-is)

The `executable/` directory contains the original 2005/2007 distribution
binaries, organized by platform.

| Platform | Binary                        | Notes                                                                                                                                  |
| -------- | ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Linux    | `executable/linux/simped`     | ELF 32-bit Intel. On 64-bit Linux the runtime libraries `libc6-i386` (Debian/Ubuntu) or `glibc.i686` (Red Hat/Fedora) must be installed |
| macOS    | `executable/macos/simped`     | Mach-O x86_64. Runs on Apple Silicon under Rosetta 2. Native build via `make` is preferred                                              |
| Windows  | `executable/windows/simped.exe` + `cygwin1.dll` | Both files must reside in the same directory                                                                |

Copy the binary onto a directory in your `PATH` to invoke it as `simped` from
any location.

---

## Running SimPed

```sh
simped <parameter_file>
```

The parameter file declares the input pedigree filename and the output
filename on its first line. Both paths are resolved relative to the current
working directory.

```
$ cd examples && ../simped input1.dat
  ********************************************************
  *                                                      *
  *               Program   SimPed                       *
  *                                                      *
  ********************************************************

Usage: simped parameter_file

Constants in effect:
  Maximum number of pedigrees                   100
  Maximum number of individuals                 3000
  Maximum number of loci                        10000
  Maximum number of haplotypes                  1000
  Maximum number of alleles at a locus          20

Input pedigree file: pedin1.pre
Output file: pedfile1.pre
```

---

## Input files

SimPed reads two whitespace-separated text files: a *pedigree file* in
standard LINKAGE format and a *parameter file* controlling the simulation.

### Pedigree file

LINKAGE-style. Each line is one individual.

| Column | Meaning                                                                |
| ------ | ---------------------------------------------------------------------- |
| 1      | Pedigree identifier                                                    |
| 2      | Individual ID                                                          |
| 3      | Father's ID (`0` if founder)                                           |
| 4      | Mother's ID (`0` if founder)                                           |
| 5      | Sex (`1` = male, `2` = female)                                         |
| 6–N    | Affection status and/or quantitative trait columns (may be omitted)    |
| N+1    | Marker mask: which marker genotypes to emit for this individual        |

The marker mask resolves as follows:

- Single `0`: no genotypes emitted, all markers blanked to `0 0`.
- Single `1`: all markers emitted.
- A string of `0`/`1` of length equal to the total number of markers: per-marker mask, where `1` emits the simulated genotype and `0` blanks it.

Pedigree and individual IDs may be numeric, alphabetic, or mixed. Multiple
pedigrees may share one file, delimited by a change in the pedigree
identifier column.

### Parameter file

The header is seven fixed lines; from line 8 onward the file is a sequence of
haplotype and/or genotype blocks.

| Line | Contents                                                                                                                                                          |
| ---- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1    | `<pedigree_file> <output_file>`                                                                                                                                   |
| 2    | Three random number generator seeds, each less than 30000                                                                                                         |
| 3    | Number of trait/affection columns between sex (column 5) and the marker mask                                                                                      |
| 4    | Number of replicates                                                                                                                                              |
| 5    | Total number of marker loci, then the number of times the block pattern (lines 8+) is repeated to fill those loci                                                 |
| 6    | Inter-marker distance type: `1` recombination fractions, `2` Kosambi map distances (Morgans), `3` Haldane map distances (Morgans)                                 |
| 7    | `<count> <values…>`: number of inter-marker distances supplied followed by the values themselves. See "Inter-marker distance expansion" below                     |

**Line 5 pattern repetition.** `totalLoci / iteration` must be a whole
number. If a block pattern covers 13 markers and 1300 loci are wanted, set
line 5 to `1300 100`; the blocks from line 8 onward are replayed 100 times.

**Inter-marker distance expansion.** The `<count>` field on line 7 admits
three forms:

- `count = 1`: one value, used for every inter-marker interval.
- `count = totalLoci − 1`: every interval supplied explicitly.
- `1 < count < totalLoci − 1`: the supplied values cycle as a repeating pattern across the `totalLoci − 1` intervals.

Kosambi and Haldane map distances are converted internally to recombination
fractions using $\theta = \tfrac{1}{2}(e^{4d}-1)/(e^{4d}+1)$ (Kosambi) and
$\theta = \tfrac{1}{2}(1-e^{-2d})$ (Haldane), with $d$ in Morgans.

#### Haplotype block (line begins with `1`)

```
1 <n_haplotypes> <n_loci> <n_repeats>
<freq_1> <freq_2> ... <freq_n_haplotypes>
<allele_1> <allele_2> ... <allele_n_loci>     ← haplotype 1 alleles
<allele_1> <allele_2> ... <allele_n_loci>     ← haplotype 2 alleles
...
```

Loci within a block are in linkage disequilibrium and share one set of
haplotype frequencies; distinct blocks are in linkage equilibrium with one
another.

#### Genotype block (line begins with `2`)

```
2 <n_loci> <n_repeats>
<n_alleles_locus_1> <freq_a1> <freq_a2> ... <freq_a(n-1)>
<n_alleles_locus_2> <freq_a1> <freq_a2> ... <freq_a(n-1)>
...
```

Only `n_alleles − 1` frequencies are listed per locus; the final frequency is
inferred. Genotypes are simulated under Hardy-Weinberg equilibrium.

Haplotype and genotype blocks may alternate freely. The total number of
marker loci across all blocks in one round of the pattern must equal
`totalLoci / iteration` from line 5.

---

## Worked examples

Each example below has a parameter file `input<k>.dat` and a pedigree file
`pedin<k>.pre` in `examples/`. Examples 1 to 4 reuse one nuclear pedigree of
six individuals (parents 1 and 2, children 3 to 6) with body mass index as
the quantitative trait. Example 5 demonstrates a multi-family pedigree file.

### Example 1: pure haplotypes

10 marker loci forming 8 haplotypes, all in linkage disequilibrium, 1000
replicates, recombination fractions supplied directly.

```
pedin1.pre pedfile1.pre  << name of pedigree file, name of output file
23221 1601 21001         << three random seeds
0                        << number of trait columns
1000                     << number of replicates
10 1                     << total markers, pattern repeats
1                        << recombination fractions on line 7
9 0.001 0.0 0.005 0.002 0.0 0.01 0.002 0.0 0.001
1 8 10 1                 << haplotype block: 8 haplotypes, 10 loci, repeat 1×
0.35 0.15 0.10 0.08 0.08 0.08 0.08 0.08          << haplotype frequencies
1 1 1 1 1 1 1 1 1 1                              << haplotype 1
2 2 2 2 1 1 1 1 1 1
1 1 2 2 2 2 2 1 1 2
1 2 1 2 1 2 1 2 1 2
1 1 1 1 2 2 2 2 2 1
2 2 2 1 1 1 1 2 2 2
1 1 2 2 1 1 2 2 1 1
2 2 2 1 1 2 2 2 2 1
```

### Example 2: haplotypes with Kosambi map distances, per-individual mask, multiallelic markers

10 marker loci forming 10 haplotypes, 10000 replicates, BMI included as a
quantitative trait. Markers 2 and 9 have more than two alleles. The pedigree
file uses per-individual marker masks to blank specific markers for
individuals 4 and 6:

```
1 1 0 0 1 22.3  0                          ← individual 1 deceased, all blanked
1 2 0 0 2 30.3  1                          ← all markers emitted
1 3 1 2 1 29.6  1
1 4 1 2 2 24.5  1 1 1 1 0 1 1 0 1 1        ← markers 5 and 8 blanked
1 5 1 2 1 35.7  1
1 6 1 2 2 23.5  1 1 0 1 1 1 1 1 0 0        ← markers 3, 9, 10 blanked
```

Parameter file `input2.dat` opens with:

```
pedin2.pre pedfile2.pre  << name of pedigree file, name of output file
21998 26002 11981        << three random seeds
1                        << one trait column (BMI)
10000                    << number of replicates
10 1                     << total markers, pattern repeats
2                        << Kosambi distances on line 7
9 0.001 0.02 0.005 0.02 0.001 0.001 0.002 0.01 0.001
1 10 10 1                << haplotype block: 10 haplotypes, 10 loci
0.30 0.25 0.10 0.10 0.05 0.04 0.04 0.04 0.04 0.04
1 3 1 1 1 1 1 1 4 1                        ← haplotype 1 (alleles 3, 4 multiallelic)
...
```

### Example 3: pure genotypes under Hardy-Weinberg

20 marker loci in linkage equilibrium, mixing diallelic, triallelic, and
hexallelic markers, 1000 replicates. Allele-frequency blocks are repeated to
fill the 20 loci.

```
pedin3.pre pedfile3.pre  << name of pedigree file, name of output file
1344 2673 12228          << three random seeds
0                        << number of trait columns
1000                     << number of replicates
20 2                     << total markers, pattern repeats 2× to fill 20 loci
1                        << recombination fractions on line 7
10 0.01 0.05 0.05 0.09 0.01 0.01 0.02 0.05 0.01 0.05
2 3 2                    << genotype block: 3 loci, pattern repeated 2×
2 0.3                    << locus 1: 2 alleles, first allele freq 0.3
2 0.4                    << locus 2
2 0.2                    << locus 3
2 2 1                    << genotype block: 2 loci
3 0.7 0.1                << locus, 3 alleles
6 0.1 0.1 0.1 0.1 0.1    << locus, 6 alleles
2 1 2                    << genotype block: 1 locus, pattern repeated 2×
2 0.5
```

### Example 4: mixture of haplotypes and genotypes

220 marker loci total, generated from a 22-marker block pattern repeated 10
times. The block contains one genotype block (6 markers) followed by two
haplotype blocks (5 markers and 2 markers). This case exercises every code
path: pattern repetition, alternating block types, masking via the pedigree
file, and trait columns.

### Example 5: multi-family pedigree file

`examples/pedin5.pre` interleaves two pedigrees (IDs `3` and `1`) in one
file. SimPed processes them independently; individual IDs need not be
sequential within a family.

---

## Output file

Output follows LINKAGE pedigree format. One line per individual per
replicate:

```
<replicate_pedigree_id>  <ind_id>  <father>  <mother>  <sex>  [trait_columns…]  <a1_l1> <a2_l1>  <a1_l2> <a2_l2>  …
```

Replicates are concatenated; the pedigree-id column is re-numbered across
replicates so the file remains a valid multi-pedigree LINKAGE file. The
first two lines of `pedfile1.pre` from Example 1 (truncated to 20 markers)
read:

```
1  1  0  0  1   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1  2  0  0  2   2 2 2 2 2 2 1 1 1 1 2 1 2 1 2 2 2 2 1 2
```

For a post-MAKEPED file with first-offspring, paternal/maternal sibling, and
proband columns, run the output through `makeped` from the LINKAGE package.

---

## Algorithm

For each replicate the program executes the following.

1. **Founder assignment.** For each founder (both parental IDs `0`) and each
   block of marker loci specified from line 8 onward, two haplotypes are
   drawn independently. In a haplotype block, a uniform random number is
   compared to the cumulative haplotype-frequency vector to select one of
   the supplied haplotype rows. In a genotype block, alleles at each locus
   are drawn independently from the supplied allele-frequency vector.
2. **Meiotic drop-down.** Within each pedigree, every individual whose
   parents are both filled is visited in turn. For each parental
   contribution (paternal and maternal), one of the parent's two haplotypes
   is selected uniformly at random at locus 1. At each subsequent locus *j*,
   a uniform draw is compared to the inter-marker recombination fraction
   $\theta_{j-1}$; with probability $\theta_{j-1}$ the haplotype of origin
   switches, otherwise it persists. The process repeats over all `totalLoci`
   loci. A while-loop over the family ensures that individuals are filled in
   topological order regardless of input ordering.
3. **Masking and output.** The pedigree file is read a second time; the
   marker mask column is consulted per individual to decide whether each
   locus is emitted as the simulated genotype or as `0 0`. The pedigree
   identifier is re-numbered to remain unique across replicates.

The random number generator is a three-component Wichmann-Hill style
generator seeded by line 2 of the parameter file. Identical seeds reproduce
identical output.

---

## Compile-time limits

Hard limits are set as `#define` constants in `src/simped.c`. Adjust and
rebuild for larger inputs.

| Constant         | Default | Meaning                              |
| ---------------- | ------- | ------------------------------------ |
| `MAX_PED`        | 100     | maximum number of pedigrees          |
| `MAX_IND`        | 3000    | maximum number of individuals        |
| `MAX_LOCUS`      | 10000   | maximum number of loci               |
| `MAX_HAPLOTYPES` | 1000    | maximum number of haplotypes         |
| `MAX_ALLELE`     | 20      | maximum number of alleles at a locus |

---

## Changes from the upstream source

The 2005 source did not compile on macOS or under recent gcc/clang because
of a name collision and several implicit-int constructs that ceased to be
valid C in C99. Behavior is unchanged; the diff is mechanical.

- The global `index` was renamed to `founder_marker_idx`. The original name
  collides with the POSIX `index()` function declared in `<strings.h>`,
  which is transitively included on macOS and BSD platforms.
- Explicit return types were added to `main()` and `writeResult()`.
- Forward declarations were added for all internal functions.
- Field-width caps were added to `scanf` formats and return values of
  `fopen` and `fscanf` are now checked.
- The haplotype-frequency cumulative-sum loop no longer reads
  `possibilities[-1]` on its first iteration.
- A small number of stack-allocated variable-length arrays were replaced
  with fixed-size buffers tied to the existing `MAX_LOCUS` and `MAX_ALLELE`
  constants. This removes a portability hazard without changing capacity.

Output from the rebuilt binary is byte-identical to the original
precompiled binaries on every example in `examples/`.

---

## License and contact

The upstream distribution did not include an explicit software license.
The original program and documentation is © Suzanne M.
Leal (2005). For redistribution terms or questions about the underlying
algorithm, contact the original authors.
