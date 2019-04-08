# mirvestigator-rust - Implementation of the miRvestigator algorithm in Rust

## Description

This is an implementation of miRvestigator in Rust. It is designed
as a standalone tool improve a few properties of the original version:

  - miRvestigator has generally been used within a larger context, either
    FIRM or SYGNAL or the miRvestigator website, but it makes sense to
    maintain it as a separate tool that can be called separately from
    various different contexts
  - Running data sets with a larger number of PSSMs leads to very long
    run times, sometimes a few days for a few thousand PSSMs, this is largely
    to 2 aspects:
        1. the computationally heavy parts are done in pure Python
        2. even though the nature of the data lends itself naturally
           to running computations concurrently there is a no parallelization

By choosing Rust a native implementation was realized that can be easily maintained
and enhanced if necessary. The port was relatively straightforward.

## Building

Development builds:

```
$ cargo build
```

Release builds:

```
$ cargo build --release
```

The difference in performance is pretty dramatic between development
and release version, so for production always compile for release

## Running

```
$ mirvestigator <SEQS> <PSSMS> <MIRBASE>
```

The parameters are as follows

  - seqs.txt file:

    a file that contains one DNA sequence per row, this will be used to optimize
    internal filtering

  - PSSM JSON file:

    A list of PSSM objects, the only required fields are

      - "name": the name of the PSSM
      - "matrix" the PSSM, which is a list of lists, n rows of 4 columns,
        assuming the alphabet A, C, G, T

  - mature.fa.gz file:
    a file downloaded from mirbase.org containing the miRNAs



The tool will write a list of entries to stdout in the format

```
<PSSM name>,<mirnas>,<match type>
```

mirnas is a list of micro RNAs that is separated with underscores.
The output format is compatible with the output format of the original
miRvestigator.

So one could run the tool and redirect its output into a file

```
$ mirvestigator <SEQS> <PSSMS> <MIRBASE> > results.csv
```

