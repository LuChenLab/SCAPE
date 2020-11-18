# Single Cell Alternative Polyadenylation using Expectation-maximization (SCAPE)



## Usage

```
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  A pure python package for analysing alternative polyadenylation at single
  cell levels. Current version: 1.0.0

Options:
  --help  Show this message and exit.

Commands:
  apamix
  prepare
```

### prepare

prepare utr and intron region for inferring alternative polyadenylation

```
Usage: main.py prepare [OPTIONS]

Options:
  --gtf TEXT     The gtf (gz or not) file for preparing utr and intron region.
                 [required]

  --prefix TEXT  The prefix of output bed file.
  --help         Show this message and exit.
```

### apamix

```
Usage: main.py apamix [OPTIONS]

Options:
  --bed TEXT       The target regions (bed format) used for mix inference
                   [required]

  --bam TEXT       The bam file (sorted and indexed)  [required]
  -o, --out TEXT   The output path  [required]
  --cores INTEGER  Num (cores) of region are infering at once
  --cb TEXT        The cell barcode file, one cb for one line.  [required]
  -v, --verbose    Verbose mode
  --help           Show this message and exit.
```