# Bioinformatic_tools

## What is it?

The package is a toolkit for working with DNA sequences.
It contains two main functions `run_dna_rna_tools()` and
`filter_fastq()` and their dependencies.

## Main script functions

### `run_dna_rna-tools()`

Function accepts **at least two** `str` arguments.
All arguments except for the last one are DNA or RNA sequences.
The last argument is the action required for the input 
sequences: `'transcribe'`, `'reverse'`, `'complement'` or `
'reverse complement'`
Requires the import of **dna_rna_aux_mod.py**

**Arguments**: `str`

**Returns**: `str | list`

#### Usage examples
```
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```

### `filter_fastq()`

Function accepts **FASTQ dictionary** as first argument and discards
sequences which do not meet criteria for *length, GC content* and *quality*.
Requires the import of **filter_fastq_aux_mod.py**

**Arguments**:

Function filter_fastq accepts 4 arguments:
`seqs`, `gc_bounds`, `length_bounds`, `quality_threshold`:

* `seqs` - a dictionary consisting of fastq sequents. The structure is as follows.
The key is a string, the name of the sequence. The value is a tuple of two
strings: sequence and quality.

* `gc_bounds` - the interval of the GC composition (in percent)
(by default is (0, 100), i.e. all reads are saved). If a single number
is passed to the argument, it is assumed that this is the upper bound.
Examples: gc_bounds = (20, 80) - save only reads with GC composition
from 20 to 80%, gc_bounds = 44.4 - save reads with GC composition less
than 44.4%.

* `length_bounds` - the length interval for filtering, everything
is the same as for `gc_bounds`, but by default it is (0, 2**32).

* `quality_threshold` - the threshold value of the average read quality for
filtering. Is 0 by default (phred33 scale). Reads with average
quality for all nucleotides below the threshold are discarded.

**Returns**: `dict`

## Auxiliary functions

### `is_dna()`

Check whether the argument is DNA sequence. Works with upper
and lower case letters

### `is_rna()`

Check whether the argument is RNA sequence. Works with upper
and lower case letters

### `transcribe()`

Returns transcribed argument sequence

### `reverse()`

Returns reverse argument sequence

### `complement()`

Returns complement argument sequence

### `reverse_complement()`

Returns reverse complement argument sequence

### `gc_content()`

Returns percent of G + C in argument sequence

### `seq_quality()`

Accepts `'quality'` `str` line from form FASTQ dictionary in Phred33 code
and returns mean `float` quality index for the whole sequence
