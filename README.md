# Bio_files_processor

This `.py` script contains two functions: `convert_multiline_fasta_to_oneline()`
and `parse_blast_output()`

## Main script functions

### `convert_multiline_fasta_to_oneline()`

Accepts 2 arguments as input: `input_fasta` and `output_fasta`.
Reads the fasta file submitted as input_fasta, in which the
sequence (DNA/RNA/protein/ ... ) can be split into several
lines, after which it saves it to a new fasta file in which
each sequence fits into one line

    Arguments: input_fasta: 'input_file_name',
    output_fasta: 'output_file_name'

    Returns None, creates output .fasta file

### `parse_blast_output()`

Accepts 2 arguments as input: input_file and output_file.
The function reads the specified file, selects the first
row from the Description column for each QUERY (Sequences
producing significant alignments). Saves the list of obtained
proteins to a new file in one column sorted alphabetically.

    Arguments: input_file: 'input_file_name',
    output_file: 'output_file_name'

    Returns None, creates output .txt file

# Bioinformatic_tools

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

Function accepts **FASTQ file** as first argument and discards
sequences which do not meet criteria for *length, GC content* and *quality*
Remaining reads are written in output file in **`/filtered` folder**
as FASTQ dictionary
Requires the import of **filter_fastq_aux_mod.py**

**Arguments**:

Function filter_fastq accepts 5 arguments:
`input_fastq`, `output_fastq`, `gc_bounds`, `length_bounds`, `quality_threshold`:

* `input_fastq` - input file in FASTQ format

* `output_fastq` - file containing FASTQ dictionary
in format `{name: ("sequence, quality")}`

* `gc_bounds` - the interval of the GC composition (in percent)
(by default is (0, 100), i.e. all reads are saved). If a single number
is passed to the argument, it is assumed that this is the upper bound.
Examples: gc_bounds = (20, 80) - save only reads with GC composition
from 20 to 80%, gc_bounds = 44.4 - save reads with GC composition less
than 44.4%

* `length_bounds` - the length interval for filtering, everything
is the same as for `gc_bounds`, but by default it is (0, 2**32)

* `quality_threshold` - the threshold value of the average read quality for
filtering. Is 0 by default (phred33 scale). Reads with average
quality for all nucleotides below the threshold are discarded

**Returns**: None

## Auxiliary functions

+ `is_dna()` check whether the argument is DNA sequence. Works with upper
and lower case letters
+ `is_rna()` check whether the argument is RNA sequence. Works with upper
and lower case letters
+ `transcribe()` returns transcribed argument sequence
+ `reverse()` returns reverse argument sequence
+ `complement()` returns complement argument sequence
+ `reverse_complement()` returns reverse complement argument sequence
+ `gc_content()` returns percent of G + C in argument sequence
+ `seq_quality()` accepts `'quality'` `str` line from form FASTQ dictionary in Phred33 code
and returns mean `float` quality index for the whole sequence
+ `fastq_to_dict()` accepts FASTQ file as input. Returns dictionary containing reads in format
`{name: (sequence, quality)}`
+ `write_in_file()` accepts FASTQ reads in format `{name: (sequence, quality)}` as input.
Returns file containing this dictionary.
