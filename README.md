# The library overview

The `seqlibNaRo` library provides functionality for handling DNA sequences using the IUPAC nucleotide alphabet. It includes support for reading sequences from FASTA files, reverse complementing, GenBank-style slicing (1-based inclusive), and sequence manipulation through addition and object copying.

# Library installation using _pip_

Installation of the `seqlibNaRo`:

```Bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple seqlibNaRo
```

# A quick example of the library usage

Learn how to use the `seqlibNaRo` library with examples provided below.

```Python
from seqlibNaRo import DNASeq

# Example 1: Load from file
seqs = DNASeq.from_file("example.fasta")
print(seqs["seq1"])  # FASTA formatted output

# Example 2: Get reverse complement
rev = seqs["seq1"].revcmpl()
print(rev)

# Example 3: Slice a region (GenBank notation, inclusive)
fragment = seqs["seq1"][10:20]
print(fragment)

# Example 4: Add two sequences
combined = seqs["seq1"] + seqs["seq2"]
print(combined)
```

# The `DNASeq` class in details

The DNASeq class represents a single DNA sequence.
Methods:

- `from_file(file_path: str) -> List[DNASeq]` – loads multiple sequences from a FASTA file and returns a dictionary of DNASeq objects.

- `__str__()` – returns the sequence in FASTA format.

- `__repr__()` – returns a short preview of the sequence.

- `__len__()` – returns the length of the sequence.

- `__add__(self, other)` – concatenates two sequences.

- `revcmpl(self)` – returns the reverse complement as a new DNASeq object.

- `copy(self)` – returns a deep copy of the sequence object.

## Initialisation

`DNASeq(seqid: str, title: str, seq: str)`

## Indexing and slicing

You can slice using GenBank-style notation:

```Python
subseq = seq[5:10]  # includes both position 5 and 10 (1-based)
revsubseq = seq[10:5]  # returns reverse-complement of 5–10 region
```

Single index also follows GenBank notation:

```Python
nuc = seq[1]  # returns the first base
```

## Adding two sequences

Sequences can be concatenated with +:

```Python
new_seq = seq1 + seq2
```
