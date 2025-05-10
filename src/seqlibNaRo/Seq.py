from typing import List, Self
import numpy as np

class DNASeq:
    """
    A class to represent a DNA sequence using IUPAC nucleotide code.
    Supports reverse-complementation, slicing using GenBank notation,
    FASTA formatting, and object-based operations like addition.
    """

    ALPH = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N', '-': '-'
    }

    def __init__(self, seqid: str, title: str, seq: str) -> None:
        """Initializes a new DNASeq object with given ID, title, and sequence."""
        self.seqid = seqid
        self.title = title
        self.seq = seq

    @classmethod
    def from_file(cls, file_path: str) -> List[Self]:
        """
        Parses a FASTA file and returns a dictionary of DNASeq objects.

        Raises:
            Exception: If sequence contains invalid characters or seqid is not unique.
        """
        seqs, seqid = {}, ''
        with open(file_path) as f:
            lines = f.readlines()
        for line in lines:
            if line[0] == '>':
                header = line[1:-1].split(' ', 1)
                seqid = header[0]
                title = header[1] if len(header) > 1 else ''
                if seqid in seqs:
                    raise Exception(f'Non-unique: {seqid}')
                seqs[seqid] = cls(seqid, title, [])
            else:
                seq_to_append = line[:-1]
                if set(seq_to_append).issubset(cls.ALPH.keys()):
                    seqs[seqid].seq.append(seq_to_append)
                else:
                    raise Exception(f'Nucleotide sequence "{seq_to_append}" not in alphabet!')
        for seq in seqs.values():
            seq.seq = ''.join(seq.seq)
        return seqs

    def __repr__(self) -> str:
        """Returns a short representation of the sequence (first 10 nucleotides)."""
        return f'{self.seq[:10] if len(self.seq) > 10 else self.seq}...'

    def __len__(self) -> int:
        """Returns the length of the sequence."""
        return len(self.seq)

    def __str__(self) -> str:
        """Returns the sequence in FASTA format (60-character lines)."""
        step = 60
        lines = '\n'.join([self.seq[i:i+step] for i in range(0, len(self.seq), step)])
        title = f' {self.title}' if self.title != '' else ''
        return f'>{self.seqid}{title}\n{lines}'

    def revcmpl(self) -> Self:
        """Returns a new DNASeq object that is the reverse complement of the current one."""
        seq = self.seq[::-1]
        seq_revcmpl = ''.join([self.ALPH[nuc] for nuc in seq])
        seqid = f'{self.seqid}_revcmpl'
        return DNASeq(seqid, self.title, seq_revcmpl)

    def __getitem__(self, key):
        """
        Supports GenBank-style slicing (1-based, inclusive, reverse-aware).
        Returns a single character (str) or a sliced DNASeq object.
        """
        if isinstance(key, slice):
            start, end, step = key.start, key.stop, key.step
            if step is not None:
                raise KeyError('Step is not allowed in GenBank notation')
            if start is None or end is None:
                raise KeyError('Start and end indices are required')
            if not np.issubdtype(type(start), np.integer) or not np.issubdtype(type(end), np.integer):
                raise TypeError('Start and end must be integers')
            if start <= 0 or end <= 0:
                raise KeyError('Minimal value for start and end is 1')
            strand = 1 if start <= end else -1
            start, end = sorted([start, end])
            start -= 1
            seqid = f'{self.seqid}_loc({start+1}_{end})'
            sub_seq = DNASeq(seqid, self.title, self.seq[start:end])
            if strand == -1:
                sub_seq = sub_seq.revcmpl()
            return sub_seq
        else:
            if not np.issubdtype(type(key), np.integer):
                raise TypeError('Index must be an integer')
            if key <= 0:
                raise KeyError('Minimal value of index is 1')
            return self.seq[key - 1]

    def __add__(self, other_seq: Self) -> Self:
        """
        Concatenates two DNASeq objects and returns a new DNASeq object.
        """
        return DNASeq(f'{self.seqid}_{other_seq.seqid}',
                      f'{self.title}_{other_seq.title}',
                      f'{self.seq}{other_seq.seq}')

    def copy(self) -> Self:
        """Returns an exact copy of the current DNASeq object."""
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result
