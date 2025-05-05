from typing import List, Self
import numpy as np

# And there it is, our DNASeq class!
class DNASeq:
 
    #===SECTION=1========================================================
    
    # A dictionary with the IUPAC ambiguous DNA alphabet
    # (https://www.bioinformatics.org/sms/iupac.html)
    # and a gap sign. Each key is paired with a value
    # that is the key's complement. This dictionary
    # is defined as a class variable. Revise
    # how to access values of such variables,
    # you will need it in the coming sections.
    # TODO: Fill up the missing letters,
    #       make them capitals.
    
    ALPH = {
       'A' : 'T',   'T' : 'A',   'G' : 'C',   'C' : 'G',
       'K' : 'M',   'M' : 'K',   'B' : 'V',   'D' : 'H',
       'H' : 'D',   'V' : 'B',   'N' : 'N',   '-' : '-'
    }
    
    # TEST the code in the Section 1.
    
    
    #===SECTION=2========================================================
    
    # The special method __init__() initialising the
    # initial state of a newly created object (a class instance).
    # TODO: Define the special method __init__() proper for
    #       initialisation of a new object of DNASeq type,
    #       which is to be created as following:
    #       seq = DNASeq(seqid, title, seq)
    
    def __init__(self, seqid: str, title: str, seq: str) -> None:
        # TODO: Assign values of the arguments passed
        #       while creating a new object to the
        #       object attributes that will have
        #       the same names.
        
        self.seqid = seqid
        self.title = title
        self.seq = seq
        
    # TEST the code in the Section 2.
        
        
    #===SECTION=3========================================================
        
    # Now we need a class method that will help us to
    # deserialise sequences from a FASTA format file to
    # a collection of DNASeq objects.
    # TODO: Define a class method from_file() that next to
    #       the automatically passed reference to a class
    #       will accept an argument with a file path, and
    #       will work as following:
    #       seqs = DNASeq.from_file('some/path/some_file.fasta')
    
    @classmethod
    def from_file(cls, file_path: str) -> List[str]:
        # TODO: Create an empty seqs dictionary and a variable for
        #       a current seqid.
        seqs, seqid = {}, ''
        
        # TODO: Open the input file for reading
        with open(file_path) as f:
            lines = f.readlines()

            
        # TODO: Start iteratively parsing the opened file
        #       line by line. Refer to the lecture presentation
        #       for details.
        for line in lines:
            # TODO: Parse sequences form the input file.
            # TODO: An extra step: convert every sequence line
            #       to a set and use a proper method for Python sets,
            #       which will allow you to see if the
            #       alphabet of seq fragment is a subset of ALPHABET,
            #       which is one of previously defined
            #       class variables. If not raise an exception.
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
                    raise Exception('Nucleotide sequence "{seq_set}" not in alphabet!')
        for seq in seqs.values():
            seq.seq = ''.join(seq.seq)
            
        # TODO: Return the reference to the collection of new objects
        #       of the cls class.
        return seqs
    
    # TEST the code in the Section 3.
    
        
    #===SECTION=4========================================================
    
    # The special method __repr__() returns a string representation of
    # an object, we will define this representation as 10 first letters
    # of the sequence stored within an object and three dots '...'
    # UNLESS the sequence length is shorter or equal to 10 letters.
    
    def __repr__(self) -> str:
        return f'{self.seq[:10] if len(self.seq) > 10 else self.seq}...'
    
    # TEST the code in the Section 4.
    
    
    #===SECTION=5========================================================
        
    # TODO: Define the proper special method that will return the length of
    #       the sequence contained within a given DNASeq object when the reference
    #       to that object is passed to the built-in len() function.
    
    def __len__(self) -> int:
        return len(self.seq)
    
    # TEST the code in the Section 5.
    
    
    #===SECTION=6========================================================
    
    # TODO: Implement the special method __str__(), which returns
    #       a string whenever a DNASeq object is being converted to one,
    #       eg. str(obj) or print(obj).
    #       Make it return the sequence contained within an object 
    #       as FASTA formatted sequence string.
    #       Let the sequence be divided into 60-character lines.
    
    def __str__(self):
        step = 60
        lines = '\n'.join(
           [self.seq[i:i+step] for i in range(0, len(self.seq), step)]
        )

        title = f' {self.title}' if self.title != '' else ''

        fasta = f'>{self.seqid} {title}\n{lines}'

        return fasta
    
    # TEST the code in the Section 6.
    
    
    #===SECTION=7=======================================================
    
    # Let's create a custom method and call it revcmpl().
    # It will simply return a new object of DNASeq type, which
    # will contain a reverse-complement sequence to the
    # one contained within the original object the method is called from.
    # To make both object distinguishable, create a new
    # seqid for the sequence by adding '_revcmpl' suffix
    # to the seqid of the original one.
    
    def revcmpl(self):
        # TODO: Using string slices and step equal to -1,
        # reverse the sequence.
        
        seq = self.seq[::-1]
        
        # TODO: Using string method join(), the class dictionary ALPH and a
        #       list comprehension expression, translate the reversed sequence and
        #       convert into a string
        
        seq_revcmpl = ''.join([self.ALPH[nuc] for nuc in seq])
        
        # TODO: Create seqid variable and assign to it the object's seqid
        #       together with the suffix '_revcmpl'.
        
        seqid = f'{self.seqid}_revcmpl'
        
        # TODO: Create a new object of the DNASeq type using the new seqid,
        #       title contained in the object as well as
        #       reversed and translated sequence, return the reference
        #       to new object.
        
        return DNASeq(seqid, self.title, seq_revcmpl)
    
    # TEST the code in the Section 7.
    
    
    #===SECTION=8========================================================
    
    # Beware! This is the hardest part to go through.
    # We will implement the special method __getitem__()
    # that allows to program what happens when an object
    # is indexed or sliced (like a list or a string),
    # eg. obj[3] or obj[4:5].
    #
    # We know that indexing and slicing in Python works as follows:
    # - indexing starts at 0
    # - when a slice is taken [3:4] the first index is included
    #   the last not: eg. s = 'abcde', s[1:3] -> 'bc'
    #                          01234
    # - slice [3:1] will be an empty sequence
    #
    # When it comes to biological sequences, those are indexed
    # according to GenBank notation, which is a bit different:
    # - indexing starts at 1
    # - when a fragment is requested, both indices are inclusive:
    #   seq = 'ATGCTACG', seq[1:3] -> 'ATG'
    #          12345678
    # - the start index greater than stop index indicates
    #   a reverse complement (the complementary strand):
    #   seq[3:1] -> 'CAT'
    #
    # Now we will try to implement this behaviour in case of
    # our objects of DNASeq type.
    
    def __getitem__(self, key):
        
        if isinstance(key, slice):
            # If the key is a slice object, it has three properties:
            # start, stop and step, eg. list[start:stop:step].
            # If any is missing its value is None, eg. list[start:end].
            # First two will be equivalent of beginning and end
            # in GenBank notation.
            
            # The whole point here is to translate GenBank indices
            # into Python ones and index or slice the sequence
            # contained in the DNASeq object, which is surely of
            # Python string format, so it cannot be directly indexed
            # or sliced with GenBank indices without translating them.
            
            # To make the code nicer, let's assign values of slice
            # properties to single variables.
            
            start, end, step = key.start, key.stop, key.step
            
            # Let's test a few possibilities and exclude those
            # that cannot be translated into GenBank notation.
            
            if step is not None:
                # There is not an equivalent of step in GenBank notation.
                # If step is defined anyway, raise a KeyError.
                
                raise KeyError('Step is not allowed in GenBank notation')
            
            if start is None:
                # Start must be provided.
                
                raise KeyError('Start index is required')
            
            if end is None:
                # As well as the end.
                
                raise KeyError('End index is required')
                
            if not np.issubdtype(type(start), np.integer) or \
              not np.issubdtype(type(end), np.integer):
                # Start and end must be defined as integer values,
                # otherwise raise a KeyError.
                
                raise TypeError('Start and end must be integers')
            
            if start <= 0 or end <= 0:
                # Both must be greater than 0 as in GenBank notation
                # indexing starts at 1.
                
                raise KeyError('Minimal value for start and end is 1')
                
            # Now we are sure there is only start and end in the slice,
            # and that both are integer type values and equal or greater to 1.
            # Let's move on then.
            
            # TODO: Create a variable strand and set it to 1 if start is
            #       less or equal to end, otherwise to -1.
                
            strand = 1 if start <= end else -1
            
            # TODO: If strand is equal to -1, swap the values of start and end,
            #       by using tuple unpacking, so start is less than end (they are ordered).
            #       You may also do it without actually checking whether one is grater
            #       than the other, using build-in sorted() of numpy np.sort() function and
            #       tuple unpacking expression.
            
            start, end = sorted([start, end])
            # if strand == -1:
            #     temp = start
            #     start = end
            #     end = temp
            
            # TODO: Create a new seqid by adding to the existing one
            #       the suffix '_loc(start_end)', where start and end are
            #       values of our variables.
            
            seqid = f'{self.seqid}_loc({start}_{end})'
            
            # TODO: Decrease start by 1. If start in GenBank notation is 1,
            #       it is and equivalent of 0 in Python indexing (etc.),
            #       so it needs to be decreased by 1 to be translated
            #       from GenBank to Python.
            #
            #       You must NOT decrease end, as in GenBank it is inclusive,
            #       but in Python exclusive, so it is decreased
            #       somewhat automatically.
            
            start -= 1
            
            # TODO: Create a new object of DNASeq type by using new seqid
            #       you created a few lines before, title of the existing object,
            #       a slice of the sequence contained in that object by using
            #       the adjusted start index and the end index (string slicing).
            
            sub_seq = DNASeq(seqid, self.title, self.seq[start:end])
            
            # TODO: If strand is -1, assign to the same reference ("variable")
            #       a reverse complement of the new sequence object by using
            #       the method revcmpl() you implemented before.
            
            if strand == -1:
                sub_seq = sub_seq.revcmpl()
                
            # TODO: return the reference to the new sequence object.
                
            return sub_seq
                
        else:
            
            # If we are here, it means that key is not a slice object.
            # Then we allow it be an integer greater than 0,
            # compliant with the GenBank notation.
            
            if not np.issubdtype(type(key), np.integer):
                raise TypeError('Index must be an integer')
                
            if key <= 0:
                raise KeyError('Minimal value of index is 1')
                
            # In case it is just one letter not a slice, we
            # will return simply a letter from the sequence
            # (a string value), NOT a new DNASeq object.
            # We need to remember to decrease the key value
            # by one to translate it from GenBank to Python.
                
            return self.seq[key - 1]
    
    # TEST the code in the Section 8.
        
        
    #===SECTION=9========================================================
    
    # Implement the special method __add__(), which will be
    # invoked when two objects of the DNASeq type are being added,
    # eg. seq_3 = seq_1 + seq_2
    #
    # Make the result of addition a new object of type(self) type
    # with its properties values as follows:
    # - seqid will be seqids of added objects,
    #   separated by an underscore
    # - title will be titles of added objects,
    #   separated by an underscore
    # - seq will simply be a concatenation of both sequences
    #   stored within the added objects
    #
    # return the reference to the new object.
    
    def __add__(self, other_seq: Self) -> Self:
        return DNASeq(f'{self.seqid}_{other_seq.seqid}',
                      f'{self.title}_{other_seq.title}',
                      f'{self.seq}{other_seq.seq}')

    
    # TEST the code in the Section 9.
    
    
    #===SECTION=10========================================================
        
    # Define a custom method copy() that will return an exact copy
    # of the existing object (as a type(self) type object).

    def copy(self) -> Self:
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result
    
    # TEST the code in the Section 10.
    