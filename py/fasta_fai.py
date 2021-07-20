# Read fasta file

import os
import math

class Index:
    """
    Fasta index
    """
    
    def __init__(self, fai):
        self.names = []
        self.lengths = {}
        self.offsets = {}
        self.n_bps = {}
        self.n_chars = {}
        for line in fai:
            l = line.strip().split()
            nm = l[0]
            self.names.append(nm)
            self.lengths[nm] = int(l[1])
            self.offsets[nm] = int(l[2])
            self.n_bps[nm] = int(l[3])
            self.n_chars[nm] = int(l[4])
    
    def name(self, i):
        if (i > len(self.names)) or (i < 0):
            raise ValueError('i must be an integer from 0 to number of contigs')
        return self.names[i]
    
    def length(self, seqname):
        if seqname not in self.names:
            raise ValueError('seqname must be in names')
        return self.lengths[seqname]
    
    def offset(self, seqname):
        if seqname not in self.names:
            raise ValueError('seqname must be in names')
        return self.offsets[seqname]
    
    def n_bp(self, seqname):
        if seqname not in self.names:
            raise ValueError('seqname must be in names')
        return self.n_bps[seqname]
    
    def n_char(self, seqname):
        if seqname not in self.names:
            raise ValueError('seqname must be in names')
        return self.n_chars[seqname]
    

class Reader:
    """ 
    Read from an indexed fasta file
    """

    def __init__(self, file):
        """
        Initialize Reader with indexed fasta file

        Attributes
        ----------
        file : string
            File name of fasta file. The index with .fai appended 
            must exist.
        """
        self.fs = open(file, 'rt')
        file_fai = file + ".fai"
        fai_fs = open(file_fai, 'rt')
        l = fai_fs.readlines()
        fai_fs.close()
        self.ix = Index(l)

    def get_seq(self, seq_name, pos, length = 1):
        """
        Get sequence from indexed fasta file
        
        Attributes
        ----------
        seq_name : string
            Name of sequence
        pos : int
            Start position of sequence (0-based)
        length : int
            Length of sequence
        """
        offset = self.ix.offset(seq_name)
        c_l = self.ix.length(seq_name)
        n_bp = self.ix.n_bp(seq_name)
        n_char = self.ix.n_char(seq_name)
        if length <= 0:
            raise ValueError('length must be a positive integer')
        if pos < 0:
            raise ValueError('pos must be a non-negative integer')
        if (pos + length) > c_l:
            raise ValueError('Invalid position and length: too long')
        n_l = n_char * math.floor(pos / n_bp)
        pc = pos % n_bp
        fp = offset + n_l + pc
        g = self.fs.seek(fp)
        seq = ''
        while length > 0:
            bp = self.fs.read(1)
            if bp == '\n':
                continue
            seq = seq + bp
            length = length - 1
        return seq
    
    def close(self):
        self.fs.close()


