
import pysam
import fasta_fai
import gzip

STRAND = ("AMB", "FWD", "RVR")
REFSEQ = ("NA", "REF", "ALT")
VTYPE = ("NA", "SNP", "INS", "DEL")

COMPLS = {'A' : 'T', 
        'T' : 'A', 
        'G' : 'C', 
        'C' : 'G', 
        '.' : '.'}

def compl(seq, reverse = False):
    if reverse:
        seq = seq[::-1]
    ret = [COMPLS[x] for x in seq]
    ret = ''.join(ret)
    return ret


def is_ambg(ref, alt):
    if ref == 'A' and alt == 'T':
        return True
    if ref == 'T' and alt == 'A':
        return True
    if ref == 'G' and alt == 'C':
        return True
    if ref == 'C' and alt == 'G':
        return True
    return False

def check_ref_strand(vcf_file, fasta_file, out_file, verbose = False):
    """
    Check strand, reference allele, and strand ambiguity of variants 
    in a vcf file. For each variant, outputs:
        - if strand is forward, reverse, or ambiguous
        - if reference seq allele is reference or alternate
        - if SNP, insertion, or deletion
    
    Attributes
    ----------
    vcf_file : str
        VCF/BCF file name.
    fasta_file : str
        Fasta file name. Must be unzipped and have an fai index.
    out_file : str
        Output file name.
    """
    rdr = fasta_fai.Reader(fasta_file)
    vcf_in = pysam.VariantFile(vcf_file, mode = 'r')
    
    if out_file[-3:] == ".gz":
        out_fs = gzip.open(out_file, mode = 'wt')
    else:
        out_fs = open(out_file, mode = 'w')
    
    outl = ["CHR", "POS", "ID", "REF", "ALT", "STRAND", "REFSEQ", "TYPE"]
    out = '\t'.join(outl) + '\n'
    out_fs.write(out)
    counter = 0
    for rec in vcf_in.fetch():
        counter = counter + 1
        if verbose and counter % 1000 == 0:
            print(counter, "records")
        vid = rec.id
        chr_name = rec.chrom
        ref = rec.alleles[0]
        if (len(rec.alleles) == 1):
            alt = '.'
        else:
            alt = rec.alleles[1]
        pos = rec.pos - 1
        refseq_base = rdr.get_seq(chr_name, pos, len(ref))
        strnd = 0
        refseq = 0
        vtype = 0
        if ref == refseq_base:
            strnd = 1
            refseq = 1
        elif alt == refseq_base:
            strnd = 1
            refseq = 2
        else:
            refc = compl(ref, reverse = True)
            altc = compl(alt, reverse = True)
            if refc == refseq_base:
                strnd = 2
                refseq = 1
            elif altc == refseq_base:
                strnd = 2
                refseq = 2
        if is_ambg(ref, alt):
            strnd = 0
        if len(ref) == 1 and len(alt) == 1:
            vtype = 1

        if len(ref) > 1 and len(alt) == 1:
            vtype = 2
        if len(ref) == 1 and len(alt) > 1:
            vtype = 3
        
        outl = [chr_name, str(pos+1), vid, ref, alt, STRAND[strnd], REFSEQ[refseq], VTYPE[vtype]]
        out = '\t'.join(outl) + '\n'
        out_fs.write(out)
    rdr.close()
    vcf_in.close()
    out_fs.close()
    if verbose:
        print("Finished", counter, "records")


def get_refseq_allele(vcf_file, fasta_file, out_file, verbose = False):
    """
    Output a reference VCF file that contains the ref/alt allele coding 
    with the ref allele matching the reference sequence. This assumes that 
    all alleles are in the forward strand. If neither the ref or alt 
    allele is not found in the reference sequence, the alleles are output 
    as missing.
    
    Attributes
    ----------
    vcf_file : str
        VCF/BCF file name.
    fasta_file : str
        Fasta file name. Must be unzipped and have an fai index.
    out_file : str
        Output VCF file name.
    """

    rdr = fasta_fai.Reader(fasta_file)
    vcf_in = pysam.VariantFile(vcf_file, mode = 'r')
    
    vcf_out = pysam.VariantFile(out_file, mode = 'w')
    for r in vcf_in.header.records:
        vcf_out.header.add_record(r)
    
    counter = 0
    for rec in vcf_in.fetch():
        counter = counter + 1
        if verbose and counter % 1000 == 0:
            print(counter, "records")
        rec_out = vcf_out.new_record()
        rec_out.id = rec.id
        rec_out.pos = rec.pos
        rec_out.chrom = rec.chrom
        o = list(range(len(rec.alleles)))
        orv = list(reversed(o))
        ref_i = 0
        for i in orv:
            a = rec.alleles[i]
            refseq_base = rdr.get_seq(rec.chrom, rec.pos-1, len(a))
            if a == refseq_base:
                ref_i = i
        o[ref_i] = 0
        o[0] = ref_i
        alleles = list()
        for i in o:
            alleles.append( rec.alleles[i] )
        if (len(alleles) == 1):
            alleles.append('.')
        rec_out.alleles = tuple(alleles)
        alleles_set = set(rec.alleles)
        vcf_out.write(rec_out)
    rdr.close()
    vcf_in.close()
    vcf_out.close()
    if out_file[-3:] == ".gz":
        pysam.tabix_index(out_file, preset = "vcf", force = True)

def garb():
        alleles = [rec.alleles]
        ref_out = "."
        alt_l = []
        alt_out = "."
        pos = rec.pos - 1
        for a in rec.alleles:
            refseq_base = rdr.get_seq(chr_name, pos, len(a))
            if a == refseq_base:
                ref_out = a
                alt_set = alleles_set.difference(ref_out)
                if len(alt_set) > 0:
                    alt_out = ','.join(list(alt_set))
                else:
                    alt_out = '.'
                break

        outl = [chr_name, str(pos+1), vid, ref_out, alt_out, "100", ".", "."]
        out = '\t'.join(outl) + '\n'
        out_fs.write(out)

def garb2():
    rdr.close()
    vcf_in.close()
    out_fs.close()
    if verbose:
        print("Finished", counter, "records")

