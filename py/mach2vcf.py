
import sys
import argparse
import gzip
import pandas as pd
import numpy as np
import datetime
from multiprocessing import Pool
import fasta_fai

#######################################
# Functions
#######################################

def read_dose(fn):
    if fn[-3:] == ".gz":
        dosefs = gzip.open(fn, 'rt')
    else:
        dosefs = open(fn, 'r')
    samples = []
    geno = []
    for line in dosefs:
        linel = line.strip().split('\t')
        id = linel.pop(0)
        typ = linel.pop(0)
        samples.append(id)
        linel = [2 - float(x) for x in linel]
        geno.append(linel)
    return samples, geno


def dose2gt(x, aa = "0/0", ab = "0/1", bb = "1/1", rev = False):
    if x < 0.5:
        if rev:
            return bb
        else:
            return aa
    elif x >= 0.5 and x < 1.5:
        return ab
    else:
        if rev:
            return aa
        else:
            return bb

def dose2gp(x, nr = 6, rev = False):
    gp = [0,0,0]
    if x < 1:
        gp[0] = 1 - x
        gp[1] = x
    elif x >= 1:
        gp[1] = 2 - x
        gp[2] = 1 - gp[1]
    gp = [str(round(i, nr)) for i in gp]
    if rev:
        gp = gp[::-1]
    gpr = ','.join(gp)
    return gpr

def gp2pl(x, nr = 6):
    lim = 10 ** -nr
    gp = [float(i) for i in x.split(',')]
    gp = [lim if i < lim else i for i in gp]
    pl = [-10 * np.log10(i) for i in gp]
    pl = [str(int(round(i))) for i in pl]
    plr = ','.join(pl)
    return plr

def dose2vcf(x, nr = 6, rev = False):
    gt = dose2gt(x, rev = rev)
    gp_pl = dose2gp_pl(x, nr = nr, rev = rev)
    r = ':'.join([gt] + gp_pl)
    return r

def dose2vcf2(x, nr = 4, nl = 6,
        aa = "0/0", ab = "0/1", bb = "1/1"):
    gt = ""
    if x < 0.5:
        gt = aa
    elif x >= 1.5:
        gt = bb
    else:
        gt = ab
    lim = 10 ** -nl
    plim = -10 * np.log10(lim)
    gp = [0,0,0]
    if x < 1:
        gp[0] = 1 - x
        gp[1] = x
    elif x >= 1:
        gp[1] = 2 - x
        gp[2] = 1 - gp[1]
    gp = [round(i, nr) for i in gp]
    pl = [plim if i < lim else -10 * np.log10(i) for i in gp]
    #
    gp = [str(i) for i in gp]
    pl = [str(int(round(i))) for i in pl]
    #
    gpr = ','.join(gp)
    plr = ','.join(pl)
    r = ':'.join([gt, gpr, plr])
    return r

dose2vcf2v = np.vectorize(dose2vcf2)

#######################################
# Parse and set up argument parameters
#######################################

parser = argparse.ArgumentParser()
parser.add_argument("--dose", "-d", help="Mach input dose file", required=True)
parser.add_argument("--info", "-i", help="Mach input INFO file.", required=True)
parser.add_argument("--fasta", "-f", help="Fasta file.", required=True)
parser.add_argument("--fai", help="Fasta index file.", required=True)
parser.add_argument("--fwd_swap", help="Swap ref/alt alleles if necessary to match reference sequence.", action='store_true')
parser.add_argument("--chrm", "-c", help="Chromosome name(s) to include in header. Multiple values should be comma separated.", required=True)
parser.add_argument("--out", "-o", help="Output file name.", required=True)

args = parser.parse_args()

info = pd.read_table(args.info, header = 0, compression = "gzip")
fai = pd.read_table(args.fai, header = None, index_col = 0)
chrms = args.chrm
chrms = chrms.split(',')

today = datetime.date.today()
todaystr = today.strftime("%Y%m%d")

# n_threads = int(args.threads)
n_threads = 1

print("Swapping to match forward strand:", args.fwd_swap)

print("Reading dosages", flush = True)
samples, geno = read_dose(args.dose)
print("Done", flush = True)


geno_a = np.array(geno)
del geno

testing = False
ntest = 1000
if testing:
    geno_a = geno_a[:,:ntest]
    info = info.iloc[:ntest,:]

nv = geno_a.shape[1]
ns = geno_a.shape[0]

#print("Converting dosages")
if n_threads > 1:
    # Parellel version
    if nv < n_threads:
        chunks = 1
    else:
        chunks = round(nv / n_threads)
    def dose2vcf2vw(i):
        return dose2vcf2v(geno_a[:,i])
    pool = Pool(processes = n_threads)
    geno_d = pool.imap(dose2vcf2vw, range(nv), chunks) # return generator
    pool.close()
    pool.join()
    geno_d = [x for x in geno_d] # nv-length list
    geno_d = np.array(geno_d) # nv by ns array
    geno_d = np.transpose(geno_d)
#else:
#    geno_d = dose2vcf2v(geno_a)

#print("Done")

# del geno_a

ofs = open(args.out, 'w')

#######################################
# headers
#######################################

# Headers
hdr_vcf = "##fileformat=VCFv4.3"
hdr_date = '##fileDate=' + todaystr
hdr_gt = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
hdr_gp = '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probability, between 0 and 1">'
hdr_pl = '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">'
info_rsq = '##INFO=<ID=RSQ,Number=1,Type=Float,Description="Imputation R squared, as output by Minimac">'
info_af = '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency for each ALT allele">'
info_lo = '##INFO=<ID=LOOSQ,Number=1,Type=Float,Description="Leave-one-out R squared, as output by Minimac">'
info_geno = '##INFO=<ID=GN,Number=0,Type=Flag,Description="Was the variant genotyped?">'
hdr_contigs = []
for chrm in chrms:
    chrmlen = fai.loc[chrm,1]
    hdr_contigs.append('##contig=<ID=' + chrm + ',length=' + str(chrmlen) + '>')

hdr_contig = '\n'.join(hdr_contigs)

vcf_hdr_l = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
for k in samples:
    vcf_hdr_l.append(k)

vcf_hdr = '\t'.join(vcf_hdr_l)

hdr_l = [hdr_vcf, hdr_date, hdr_gt, hdr_gp, hdr_pl, 
        info_rsq, info_af, info_lo, info_geno, hdr_contig, vcf_hdr]
hdr_all = '\n'.join(hdr_l)

#######################################
# variant records
#######################################

rdr = fasta_fai.Reader(args.fasta)

frmt = 'GT:GP:PL'

t = hdr_all + '\n'
ch = ofs.write(t)

for v in range(nv):
    if v > 0 and v % 1e3 == 0:
        print(v, " variants", flush = True)
    #
    t_name = info.iloc[v, 0]
    t_name_l = t_name.split(":")
    t_chrm = t_name_l[0]
    t_pos = t_name_l[1]
    t_ref = info.iloc[v, 1]
    t_alt = info.iloc[v, 2]
    t_qual = "."
    t_filter = "."
    #
    t_info_l = []
    t_info_l.append("RSQ=" + str(info['Rsq'][v]))
    t_info_l.append("AF=" + str(info['MAF'][v]))
    if info['LooRsq'][v] != "-":
        t_info_l.append("LOOSQ=" + str(info['LooRsq'][v]))
    if info['Genotyped'][v] != "-":
        t_info_l.append("GN")
    t_info = ';'.join(t_info_l)
    #
    #t_gt = geno_d[:,v].tolist()
    # swap alleles
    if args.fwd_swap:
        refseq_base = rdr.get_seq(t_chrm, int(t_pos) - 1, len(t_ref))
        if t_ref == refseq_base:
            pass
        else:
            if t_alt != refseq_base:
                print("Warning:", "no alleles for", t_name, "match the reference sequence. Did not swap.")
            else:
                geno_a[:,v] = 2 - geno_a[:,v]
                reft = t_ref
                t_ref = t_alt
                t_alt = reft
    t_gt = dose2vcf2v(geno_a[:,v])
    t_gt = t_gt.tolist()
    t_name = t_chrm + ":" + t_pos + ":" + t_ref + "_" + t_alt
    #
    t_l = [t_chrm, t_pos, t_name, t_ref, t_alt, t_qual, t_filter, t_info, frmt] + t_gt
    t = '\t'.join(t_l) + '\n'
    ch = ofs.write(t)

ofs.close()

