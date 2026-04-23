import sys
import argparse

parser = argparse.ArgumentParser(prog='prepare-vcf.py', description="cat <vcf-file> | python3 prepare-vcf.py ")
parser.add_argument('--missing', metavar='MISSING', type=float, default=0.0,
                    help="Maximum allowed fraction of missing alleles per position (applied to on-reference contigs).")
parser.add_argument('--alt-missing', metavar='MISSING', type=float, default=None,
                    help="Maximum allowed fraction of missing alleles per position for *_alt contigs. "
                         "Default: same as --missing. Use a higher value (e.g. 0.9) to keep "
                         "off-reference variants that are sparse by construction (present in "
                         "only a few panel haplotypes) while still dropping the 'literally "
                         "empty across every haplotype' long tail that explodes PanGenie's "
                         "kmer tables.")
parser.add_argument('--alt-suffix', metavar='SUFFIX', default='_alt',
                    help="CHROM suffix marking off-reference contigs that get the --alt-missing "
                         "threshold instead of --missing (default: _alt). Empty string disables the "
                         "alt-vs-onref distinction (both use --missing).")
args = parser.parse_args()

alt_missing_threshold = args.alt_missing if args.alt_missing is not None else args.missing

total_records = 0
total_alleles = 0
missing_records = 0
missing_alleles = 0
ns_records = 0
ns_alleles = 0
written_records = 0
written_alleles = 0
alt_exempt_records = 0

for line in sys.stdin:
	if line.startswith('#'):
		# store header lines
		print(line.strip())
		continue
	total_records += 1
	fields = line.strip().split()
	n_alt_alleles = len(fields[4].split(','))
	total_alleles += n_alt_alleles
	if any([c not in 'CAGTcagt,' for c in fields[4]]):
		ns_records += 1
		ns_alleles += n_alt_alleles
		continue
	n_paths = len(fields[9:]) * 2
	n_missing = 0
	n_total = 0
	for gt in fields[9:]:
		n_total += 2
		if gt == './.':
			n_missing += 2
		elif gt == '.':
			n_missing += 2
		elif '|' in gt:
			alleles = gt.split('|')
			n_missing += alleles.count('.')
		else:
			raise Exception('VCF contains unphased positions.')
	frac_missing = n_missing / n_total
	is_alt_contig = args.alt_suffix != '' and fields[0].endswith(args.alt_suffix)
	threshold = alt_missing_threshold if is_alt_contig else args.missing
	if frac_missing > threshold:
		missing_records += 1
		missing_alleles += n_alt_alleles
		continue
	# Track how many alt records were retained only because the looser alt
	# threshold applied (useful sanity-check signal in the stderr summary).
	if is_alt_contig and frac_missing > args.missing:
		alt_exempt_records += 1
	print(line.strip())
	written_records += 1
	written_alleles += n_alt_alleles

# print statistics
sys.stderr.write('skipped ' + str(missing_records) + ' (' + str(missing_alleles) + ') records (alleles) for which fraction of missing alleles exceeds threshold.\n')
sys.stderr.write('kept ' + str(alt_exempt_records) + ' records on *' + args.alt_suffix +
                 ' contigs that would have been filtered by --missing (alt threshold='
                 + str(alt_missing_threshold) + ').\n')
sys.stderr.write('skipped ' + str(ns_records) + ' (' + str(ns_alleles) + ') records (alleles) for which alternative alleles contained Ns\n.')
sys.stderr.write('kept ' + str(written_records) + ' (' + str(written_alleles) + ') records (alleles) of ' + str(total_records) + ' (' + str(total_alleles) + ').\n')
