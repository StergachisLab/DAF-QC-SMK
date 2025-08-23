import pysam
import numpy as np
import pandas as pd
from pyfaidx import Fasta


ct_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align/BRIP1.CT.bam"
ga_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align/BRIP1.GA.bam"
#ct_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/TCF4.25.1.28/align/TCF4.CT.bam"
#ga_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/TCF4.25.1.28/align/TCF4.GA.bam"
ct_bam = pysam.AlignmentFile(ct_bampath, "rb")
ga_bam = pysam.AlignmentFile(ga_bampath, "rb")
#bam_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/align/htt_test.filtered.bam"
#bam = pysam.AlignmentFile(bam_path, "rb")
#seq_metrics = pd.read_csv("/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test.detailed_seq_metrics.reads.tbl.gz", header=0, index_col=0, sep="\t")
fasta_path = "/mmfs1/gscratch/stergachislab/assemblies/simple-names/hg38.fa"
output_bam_path = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align/BRIP1.phased_reads.bam"

fasta = Fasta(fasta_path)


#chrom='chr18'
#start = 55584419
#end = 55587077
chrom = 'chr17'
start = 61860066
end = 61865269-50

threshold = 0.8
het_tolerance = 0.3
#mask_region = "chr18:55586113-55586229"


base_dict = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    'I': 4,
}


ref_matrix = np.zeros(end - start, dtype=int)


for i, base in enumerate(fasta[chrom][start:end].seq):
    if base == 'N':
        ref_matrix[i] = 5
    else:
        ref_matrix[i] = base_dict[base]


#basecounts = np.zeros((10, end - start), dtype=int)


ct_count = ct_bam.count(contig=chrom, start=start, end=end)
ga_count = ga_bam.count(contig=chrom, start=start, end=end)
#print(f"CT reads count: {ct_count}")
#print(f"GA reads count: {ga_count}")

ct_basecounts = np.array(ct_bam.count_coverage(chrom, start, end))/ct_count
ga_basecounts = np.array(ga_bam.count_coverage(chrom, start, end))/ga_count

#ct_basecounts = ct_basecounts/np.sum(ct_basecounts, axis=0, keepdims=True)
#ga_basecounts = ga_basecounts/np.sum(ga_basecounts, axis=0, keepdims=True)

# Identify and ignore reference match postions
a_mask= (ct_basecounts[0] > threshold) & (ga_basecounts[0] > threshold) & (ref_matrix == 0)
c_mask= (ga_basecounts[1] > threshold) & (ref_matrix == 1)
g_mask= (ct_basecounts[2] > threshold) & (ref_matrix == 2)
t_mask= (ct_basecounts[3] > threshold) & (ga_basecounts[3] > threshold) & (ref_matrix == 3)

ref_pos = a_mask | c_mask | g_mask | t_mask
variant_candidates = ~ref_pos

homo_variants = {}
het_variants = {}

# Find homozygous variants
homo_variants["a"] = variant_candidates & (ct_basecounts[0] > threshold) & (ga_basecounts[0] > threshold)
homo_variants["c"] = variant_candidates & (ga_basecounts[1] > threshold)
homo_variants["g"] = variant_candidates & (ct_basecounts[2] > threshold)
homo_variants["t"] = variant_candidates & (ct_basecounts[3] > threshold) & (ga_basecounts[3] > threshold)


# Find heterozygous variants
het_variants['ac']= (variant_candidates &
                    (ga_basecounts[0] + ga_basecounts[1] > threshold) &
                    (ga_basecounts[0] > 0.5 - het_tolerance) &
                    (ga_basecounts[0] < 0.5 + het_tolerance) &
                    (ga_basecounts[1] > 0.5 - het_tolerance) &
                    (ga_basecounts[1] < 0.5 + het_tolerance))

het_variants['ag'] = (variant_candidates & 
                    (ct_basecounts[0] + ct_basecounts[2] > threshold) &
                    (ct_basecounts[0] > 0.5 - het_tolerance) &
                    (ct_basecounts[0] < 0.5 + het_tolerance) &
                    (ct_basecounts[2] > 0.5 - het_tolerance) &
                    (ct_basecounts[2] < 0.5 + het_tolerance))

het_variants['at'] = (variant_candidates &
                    (ct_basecounts[0] + ct_basecounts[3] > threshold) &
                    (ga_basecounts[0] + ga_basecounts[3] > threshold) &
                    (ct_basecounts[0] > 0.5 - het_tolerance) &
                    (ct_basecounts[0] < 0.5 + het_tolerance) &
                    (ct_basecounts[3] > 0.5 - het_tolerance) &
                    (ct_basecounts[3] < 0.5 + het_tolerance) &
                    (ga_basecounts[0] > 0.5 - het_tolerance) &
                    (ga_basecounts[0] < 0.5 + het_tolerance) &
                    (ga_basecounts[3] > 0.5 - het_tolerance) &
                    (ga_basecounts[3] < 0.5 + het_tolerance))

het_variants['cg'] = (variant_candidates &
                    (ct_basecounts[1] + ct_basecounts[2] + ct_basecounts[3] > threshold) &
                    (ga_basecounts[0] + ga_basecounts[1] + ga_basecounts[2] > threshold) &
                    (ct_basecounts[2] > 0.5 - het_tolerance) &
                    (ct_basecounts[2] < 0.5 + het_tolerance) &
                    (ga_basecounts[1] > 0.5 - het_tolerance) &
                    (ga_basecounts[1] < 0.5 + het_tolerance))

het_variants['ct'] = (variant_candidates &
                    (ga_basecounts[1] + ga_basecounts[3] > threshold) &
                    (ga_basecounts[1] > 0.5 - het_tolerance) &
                    (ga_basecounts[1] < 0.5 + het_tolerance) &
                    (ga_basecounts[3] > 0.5 - het_tolerance) &
                    (ga_basecounts[3] < 0.5 + het_tolerance))

het_variants['gt'] = (variant_candidates &
                    (ct_basecounts[2] + ct_basecounts[3] > threshold) &
                    (ct_basecounts[2] > 0.5 - het_tolerance) &
                    (ct_basecounts[2] < 0.5 + het_tolerance) &
                    (ct_basecounts[3] > 0.5 - het_tolerance) &
                    (ct_basecounts[3] < 0.5 + het_tolerance))


# Remove values in het_variants dictionary that are all False
het_variants = {k: v for k, v in het_variants.items() if np.any(v)}
het_variant_list = []
for key in het_variants:
    coordinates= np.array(np.where(het_variants[key])) + start
    coordinates_list = coordinates[0].tolist()
    coordinates_list = [(coordinate, key) for coordinate in coordinates_list]
    het_variant_list.extend(coordinates_list)

het_variant_list.sort(key=lambda x: x[0])
het_variant_pos = [pos[0] for pos in het_variant_list]
het_variant_types = [pos[1] for pos in het_variant_list]
#het_variant_pos = np.array(np.where(np.any(list(het_variants.values()), axis=0))) + start


bases=[]
read_names = []

for fiber in ga_bam.fetch(chrom, start, end):
    pairs = fiber.get_aligned_pairs(with_seq=True)
    fiber_sequence = fiber.query_sequence
    fiber_bases = []
    for query_idx, ref_idx, ref_base in pairs:
        if ref_idx is not None and ref_idx in het_variant_pos:
            if query_idx is not None:
                fiber_bases.append(fiber_sequence[query_idx])
            else:
                fiber_bases.append(np.nan)
    bases.append(fiber_bases)
    read_names.append(fiber.query_name)





ga_variant_df = pd.DataFrame(bases, index=read_names, columns=het_variant_pos)

#haplotype_counts = variant_df[het_variant_pos[0][1:]].value_counts()




ga_variant_df_modified = ga_variant_df.copy()
for col in ga_variant_df_modified.columns:
    variant_type = het_variant_types[het_variant_pos.index(col)]
#    if 'c' in variant_type and variant_type != 'ct':
#        ga_variant_df_modified[col] = ga_variant_df_modified[col].replace({'T':'Y', 'C':'Y'})
    if 'g' in variant_type and variant_type != 'ag':
        ga_variant_df_modified[col] = ga_variant_df_modified[col].replace({'A':'R', 'G':'R'})
    elif variant_type == 'ag':
        ga_variant_df_modified[col] = np.nan

ga_grouped = ga_variant_df_modified.groupby(het_variant_pos)
ga_haplotype_indices = ga_grouped.groups
ga_haplotype_counts = {key:len(ga_haplotype_indices[key]) for key in ga_haplotype_indices}

ga_hap1, ga_hap2 = sorted(ga_haplotype_counts, key=ga_haplotype_counts.get, reverse=True)[:2]
ga_hap1, ga_hap2 = list(ga_hap1), list(ga_hap2)
print("Top GA haplotypes:", ga_hap1, ga_hap2)
for h1,h2 in zip(ga_hap1, ga_hap2):
    if h1 == h2 and not pd.isna(h1):
        print("Conflict at position:", h1)
        ga_hap1 = np.nan
        ga_hap2 = np.nan
ga_hap_labels = {tuple(ga_hap1):1, tuple(ga_hap2):2}



for fiber in ct_bam.fetch(chrom, start, end):
    pairs = fiber.get_aligned_pairs(with_seq=True)
    fiber_sequence = fiber.query_sequence
    fiber_bases = []
    for query_idx, ref_idx, ref_base in pairs:
        if ref_idx is not None and ref_idx in het_variant_pos:
            if query_idx is not None:
                fiber_bases.append(fiber_sequence[query_idx])
            else:
                fiber_bases.append(np.nan)
    bases.append(fiber_bases)
    read_names.append(fiber.query_name)

ct_variant_df = pd.DataFrame(bases, index=read_names, columns=het_variant_pos)
ct_variant_df_modified = ct_variant_df.copy()
for col in ct_variant_df_modified.columns:
    variant_type = het_variant_types[het_variant_pos.index(col)]
    if 'c' in variant_type and variant_type != 'ct':
        ct_variant_df_modified[col] = ct_variant_df_modified[col].replace({'T':'Y', 'C':'Y'})
#    elif 'g' in variant_type and variant_type != 'ag':
#        ct_variant_df_modified[col] = ct_variant_df_modified[col].replace({'A':'R', 'G':'R'})
    elif variant_type == 'ct':
        ct_variant_df_modified[col] = np.nan

ct_grouped = ct_variant_df_modified.groupby(het_variant_pos)
ct_haplotype_indices = ct_grouped.groups
ct_haplotype_counts = {key:len(ct_haplotype_indices[key]) for key in ct_haplotype_indices}

ct_hap1, ct_hap2 = sorted(ct_haplotype_counts, key=ct_haplotype_counts.get, reverse=True)[:2]
ct_hap1, ct_hap2 = list(ct_hap1), list(ct_hap2)
print("Top CT haplotypes:", ct_hap1, ct_hap2)
for h1,h2 in zip(ct_hap1, ct_hap2):
    if h1 == h2 and not pd.isna(h1):
        print("Conflict at position:", h1)
        ct_hap1 = np.nan
        ct_hap2 = np.nan


ct_hap_labels = {tuple(ct_hap1):3, tuple(ct_hap2):4}


# Combine GA and CT haplotypes
ga_hap1_corrected = ['G' if a == 'R' else a for a in ga_hap1]
ga_hap2_corrected = ['G' if a == 'R' else a for a in ga_hap2]
ct_hap1_corrected = ['C' if a == 'Y' else a for a in ct_hap1]
ct_hap2_corrected = ['C' if a == 'Y' else a for a in ct_hap2]
hap_match_scores = {'h11':sum((a==b) and not pd.isna(a) and not pd.isna(b) for a,b in zip(ga_hap1_corrected , ct_hap1_corrected)),
                    'h12':sum((a==b) and not pd.isna(a) and not pd.isna(b) for a,b in zip(ga_hap1_corrected , ct_hap2_corrected))
                    }

if hap_match_scores['h11'] == 0 and hap_match_scores['h12'] == 0 or hap_match_scores['h11'] == hap_match_scores['h12']:
    combined_haplotyping = False
else:
    combined_haplotyping = True
    best_match = max(hap_match_scores, key=hap_match_scores.get)
    if best_match == 'h11':
        final_hap1 = [a if not pd.isna(a) else b for a,b in zip(ga_hap1_corrected , ct_hap1_corrected)]
#        final_hap1 = ['C' if a == 'Y' else a for a in final_hap1]
#        final_hap1 = ['G' if a == 'R' else a for a in final_hap1]
        final_hap2 = [a if not pd.isna(a) else b for a,b in zip(ga_hap2_corrected , ct_hap2_corrected)]
#        final_hap2 = ['C' if a == 'Y' else a for a in final_hap2]
#        final_hap2 = ['G' if a == 'R' else a for a in final_hap2]
        ct_hap_labels = {tuple(ct_hap1):1, tuple(ct_hap2):2}

    else:
        final_hap1 = [a if not pd.isna(a) else b for a,b in zip(ga_hap1_corrected, ct_hap2_corrected)]
#        final_hap1 = ['C' if a == 'Y' else a for a in final_hap1]
#        final_hap1 = ['G' if a == 'R' else a for a in final_hap1]
        final_hap2 = [a if not pd.isna(a) else b for a,b in zip(ga_hap2_corrected, ct_hap1_corrected)]
#        final_hap2 = ['C' if a == 'Y' else a for a in final_hap2]
#        final_hap2 = ['G' if a == 'R' else a for a in final_hap2]
        ct_hap_labels = {tuple(ct_hap1):2, tuple(ct_hap2):1}


print("Variant positions:", het_variant_pos)
print("Variant types:", het_variant_types)
print("Final haplotypes:", final_hap1, final_hap2)


# Write phased reads to new BAM file
ga_hap1_ids = {item:ga_hap_labels[tuple(ga_hap1)] for item in ga_haplotype_indices[tuple(ga_hap1)]}
ga_hap2_ids = {item:ga_hap_labels[tuple(ga_hap2)] for item in ga_haplotype_indices[tuple(ga_hap2)]}
ct_hap1_ids = {item:ct_hap_labels[tuple(ct_hap1)] for item in ct_haplotype_indices[tuple(ct_hap1)]}
ct_hap2_ids = {item:ct_hap_labels[tuple(ct_hap2)] for item in ct_haplotype_indices[tuple(ct_hap2)]}
ga_haplotype_indices = {**ga_hap1_ids, **ga_hap2_ids}
ct_haplotype_indices = {**ct_hap1_ids, **ct_hap2_ids}

with pysam.AlignmentFile(output_bam_path, "wb", template=ct_bam) as out_bam:
    for fiber in ga_bam.fetch(chrom, start, end):
        if fiber.query_name in ga_haplotype_indices:
            fiber.set_tag('HP', ga_haplotype_indices[fiber.query_name], value_type='i')
        out_bam.write(fiber)
    for fiber in ct_bam.fetch(chrom, start, end):
        if fiber.query_name in ct_haplotype_indices:
            fiber.set_tag('HP', ct_haplotype_indices[fiber.query_name], value_type='i')
        out_bam.write(fiber)
ct_bam.close()
ga_bam.close()

    




#for key in het_variants:
#    coordinates = np.array(np.where(het_variants[key])) + start
#    for coordinate in coordinates[0]:
#        fibers = ct_bam.pileup(chrom, coordinate, coordinate+1)
#        for fiber in fibers:


    



# phase top and bottom strands separately, then reconcile phasing when possible
# Do this when a variant can be reconciled in both strands, can consider fraction of
# aminated and deaminated type base when another variant is present on the same strand



#het_variants = a_c_heterozygous | a_g_heterozygous | a_t_heterozygous | c_g_heterozygous | c_t_heterozygous | g_t_heterozygous

#potential_indels = (np.sum(ct_basecounts, axis=0) < 0.9) | (np.sum(ga_basecounts, axis=0) < 0.9)



'''

filebase = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align"
for haplotype, index in haplotype_indices.items():
    file_name = f"{filebase}/" + "".join(map(str, haplotype)) + ".GAzids2.txt"
    with open(file_name, 'w') as f:
        for read_name in index:
            f.write(f"{read_name}\n")



#coordinates using pyft

#ct_ft = pyft.Fiberbam(ct_bampath)
ga_ft = pyft.Fiberbam(ga_bampath)
#iteration = 0
bases=[]
read_names = []
for fiber in ga_ft:
#    if iteration < 10:
    if fiber.qname == "4239_LHL4ZV_3":
        var_coordinates = fiber.lift_query_positions(het_variant_pos)
        print("var_coords:", var_coordinates)
        var_coordinates_check = fiber.lift_query_positions([x+1 for x in het_variant_pos])
        print("var_coords+1:", var_coordinates_check)

        fiber_sequence = fiber.seq
        fiber_bases = [fiber_sequence[x] for x in var_coordinates if x < len(fiber_sequence)]
        bases.append(fiber_bases)
        read_names.append(fiber.qname)
#        print(fiber_bases)
#        print(var_coordinates)
#        iteration += 1
#    else:
#        break



for coord in bam.pileup(contig=chrom, start=start, end=end, truncate=True, max_depth=2000):
#    print(f"Position: {coord.pos}")
    for read in coord.pileups:
        readname = read.alignment.query_name
        if readname not in seq_metrics.index:
            continue
        strand = seq_metrics.loc[readname, 'strand']
        idx = 0 if strand == 'CT' else 5 if strand == 'GA' else None
        if idx is None:
            continue
        if read.is_del or read.is_refskip:
            base= 'I'
        else:
            base = read.alignment.query_sequence[read.query_position]
        pos_idx = coord.pos - start
        idx += base_dict[base]
        if 0 <= pos_idx < (end - start):
            basecounts[idx, pos_idx] += 1

bam.close()

#print("Base counts matrix:")
#print(basecounts)


# create a mask to apply to positions where ref base is the most common base
'''