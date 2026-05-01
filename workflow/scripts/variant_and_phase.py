import pysam
import numpy as np
import pandas as pd
from pyfaidx import Fasta


#ct_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align/BRIP1.CT.bam"
#ga_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align/BRIP1.GA.bam"
#output_bam_path = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align/BRIP1.phased_reads.bam"
#chrom = 'chr17'
#start = 61860066
#end = 61865269-50

#threshold = 0.8
#het_tolerance = 0.3

#ct_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/TCF4.25.1.28/align/TCF4.CT.bam"
#ga_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/TCF4.25.1.28/align/TCF4.GA.bam"
#output_bam_path = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/TCF4.25.1.28/align/TCF4.phased_reads.bam"
#chrom='chr18'
#start = 55584419
#end = 55587077

#threshold = 0.9
#het_tolerance = 0.2


ct_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/NAPA/align/NAPA.CT.bam"
ga_bampath = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/NAPA/align/NAPA.GA.bam"
output_bam_path = "/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/NAPA/align/NAPA.phased_reads.bam"
chrom = "chr19"
start = 47514488
end = 47518985

threshold = 0.8 
het_tolerance = 0.2


#ct_bampath = 
#ga_bampath = 
#output_bam_path = 
#chrom = 
#start = 
#end =

#threshold = 
#het_tolerance =

fasta_path = "/mmfs1/gscratch/stergachislab/assemblies/simple-names/hg38.fa"




ct_bam = pysam.AlignmentFile(ct_bampath, "rb")
ga_bam = pysam.AlignmentFile(ga_bampath, "rb")
fasta = Fasta(fasta_path)




#ignore_region = "chr18:55586113-55586229"


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



ct_count = ct_bam.count(contig=chrom, start=start, end=end)
ga_count = ga_bam.count(contig=chrom, start=start, end=end)
all_count = ct_count + ga_count


ct_basecounts = np.array(ct_bam.count_coverage(chrom, start, end, quality_threshold=0))/ct_count
ga_basecounts = np.array(ga_bam.count_coverage(chrom, start, end, quality_threshold=0))/ga_count

#ct_basecounts = np.array(ct_bam.count_coverage(chrom, start, end))
#ga_basecounts = np.array(ga_bam.count_coverage(chrom, start, end))

#ct_basecounts = ct_basecounts/np.sum(ct_basecounts, axis=0, keepdims=True)
#ga_basecounts = ga_basecounts/np.sum(ga_basecounts, axis=0, keepdims=True)

# Identify and ignore reference match postions. Use only CT or GA when ref base on opposite strand could be deaminated
a_mask= (ref_matrix == 0) & (ct_basecounts[0] > threshold) & (ga_basecounts[0] > threshold)
c_mask= (ref_matrix == 1) & (ga_basecounts[1] > threshold)
g_mask= (ref_matrix == 2) & (ct_basecounts[2] > threshold)
t_mask= (ref_matrix == 3) & (ct_basecounts[3] > threshold) & (ga_basecounts[3] > threshold)

ref_pos = a_mask | c_mask | g_mask | t_mask
variant_candidates = ~ref_pos

# create a matrix with the first row as ref positions, and rows underneath consisting of base counts from ct and ga at each position, only at variant candidate positions
#variant_matrix = np.vstack(np.array(enumerate(0,len(ref_matrix), ref_matrix, ct_basecounts, ga_basecounts])[:,variant_candidates]




homo_variants = {}
het_variants = {}

# Find homozygous variants
homo_variants["a"] = variant_candidates & (ct_basecounts[0] > threshold) & (ga_basecounts[0] > threshold)
homo_variants["c"] = variant_candidates & (ga_basecounts[1] > threshold)
homo_variants["g"] = variant_candidates & (ct_basecounts[2] > threshold)
homo_variants["t"] = variant_candidates & (ct_basecounts[3] > threshold) & (ga_basecounts[3] > threshold)

homo_variant_list = [
(coord + start, k) 
    for k,v in homo_variants.items() 
    if np.any(v)
    for coord in np.where(v)[0]
]

homo_variant_list.sort(key=lambda x: x[0])
homo_variant_pos = [v.item() for v,_ in homo_variant_list]
homo_variant_types = [t for _,t in homo_variant_list]

# Find heterozygous variants
het_variants['ac']= (variant_candidates &
                    (ga_basecounts[0:2,:].sum(axis=0) > threshold) &
                    (ga_basecounts[0] > 0.5 - het_tolerance) &
                    (ga_basecounts[0] < 0.5 + het_tolerance) &
                    (ga_basecounts[1] > 0.5 - het_tolerance) &
                    (ga_basecounts[1] < 0.5 + het_tolerance) &
                    (ct_basecounts[0] > 0.5 - het_tolerance) &
                    (ct_basecounts[0] < 0.5 + het_tolerance) )

het_variants['ag'] = (variant_candidates & 
                    (ct_basecounts[[0,2],:].sum(axis=0) > threshold) &
                    (ct_basecounts[0] > 0.5 - het_tolerance) &
                    (ct_basecounts[0] < 0.5 + het_tolerance) &
                    (ct_basecounts[2] > 0.5 - het_tolerance) &
                    (ct_basecounts[2] < 0.5 + het_tolerance))

het_variants['at'] = (variant_candidates &
                    (ct_basecounts[[0,3],:].sum(axis=0) > threshold) &
                    (ga_basecounts[[0,3],:].sum(axis=0) > threshold) &
                    (ct_basecounts[0] > 0.5 - het_tolerance) &
                    (ct_basecounts[0] < 0.5 + het_tolerance) &
                    (ct_basecounts[3] > 0.5 - het_tolerance) &
                    (ct_basecounts[3] < 0.5 + het_tolerance) &
                    (ga_basecounts[0] > 0.5 - het_tolerance) &
                    (ga_basecounts[0] < 0.5 + het_tolerance) &
                    (ga_basecounts[3] > 0.5 - het_tolerance) &
                    (ga_basecounts[3] < 0.5 + het_tolerance))

het_variants['cg'] = (variant_candidates &
                    (ct_basecounts[1:4,:].sum(axis=0) > threshold) &
                    (ga_basecounts[:3,:].sum(axis=0) > threshold) &
                    (ct_basecounts[2] > 0.5 - het_tolerance) &
                    (ct_basecounts[2] < 0.5 + het_tolerance) &
                    (ga_basecounts[1] > 0.5 - het_tolerance) &
                    (ga_basecounts[1] < 0.5 + het_tolerance))

het_variants['ct'] = (variant_candidates &
                    (ga_basecounts[[1,3],:].sum(axis=0) > threshold) &
                    (ga_basecounts[1] > 0.5 - het_tolerance) &
                    (ga_basecounts[1] < 0.5 + het_tolerance) &
                    (ga_basecounts[3] > 0.5 - het_tolerance) &
                    (ga_basecounts[3] < 0.5 + het_tolerance))

het_variants['gt'] = (variant_candidates &
                    (ct_basecounts[[2,3],:].sum(axis=0) > threshold) &
                    (ct_basecounts[2] > 0.5 - het_tolerance) &
                    (ct_basecounts[2] < 0.5 + het_tolerance) &
                    (ct_basecounts[3] > 0.5 - het_tolerance) &
                    (ct_basecounts[3] < 0.5 + het_tolerance) &
                    (ga_basecounts[3] > 0.5 - het_tolerance) &
                    (ga_basecounts[3] < 0.5 + het_tolerance))


# Remove values in het_variants dictionary that are all False
het_variant_list = [
    (coord + start, k) 
    for k,v in het_variants.items() 
    if np.any(v)
    for coord in np.where(v)[0]
]

het_variant_list.sort(key=lambda x: x[0])
het_variant_pos = [v.item() for v,_ in het_variant_list]
het_variant_types = [t for _,t in het_variant_list]
#het_variant_pos, het_variant_types = zip(*het_variant_list) if het_variant_list else ([], []) 
#het_variant_types = list(het_variant_types)



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
#ga_hap1, ga_hap2 = list(ga_hap1), list(ga_hap2)
#print("Top GA haplotypes:", ga_hap1, ga_hap2)
for h1,h2 in zip(ga_hap1, ga_hap2):
    if h1 == h2 and not pd.isna(h1):
        print("Conflict at position:", h1)
        ga_hap1 = np.nan
        ga_hap2 = np.nan
ga_hap_labels = {ga_hap1:1, ga_hap2:2}



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
#ct_hap1, ct_hap2 = list(ct_hap1), list(ct_hap2)
#print("Top CT haplotypes:", ct_hap1, ct_hap2)
for h1,h2 in zip(ct_hap1, ct_hap2):
    if h1 == h2 and not pd.isna(h1):
        print("Conflict at position:", h1)
        ct_hap1 = np.nan
        ct_hap2 = np.nan


ct_hap_labels = {ct_hap1:3, ct_hap2:4}


# Combine GA and CT haplotypes
ga_hap1_corrected = ['G' if a == 'R' else a for a in ga_hap1]
ga_hap2_corrected = ['G' if a == 'R' else a for a in ga_hap2]
ct_hap1_corrected = ['C' if a == 'Y' else a for a in ct_hap1]
ct_hap2_corrected = ['C' if a == 'Y' else a for a in ct_hap2]
hap_match_scores = {'h11':sum((a==b) and not pd.isna(a) and not pd.isna(b) for a,b in zip(ga_hap1_corrected , ct_hap1_corrected)),
                    'h12':sum((a==b) and not pd.isna(a) and not pd.isna(b) for a,b in zip(ga_hap1_corrected , ct_hap2_corrected))
                    }


final_haplotypes = []
if hap_match_scores['h11'] == 0 and hap_match_scores['h12'] == 0 or hap_match_scores['h11'] == hap_match_scores['h12']:
    combined_haplotyping = False
    final_haplotypes.append(ga_hap1_corrected)
    final_haplotypes.append(ga_hap2_corrected)
    final_haplotypes.append(ct_hap1_corrected)
    final_haplotypes.append(ct_hap2_corrected)
else:
    combined_haplotyping = True
    best_match = max(hap_match_scores, key=hap_match_scores.get)
    if best_match == 'h11':
        final_haplotypes.append(a if not pd.isna(a) else b for a,b in zip(ga_hap1_corrected , ct_hap1_corrected))
#        final_hap1 = ['C' if a == 'Y' else a for a in final_hap1]
#        final_hap1 = ['G' if a == 'R' else a for a in final_hap1]
        final_haplotypes.append(a if not pd.isna(a) else b for a,b in zip(ga_hap2_corrected , ct_hap2_corrected))
#        final_hap2 = ['C' if a == 'Y' else a for a in final_hap2]
#        final_hap2 = ['G' if a == 'R' else a for a in final_hap2]
        ct_hap_labels = {ct_hap1:1, ct_hap2:2}

    else:
        final_haplotypes.append(a if not pd.isna(a) else b for a,b in zip(ga_hap1_corrected, ct_hap2_corrected))
        final_haplotypes.append(a if not pd.isna(a) else b for a,b in zip(ga_hap2_corrected, ct_hap1_corrected))
#        final_hap1 = [a if not pd.isna(a) else b for a,b in zip(ga_hap1_corrected, ct_hap2_corrected)]
#        final_hap1 = ['C' if a == 'Y' else a for a in final_hap1]
#        final_hap1 = ['G' if a == 'R' else a for a in final_hap1]
#        final_hap2 = [a if not pd.isna(a) else b for a,b in zip(ga_hap2_corrected, ct_hap1_corrected)]
#        final_hap2 = ['C' if a == 'Y' else a for a in final_hap2]
#        final_hap2 = ['G' if a == 'R' else a for a in final_hap2]
        ct_hap_labels = {ct_hap1:2, ct_hap2:1}





# Write phased reads to new BAM file
ga_hap1_ids = {item:ga_hap_labels[ga_hap1] for item in ga_haplotype_indices[ga_hap1]}
ga_hap2_ids = {item:ga_hap_labels[ga_hap2] for item in ga_haplotype_indices[ga_hap2]}
ct_hap1_ids = {item:ct_hap_labels[ct_hap1] for item in ct_haplotype_indices[ct_hap1]}
ct_hap2_ids = {item:ct_hap_labels[ct_hap2] for item in ct_haplotype_indices[ct_hap2]}
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



print("Homo variant positions:", homo_variant_pos)
print("Homo variant types:", homo_variant_types)
print("Het variant positions:", het_variant_pos)
print("Het variant types:", het_variant_types)
#print("Final haplotypes:", final_hap1, final_hap2)
print("Final haplotypes:", ["".join(map(str, hap)) for hap in final_haplotypes])

'''
filebase = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test_data/results/LHL4ZV_3_BRIP1-C5-F.fastq/LHL4ZV_3_BRIP1-C5-F.fastq/align"
for haplotype, index in haplotype_indices.items():
    file_name = f"{filebase}/" + "".join(map(str, haplotype)) + ".GAzids2.txt"
    with open(file_name, 'w') as f:
        for read_name in index:
            f.write(f"{read_name}\n")


'''