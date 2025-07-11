import pysam
import pandas as pd
import pyabpoa as pa


bam_file = snakemake.input.bam
output_bam = snakemake.output.bam
min_read_count = snakemake.params.consensus_min_reads
print('min_read_count:', min_read_count)

#bam_file = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/temp/htt_test/align/htt_test.filtered.bam"
#output_bam = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/htt_test.consensus.bam"
#min_read_count = 3 



def bam_to_dfm(bam_file):
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    reads = []
    for read in samfile.fetch():
        if read.has_tag('du'):
            reads.append({
                'qname': read.query_name,
                'start': read.reference_start,
                'end': read.reference_end,
                'seq': read.query_sequence,
                'du': read.get_tag('du')
            })
    samfile.close()
    dfm = pd.DataFrame(reads)

    return dfm


def abpoa_MSA(seq_list):
    # abPOA implementation. 
    '''Takes list of sequences as input and returns abPOA consensus and MSA.
    Uses abPOA built-in consensus generator'''
    aligner = pa.msa_aligner()
    alignment = aligner.msa(seq_list, out_cons=True, out_msa=True)

#    for seq in alignment.cons_seq:
#        print(seq)
    consensus = alignment.cons_seq[0]

    msa_alignment = alignment.msa_seq


    return consensus, msa_alignment


def consensus_dfm(dfm, min_read_count=3):
    # TODO: determine what minimum read count should be
    # TODO: subsample max read count

    du_names=dfm['du'].unique()[1:]
    dups_consensus=[]

    for du in du_names:
        dups_df=dfm[dfm['du']==du]
        dup_read_count=len(dups_df)
        if dup_read_count < min_read_count:
            continue
    
        sequences=dups_df['seq'].tolist()

        consensus, rep_MSA = abpoa_MSA(sequences)
        dups_consensus.append({'du': du, 'consensus': consensus, 'read_count': dup_read_count})


    dups_consensus_dfm = pd.DataFrame(dups_consensus)


    return dups_consensus_dfm


def consensus_dfm_to_bam(dups_consensus_dfm, template_bam_path, output_bam_path):
    template_bam = pysam.AlignmentFile(template_bam_path, "rb")
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=template_bam)


    for index, row in dups_consensus_dfm.iterrows():
        read_name = f"{row['du']}_consensus"
        
        a = pysam.AlignedSegment()
        a.query_name = read_name
        a.query_sequence = row['consensus']
        a.flag = 4  # unmapped flag
        a.is_unmapped = True
        a.set_tag('du', row['du'])
        zm_tag=row['du'].split('/')[1]
        a.set_tag('zm', int(zm_tag))
        a.set_tag('YC', '240,187,201')
        a.set_tag('ds', int(row['read_count']))


        output_bam.write(a)
    
    output_bam.close()


def dedup_bam(bam_file, output_bam_path, min_read_count=3):
    # First create a dfm from the input table to group reads by pbmarkdup du group tags
    dup_dfm = bam_to_dfm(bam_file)
    # Use MSA to get consensus sequences (using default read numbers here)
    dups_consensus_dfm = consensus_dfm(dup_dfm, min_read_count=min_read_count)
    # Create an unaligned bam file from consensus sequences
    consensus_dfm_to_bam(dups_consensus_dfm, bam_file, output_bam_path)




dedup_bam(bam_file, output_bam, min_read_count= min_read_count)