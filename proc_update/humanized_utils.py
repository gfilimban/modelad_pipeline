import pandas as pd
import pyranges as pr
import pyfaidx
import re
import textwrap


def get_gene_t_fastq(wc,
                     fa_file,
                     ofile):
    """
    Get the fastq for the reference transcript associated with
    "gene"
    """
    import pdb; pdb.set_trace()
    if 'mouse_gene' in list(wc.keys()):
        species = 'mouse'
        gene = wc['mouse_gene']
    elif 'human_gene' in list(wc.keys()):
        species = 'human'
        gene = wc['human_gene']

    if gene == 'dummy':
        pathlib.Path(ofile).touch()
        return

    indexname = f'{ofile}.fai'
    fa = pyfaidx.Fasta(fa_file, indexname=indexname)
    t_keys = list(fa.keys())
    t_df = pd.DataFrame()
    gene_str = f'|{gene}|'
    keys = [t for t in t_keys if gene_str in t]

    if len(keys) == 0:
        raise ValueError('Try again; probably frozen .fai file')

    with open(ofile, 'w') as o:
        for t in keys:
            read_name = f'@{t}'
            read = fa[t][:].seq
            phred = ''.join(['5' for i in range(len(read))])
            o.write(read_name+'\n')
            o.write(read+'\n')
            o.write('+\n')
            o.write(phred+'\n')

def write_chr(seq, ofile, chr_name, line_lim=None):
    with open(ofile, 'w') as ofile:
        ofile.write(f'>{chr_name}\n')
        if line_lim:
            wrapped_seq = textwrap.wrap(seq, width=line_lim)
            for seq in wrapped_seq:
                ofile.write(seq+'\n')
        else:
            ofile.write(str(seq))


def get_gene_gtf_entry(gtf_file, gene):
    """
    Get the GTF entries from the gene
    """

    gtf_df = pr.read_gtf(gtf_file, as_df=True)

    gtf_df = gtf_df.loc[(gtf_df.gene_name==gene)]

    return gtf_df


def get_gene_seq(fa_file,
                gtf_file,
                gene,
                whole_chr=False,
                ofile=None,
                chr_name=None,
                slack=0):

    gtf_df = pr.read_gtf(gtf_file, as_df=True)

    gtf_df = gtf_df.loc[(gtf_df.gene_name==gene)&(gtf_df.Feature=='gene')]
    assert len(gtf_df.index) == 1

    # get start and end of gene
    start = min(gtf_df['Start'].values[0], gtf_df['End'].values[0])-slack
    end = max(gtf_df['Start'].values[0], gtf_df['End'].values[0])+slack
    ch = gtf_df['Chromosome'].values[0]
    print(ch)
    print(start)
    print(end)
    strand = gtf_df['Strand'].values[0]

    fa = pyfaidx.Fasta(fa_file)

    # just get the sequence of the gene
    if not whole_chr:
        if strand == '+':
            gene_seq = fa[ch][start:end]
        else:
            gene_seq = fa[ch][start:end].complement
    # get the sequence of the whole chromosome
    else:
        gene_seq = fa[ch][:]
    print(f'Length of pseudochrom: {len(gene_seq.seq)}')

    if ofile:
        if not chr_name:
            chr_name = gene

        write_chr(gene_seq.seq, ofile, chr_name)

    return gene_seq.seq

def replace_seq(mod_fa,
                seq):
    """
    Replace sequence in seq with sequence in mod_fa. Mod_fa must contain some overlap
    with the sequence in seq to find the insertion sites.

    Parameters:
        mod_fa (str): Path to fasta file with sequence to insert
        seq (str): Sequence to add mod_fa into
    """

    # get anchors at 5' and 3' ends of sequence
    fa = pyfaidx.Fasta(modified_fa)
    mod_seq = fa[list(fa.keys())[0]][0:-1].seq

    anchor_5 = mod_seq[:20]
    anchor_3 = mod_seq[-20:]

    start_inds = [(m.start(), m.end()) for m in re.finditer(anchor_5, seq, re.IGNORECASE)]
    assert len(start_inds) == 1

    end_inds = [(m.start(), m.end()) for m in re.finditer(anchor_3, seq, re.IGNORECASE)]
    assert len(end_inds) == 1

    assert seq[start_inds[0][0]:start_inds[0][1]].lower() == anchor_5.lower()
    assert seq[end_inds[0][0]:end_inds[0][1]].lower() == anchor_3.lower()

    # replace_seq = seq[start_inds[0][1]:end_inds[0][0]] # this will replace excluding anchors
    replace_seq = seq[start_inds[0][0]:end_inds[0][1]] # this will replace inluding anchors
    with open('ref/replace_seq.fa', 'w') as ofile:
        ofile.write(replace_seq)

    assert replace_seq[:20].lower() == anchor_5.lower()
    assert replace_seq[-20:].lower() == anchor_3.lower()

    # print(replace_seq[0])
    # print(replace_seq[-1])

    inds = [(m.start(), m.end()) for m in re.finditer(replace_seq, seq)]
    assert len(inds) == 1

    new_len = len(seq)-len(replace_seq)+len(mod_seq)
    new_seq = seq.replace(replace_seq, mod_seq)
    assert len(new_seq) == new_len
    return new_seq

def refmt_mapped_transcript_gtf(wc, locus_type, ifile, ofile):
    """
    Format a TALON gtf output from annotated transcripts mapped
    onto pseudochromosomes in order to make them appear as known.
    """
    if 'mouse_gene' in list(wc.keys()):
        species = 'mouse'
        gene = wc['mouse_gene']
    elif 'human_gene' in list(wc.keys()):
        species = 'human'
        gene = wc['human_gene']

    new_gene = f'{species}_{gene}'

    if gene == 'dummy':
        pathlib.Path(ofile).touch()
        return

    # for whole chrom, extract all GTF entries from the corresponding gene
    elif locus_type == 'human':
        if species == 'human':
            gtf_file = human_annot
        elif species == 'mouse':
            gtf_file = annot
        df = pr.read_gtf(gtf_file).df
        df = df.loc[df.gene_name == gene]

    # for other cases, use the mapped transcript TALON output
    else:
        df = pr.read_gtf(ifile, rename_attr=True).df

        keep_cols = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
           'Frame', 'gene_id', 'gene_name', 'gene_status', 'source_attr',
           'transcript_id', 'transcript_status', 'transcript_name', 'exon_number', 'exon_id',
           'exon_status']
        df = df[keep_cols]

        df.loc[df.Feature.isin(['transcript','exon']), 'transcript_status'] = 'KNOWN'
        df.loc[df.Feature=='exon', 'exon_status'] = 'KNOWN'
        df.gene_status = 'KNOWN'

    # assign these models new gene names and write to output
    df.gene_id = new_gene
    df.gene_name = new_gene

    pr.PyRanges(df).to_gtf(ofile)
