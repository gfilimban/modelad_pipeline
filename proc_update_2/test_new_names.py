import pandas as pd
import re

df = pd.DataFrame()
df['fname'] = ['/share/crsp/lab/model-ad/share/LRRNA_nanopore/lr-ad009_mTrem2-T96k_8wk/ad009-1_mux_p2_1/20240207_1438_P2S-00617-B_PAU70051_dc91e88a/trim_fastq/ad009-1_22752_mux-blk_1_t1.fastq.gz',
               '/share/crsp/lab/seyedam/share/ad_nanopore_tmp/AD002/ad002_trimfastq/ad002_10721_lig-blk_2_t1.fastq.gz']
df['basename'] = df.fname.str.rsplit('/', n=1, expand=True)[1]
df['path'] = df.fname.str.rsplit('/', n=1, expand=True)[0]

############ Dataset + flowcell df

# get flowcell
exp = '.*\/[\w-]+_(\d+)(?:_t\d+)?\.fastq(?:.gz)?'
df['flowcell'] = df.fname.str.extract(exp)

# check to make sure the same file stem isn't there more than once
# (can happen if different flow cells needed different amounts of chopping)
# df['file_stem'] = df.basename.str.rsplit('_', n=1, expand=True)[0]
exp = '.*\/([\w-]+_\d+)(?:_t\d+)?\.fastq(?:.gz)?'
df['file_stem'] = df.fname.str.extract(exp)
df['chop_num'] = df.basename.str.rsplit('.fastq', expand=True)[0].str.rsplit('_t', expand=True)[1].astype(float)
if df.file_stem.duplicated().any():
dupe_stems = df.loc[df.file_stem.duplicated(keep=False), 'basename'].tolist()
if not auto_dedupe:
  raise ValueError(f'Files {dupe_stems} seem to be duplicated. Check config file.')
else:
  print(f'Files {dupe_stems} seem to be duplicated. Automatically removing lower chop numbers')
  df = df.sort_values(by='chop_num', ascending=False)
  df = df.drop_duplicates(subset='file_stem', keep='first')

# extract the sample name
temp = df.basename.str.split('_', expand=True)[[0,1]]#.str.join('_')
df['sample_temp'] = temp[0]+'_'+temp[1]

# extract the mouse id
df['mouse_id'] = df['sample_temp'].str.split('_', expand=True)[1]

# extract the "study" name
exp = '^(ad[0-9]+)'
df['study'] = df.basename.str.extract(exp)

# extract the "library" name
exp = '^ad[0-9]+-?([\d])?'
df['library'] = df.basename.str.extract(exp)
df.loc[df.library.isnull(), 'library'] = 1
df.library = df.library.astype(int)
