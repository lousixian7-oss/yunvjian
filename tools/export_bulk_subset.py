from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
df=pd.read_csv(root/'coredata/bulk/development_expression.csv.gz',index_col=0)
assert df.shape[1]==172 and 'GSM261278' in df.columns
df.drop(columns='GSM261278').to_csv(root/'coredata/bulk/ComBat_WGCNA_171samples.csv.gz',compression='gzip')
