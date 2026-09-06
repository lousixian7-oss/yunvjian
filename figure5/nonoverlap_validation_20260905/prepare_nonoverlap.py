from pathlib import Path
import csv,gzip,tarfile,hashlib,json
root=Path(__file__).resolve().parent
out=root/'nonoverlap_validation';out.mkdir(exist_ok=True)
cel=out/'CEL';cel.mkdir(exist_ok=True)
keep=list(csv.DictReader((root/'overlap_audit/GSE16134_nonoverlapping_patient_samples.csv').open(encoding='utf-8-sig')))
ids={r['sample_id'] for r in keep};assert len(ids)==66
manifest=[]
with tarfile.open('G:/1Yunvjian/periodontitis基因/GEO/外部数据/GSE16134_RAW.tar') as tar:
 for m in tar:
  if Path(m.name).name.split('.')[0] not in ids:continue
  target=cel/Path(m.name).name
  if not target.exists():target.write_bytes(tar.extractfile(m).read())
  manifest.append({'file':target.name,'sha256_compressed':hashlib.sha256(target.read_bytes()).hexdigest()})
assert len(manifest)==66
archive=Path('G:/1Yunvjian/1A玉女煎/code/coredata/bulk')
for src,dst in [('development_expression.csv.gz','development_three_genes.csv'),('GSE16134_expression.csv.gz','full_cohort_three_genes.csv')]:
 with gzip.open(archive/src,'rt',encoding='utf-8-sig',newline='') as f, (out/dst).open('w',encoding='utf-8',newline='') as w:
  rows=csv.reader(f);writer=csv.writer(w);writer.writerow(next(rows));n=0
  for row in rows:
   if row[0] in ['AGT','CXCR4','FOS']:writer.writerow(row);n+=1
  assert n==3
(out/'CEL_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Prepared 66 original CEL files and three-gene comparison inputs.',flush=True)
