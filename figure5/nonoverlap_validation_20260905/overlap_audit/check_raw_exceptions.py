from pathlib import Path
import urllib.request,gzip,hashlib,tarfile,json,concurrent.futures
out=Path(__file__).resolve().parent
ids=['GSM261185','GSM261186','GSM261293']
def fetch(gsm):
 p=out/(gsm+'.CEL.gz')
 if not p.exists():
  url='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM261nnn/'+gsm+'/suppl/'+gsm+'.CEL.gz'
  with urllib.request.urlopen(url,timeout=120) as r:p.write_bytes(r.read())
 data=gzip.decompress(p.read_bytes())
 return gsm,hashlib.sha256(data).hexdigest(),data[:4096].hex()
source=list(concurrent.futures.ThreadPoolExecutor(3).map(fetch,ids))
print('downloaded source CEL files',flush=True)
# Hash every external CEL without altering or extracting the archive.
matches=[]
with tarfile.open('G:/1Yunvjian/periodontitis基因/GEO/外部数据/GSE16134_RAW.tar') as t:
 for member in t:
  if not member.isfile() or not member.name.lower().endswith('.gz'):continue
  raw=gzip.decompress(t.extractfile(member).read()); digest=hashlib.sha256(raw).hexdigest()
  for gsm,sha,header in source:
   if digest==sha: matches.append(dict(source=gsm,external_file=member.name,sha256=sha))
result={'source_hashes':[{'sample':s,'sha256':h} for s,h,_ in source],'exact_raw_matches':matches}
(out/'raw_CEL_exception_check.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
