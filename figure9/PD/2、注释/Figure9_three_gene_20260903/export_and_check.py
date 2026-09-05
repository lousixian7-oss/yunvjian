from pathlib import Path
import subprocess, json
from PIL import Image, ImageChops
from pypdf import PdfReader
p=Path(__file__).resolve().parent
pop=Path('C:/Users/32266/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe')
qa=p/'tmp'/'pdfs';qa.mkdir(parents=True,exist_ok=True)
records=[]
for f in sorted((p/'panels').glob('*.pdf')):
 r=PdfReader(f);assert len(r.pages)==1
 text=r.pages[0].extract_text()
 for g in ['AGT','CXCR4','FOS','PI16','COL1A1','FN1','POSTN']:assert g in text,(f.name,g)
 subprocess.run([str(pop),'-singlefile','-r','300','-png',str(f),str(f.with_suffix(''))],check=True,capture_output=True)
 im=Image.open(f.with_suffix('.png')).convert('RGB')
 assert ImageChops.difference(im,Image.new('RGB',im.size,'white')).getbbox()
 im.save(f.with_suffix('.tif'),compression='tiff_lzw',dpi=(300,300))
 subprocess.run([str(pop),'-singlefile','-scale-to','2200','-png',str(f),str(qa/f.stem)],check=True,capture_output=True)
 records.append({'file':f.name,'pages':1,'pixels':im.size,'nonblank':True})
(p/'export_QA.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print(json.dumps(records,indent=2))
