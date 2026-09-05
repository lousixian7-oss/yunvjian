from pathlib import Path
import subprocess
from PIL import Image
from pypdf import PdfReader

root = Path(__file__).resolve().parent
pdf = root / 'Figure9L_TGFb_LR_NC_PD.pdf'
reader = PdfReader(pdf)
assert len(reader.pages) == 1
text = reader.pages[0].extract_text()
assert 'TGFB1' in text and 'TGFB3' in text and 'FGF7' not in text
poppler = Path(r'C:\Users\32266\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
subprocess.run([str(poppler), '-singlefile', '-r', '600', '-png', str(pdf), str(pdf.with_suffix(''))], check=True)
im = Image.open(pdf.with_suffix('.png')).convert('RGB')
assert min(im.getextrema()[0]) < 100
im.save(pdf.with_suffix('.tif'), compression='tiff_lzw', dpi=(600, 600))
print(f'Validated: one-page PDF; PNG/TIFF {im.size}, 600 dpi')
