from pathlib import Path
import subprocess, io, shutil, json, os
from PIL import Image, ImageChops
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen.canvas import Canvas

HERE=Path(__file__).resolve().parent
OLD=HERE
FIG=HERE/'composite';FIG.mkdir(exist_ok=True)
POP=Path(os.environ.get('PDFTOPPM') or shutil.which('pdftoppm') or 'C:/Users/32266/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe')
TMP=HERE/'tmp'/'pdfs';TMP.mkdir(parents=True,exist_ok=True)
P=HERE/'panels';P.mkdir(exist_ok=True)

def old_panel(name):
    for ext in ['.pdf','.png']:
        src=OLD/'panels'/(name+ext)
        if not (P/src.name).exists():shutil.copy2(src,P/src.name)
    return name

f5=[old_panel(n) for n in ['Figure5_A_SHAP','Figure5_B_SHAP_importance','Figure5_C_internal_ROC','Figure5_D_expression','Figure5_E_external_ROC','Figure5_F_single_vs_combined_DCA','Figure5_G_nomogram','Figure5_H_calibration']]
e=old_panel('FigureS2_E_probabilities');g=old_panel('FigureS2_G_DCA')

def bounds(path):
    prefix=TMP/(path.stem+'_crop')
    subprocess.run([str(POP),'-f','1','-singlefile','-scale-to','1000','-png',str(path),str(prefix)],check=True,capture_output=True)
    im=Image.open(str(prefix)+'.png').convert('RGB')
    diff=ImageChops.difference(im,Image.new('RGB',im.size,'white')).convert('L').point(lambda x:255 if x>35 else 0)
    b=diff.getbbox()
    if b is None:raise ValueError('Blank source PDF: '+str(path))
    return b,im.size

def assemble(stem,rows,W,H):
    wr=PdfWriter();dest=wr.add_blank_page(W,H)
    labels=io.BytesIO();cv=Canvas(labels,pagesize=(W,H))
    ytop=H;label_idx=0
    for rh,entries in rows:
        x=0;bottom=ytop-rh
        for name,cw in entries:
            path=P/(name+'.pdf');reader=PdfReader(path);page=reader.pages[0]
            assert len(reader.pages)==1
            box,size=bounds(path);pw=float(page.mediabox.width);ph=float(page.mediabox.height)
            l=max(0,box[0]/size[0]*pw-3);r=min(pw,box[2]/size[0]*pw+3)
            b=max(0,(1-box[3]/size[1])*ph-3);t=min(ph,(1-box[1]/size[1])*ph+3)
            page.cropbox.lower_left=(l,b);page.cropbox.upper_right=(r,t)
            scale=min((cw-30)/(r-l),(rh-40)/(t-b))
            px=x+15+(cw-30-(r-l)*scale)/2;py=bottom+9+(rh-40-(t-b)*scale)/2
            tr=Transformation().translate(-l,-b).scale(scale).translate(px,py)
            dest.merge_transformed_page(page,tr)
            cv.setFont('Helvetica-Bold',23);cv.drawString(x+6,ytop-24,chr(65+label_idx));label_idx+=1;x+=cw
        assert abs(x-W)<.01
        ytop=bottom
    cv.save();dest.merge_page(PdfReader(labels).pages[0])
    target=HERE/(stem+'.pdf')
    with target.open('wb') as f:wr.write(f)
    subprocess.run([str(POP),'-singlefile','-r','300','-png',str(target),str(HERE/stem)],check=True,capture_output=True)
    im=Image.open(HERE/(stem+'.png')).convert('RGB');im.save(HERE/(stem+'.tif'),compression='tiff_lzw',dpi=(300,300))
    subprocess.run([str(POP),'-singlefile','-scale-to','1600','-png',str(target),str(TMP/(stem+'_QA'))],check=True,capture_output=True)
    for ext in ['.pdf','.png','.tif']:shutil.copy2(HERE/(stem+ext),FIG/(stem+ext))
    return dict(figure=stem,pages=len(PdfReader(target).pages),bytes=target.stat().st_size,pixels=im.size)

qa=[]
qa.append(assemble('Figure5_three_gene_reconciled_20260905',[(320,[(f5[0],365),(f5[1],310),(f5[2],477)]),(320,[(f5[3],620),(f5[4],532)]),(320,[(f5[5],370),(f5[6],390),(f5[7],392)])],1152,960))
qa.append(assemble('FigureS2_three_gene_reconciled_20260905',[(312,[('FigureS2_A_residual_boxplots',384),('FigureS2_B_residual_reverse_ECDF',384),('FigureS2_C_screening_ROC',384)]),(292,[('FigureS2_D_candidate_rankings',440),(e,712)]),(292,[('FigureS2_F_WGCNA_correlations',440),(g,712)])],1152,896))
(HERE/'figure_export_checks.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa,indent=2))
