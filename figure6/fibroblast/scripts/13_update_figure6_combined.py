from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

ROOT = Path(r"G:\1Yunvjian\0a26.7.7singlecell\8.10\fibroblast_revision_final")
FIG = ROOT / "results" / "figures" / "figure6_revised_GK"

combined_path = FIG / "Figure6_GK_revised_combined.png"
panel_i_path = FIG / "Figure6I_sample_level_activated_program.png"

base = Image.open(combined_path).convert("RGB")
panel_i = Image.open(panel_i_path).convert("RGB")

# Replace the former pooled-composition panel only; panels G, H, J and K are
# preserved byte-for-byte at their existing raster positions.
draw = ImageDraw.Draw(base)
replace_box = (2240, 0, base.width, 1545)
draw.rectangle(replace_box, fill="white")

target_w = base.width - 2380
target_h = int(panel_i.height * target_w / panel_i.width)
if target_h > 1370:
    target_h = 1370
    target_w = int(panel_i.width * target_h / panel_i.height)
panel_i = panel_i.resize((target_w, target_h), Image.Resampling.LANCZOS)
paste_x = 2320 + max(0, (base.width - 2320 - target_w) // 2)
paste_y = 75
base.paste(panel_i, (paste_x, paste_y))

font_path = Path(r"C:\Windows\Fonts\arialbd.ttf")
font = ImageFont.truetype(str(font_path), 74) if font_path.exists() else None
draw = ImageDraw.Draw(base)
draw.text((2260, 25), "I", fill="black", font=font)

base.save(combined_path, dpi=(300, 300), quality=95)

pdf_path = FIG / "Figure6_GK_revised_combined.pdf"
page_w = 16 * 72
page_h = page_w * base.height / base.width
c = canvas.Canvas(str(pdf_path), pagesize=(page_w, page_h))
c.drawImage(ImageReader(base), 0, 0, width=page_w, height=page_h,
            preserveAspectRatio=True, mask="auto")
c.showPage()
c.save()

print(combined_path)
print(pdf_path)
