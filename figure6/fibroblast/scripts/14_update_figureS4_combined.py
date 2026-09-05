from pathlib import Path
from shutil import copy2
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

OUT = Path(r"G:\1Yunvjian\FigureS4_fibroblast_dataset_bias")

# Panel G is now the mandatory stress-adjusted sensitivity analysis. Preserve
# one canonical name so the combined figure and standalone delivery agree.
copy2(OUT / "S4G2_stress_adjusted_activated_model.png", OUT / "S4G.png")
copy2(OUT / "S4G2_stress_adjusted_activated_model.pdf", OUT / "S4G.pdf")
copy2(OUT / "S4G2_stress_adjusted_activated_model.pdf",
      OUT / "S4G_stress_adjusted_sensitivity.pdf")
copy2(OUT / "S4H_cluster_robustness.png", OUT / "S4H.png")
copy2(OUT / "S4H_cluster_robustness.pdf", OUT / "S4H.pdf")

panel_paths = [
    OUT / "S4A.png",
    OUT / "S4B.png",
    OUT / "S4C.png",
    OUT / "S4D.png",
    OUT / "S4E.png",
    OUT / "S4F.png",
    OUT / "S4G.png",
    OUT / "S4H.png",
]
labels = list("ABCDEFGH")

canvas_w, canvas_h = 8400, 4080
margin_x, margin_y = 85, 90
gap_x, gap_y = 55, 80
cell_w = (canvas_w - 2 * margin_x - 3 * gap_x) // 4
cell_h = (canvas_h - 2 * margin_y - gap_y) // 2

combined = Image.new("RGB", (canvas_w, canvas_h), "white")
draw = ImageDraw.Draw(combined)
font_path = Path(r"C:\Windows\Fonts\arialbd.ttf")
label_font = ImageFont.truetype(str(font_path), 82) if font_path.exists() else None

for idx, (path, label) in enumerate(zip(panel_paths, labels)):
    row, col = divmod(idx, 4)
    x0 = margin_x + col * (cell_w + gap_x)
    y0 = margin_y + row * (cell_h + gap_y)
    panel = Image.open(path).convert("RGB")
    # Trim excessive standalone-figure whitespace before fitting the panel.
    bg = Image.new("RGB", panel.size, "white")
    bbox = ImageChops.difference(panel, bg).getbbox()
    if bbox:
        pad = 35
        bbox = (max(0, bbox[0] - pad), max(0, bbox[1] - pad),
                min(panel.width, bbox[2] + pad), min(panel.height, bbox[3] + pad))
        panel = panel.crop(bbox)
    # Reserve room for the bold panel label and fit the plot without distortion.
    plot_box = (cell_w - 35, cell_h - 50)
    panel = ImageOps.contain(panel, plot_box, Image.Resampling.LANCZOS)
    px = x0 + (cell_w - panel.width) // 2 + 15
    py = y0 + (cell_h - panel.height) // 2 + 15
    combined.paste(panel, (px, py))
    draw.text((x0, y0 - 22), label, fill="black", font=label_font)

png_path = OUT / "FigureS4_combined.png"
combined.save(png_path, dpi=(300, 300), quality=95)

pdf_path = OUT / "FigureS4_combined.pdf"
page_w = 28 * 72
page_h = page_w * canvas_h / canvas_w
c = canvas.Canvas(str(pdf_path), pagesize=(page_w, page_h))
c.drawImage(ImageReader(combined), 0, 0, width=page_w, height=page_h,
            preserveAspectRatio=True, mask="auto")
c.showPage()
c.save()

print(png_path)
print(pdf_path)
