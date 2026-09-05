from pathlib import Path

import pypdfium2 as pdfium
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
SPECS = [
    ("COLLAGEN", "0.603", "A"),
    ("FN1", "0.160", "B"),
    ("LAMININ", "0.134", "C"),
    ("PERIOSTIN", "0.0513", "D"),
]
POSITIONS = [(20, 0), (1260, 0), (20, 1240), (1260, 1240)]

canvas = Image.new("RGB", (2480, 2480), "white")
draw = ImageDraw.Draw(canvas)
font_bold = ImageFont.truetype(r"C:\Windows\Fonts\arialbd.ttf", 34)
font_regular = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", 22)

for (pathway, delta, label), (x, y) in zip(SPECS, POSITIONS):
    document = pdfium.PdfDocument(str(ROOT / f"Figs6_{pathway}_panel.pdf"))
    panel = document[0].render(scale=150 / 72).to_pil().convert("RGB")
    # CellChat 1.6.1 uses a fixed row height. Compress the complete, uncropped
    # pathway page to reproduce the compact matrix style of the original Fig. S6F.
    panel = panel.resize((1200, 1125), Image.Resampling.LANCZOS)

    draw.text((x + 12, y + 12), f"{label}  {pathway}", font=font_bold, fill="black")
    draw.text(
        (x + 12, y + 54),
        f"Communication probability: PD - NC = {delta}",
        font=font_regular,
        fill="black",
    )
    canvas.paste(panel, (x, y + 100))

canvas.save(
    ROOT / "Figs6.tif",
    format="TIFF",
    compression="tiff_lzw",
    dpi=(300, 300),
)
canvas.save(ROOT / "Figs6_pathway_violins.pdf", format="PDF", resolution=300)
canvas.resize((1600, 1600), Image.Resampling.LANCZOS).save(ROOT / "Figs6_preview.png")
