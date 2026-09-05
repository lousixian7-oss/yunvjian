from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas


OUT = Path(r"G:\1Yunvjian\0a26.7.7singlecell\cellchat8.14\cellchat_allcell_original_style_ECM")


def merge_fit(destination, source, box):
    x, y, box_w, box_h = box
    source_w = float(source.mediabox.width)
    source_h = float(source.mediabox.height)
    scale = min(box_w / source_w, box_h / source_h)
    tx = x + (box_w - source_w * scale) / 2
    ty = y + (box_h - source_h * scale) / 2
    destination.merge_transformed_page(
        source,
        Transformation().scale(scale).translate(tx, ty),
        over=True,
    )


# Panel C: retain the original two-circle format without the layout reset
# caused by CellChat 1.6.1.
c_writer = PdfWriter()
c_page = c_writer.add_blank_page(width=14 * 72, height=7 * 72)
pd_page = PdfReader(str(OUT / "Fig8C_PD_enhanced_circle.pdf")).pages[0]
nc_page = PdfReader(str(OUT / "Fig8C_NC_enhanced_circle.pdf")).pages[0]
merge_fit(c_page, pd_page, (0, 0, 7 * 72, 7 * 72))
merge_fit(c_page, nc_page, (7 * 72, 0, 7 * 72, 7 * 72))
c_label_pdf = OUT / "_panel_c_labels.pdf"
c = canvas.Canvas(str(c_label_pdf), pagesize=(14 * 72, 7 * 72))
c.setFont("Helvetica-Bold", 15)
c.drawCentredString(3.5 * 72, 6.70 * 72, "PD-enhanced interaction strength")
c.drawCentredString(10.5 * 72, 6.70 * 72, "NC-enhanced interaction strength")
c.save()
c_page.merge_page(PdfReader(str(c_label_pdf)).pages[0], over=True)
with (OUT / "Fig8C_allcell_PD_NC_difference_circle.pdf").open("wb") as handle:
    c_writer.write(handle)
c_label_pdf.unlink(missing_ok=True)


# Full Figure 8 layout.
page_w = 18 * 72
page_h = 23 * 72
layout = {
    "A": (0.55 * 72, 17.25 * 72, 8.0 * 72, 5.25 * 72),
    "B": (8.75 * 72, 17.30 * 72, 8.70 * 72, 4.35 * 72),
    "C": (0.55 * 72, 9.15 * 72, 16.90 * 72, 7.80 * 72),
    "D": (0.55 * 72, 0.45 * 72, 8.35 * 72, 8.35 * 72),
    "E": (9.05 * 72, 1.35 * 72, 8.40 * 72, 5.85 * 72),
}
panel_files = {
    "A": "Fig8A_allcell_communication_number_strength.pdf",
    "B": "Fig8B_allcell_incoming_outgoing_NC_PD.pdf",
    "C": "Fig8C_allcell_PD_NC_difference_circle.pdf",
    "D": "Fig8D_allcell_selected_pathway_networks.pdf",
    "E": "Fig8E_ECM_Fib_to_Activated_Fib_LR_comparison.pdf",
}

writer = PdfWriter()
page = writer.add_blank_page(width=page_w, height=page_h)
for label, file_name in panel_files.items():
    source = PdfReader(str(OUT / file_name)).pages[0]
    merge_fit(page, source, layout[label])

label_pdf = OUT / "_panel_labels.pdf"
c = canvas.Canvas(str(label_pdf), pagesize=(page_w, page_h))
c.setFont("Helvetica-Bold", 20)
for label, (x, y, _w, h) in layout.items():
    c.drawString(x - 12, y + h - 4, label)
c.save()
page.merge_page(PdfReader(str(label_pdf)).pages[0], over=True)

final_pdf = OUT / "Figure8_CellChat_allcell_ECM_Fib_original_style.pdf"
with final_pdf.open("wb") as handle:
    writer.write(handle)
label_pdf.unlink(missing_ok=True)
print(final_pdf)
