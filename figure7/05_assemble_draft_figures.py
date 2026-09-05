from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

root = Path(r"G:\1Yunvjian\0a26.7.7singlecell\8.10\fibroblast_revision_final\hdWGCNA_rerun_20260824")
fig = root / "figures"
font = ImageFont.truetype(r"C:\Windows\Fonts\arialbd.ttf", 120)


def place(canvas, image_path, box, label):
    x, y, w, h = box
    image = Image.open(image_path).convert("RGB")
    image.thumbnail((w, h), Image.Resampling.LANCZOS)
    px = x + (w - image.width) // 2
    py = y + (h - image.height) // 2
    canvas.paste(image, (px, py))
    ImageDraw.Draw(canvas).text((x + 10, y + 5), label, font=font, fill="black")


figure7 = Image.new("RGB", (6000, 6500), "white")
place(figure7, fig / "Figure7A_soft_power.png", (100, 80, 2850, 1700), "A")
place(figure7, fig / "Figure7B_dendrogram.png", (3050, 80, 2850, 1700), "B")
place(figure7, fig / "Figure7C_module_hub_genes.png", (100, 1880, 3350, 2350), "C")
place(figure7, fig / "Figure7D_sample_level_module_state_heatmap.png", (3550, 1880, 2350, 2350), "D")
place(figure7, fig / "Figure7E_AGT_selected_module_sample_correlation.png", (100, 4400, 1900, 1900), "E")
place(figure7, fig / "Figure7F_independent_module_score_by_state.png", (2050, 4400, 1900, 1900), "F")
place(figure7, fig / "Figure7G_selected_module_PD_NC_sample_level.png", (4000, 4400, 1900, 1900), "G")
figure7.save(root / "Figure7_reclustered_hdwgcna_draft.tif", compression="tiff_lzw", dpi=(400, 400))
preview7 = figure7.copy()
preview7.thumbnail((2200, 2400), Image.Resampling.LANCZOS)
preview7.save(root / "Figure7_reclustered_hdwgcna_draft_preview.png")

figs5 = Image.new("RGB", (6000, 6200), "white")
place(figs5, fig / "FigureS5A_module_feature_plots.png", (100, 80, 3650, 2300), "A")
place(figs5, fig / "FigureS5B_module_correlogram.png", (3850, 80, 2050, 2300), "B")
place(figs5, fig / "FigureS5C_hub_gene_network.png", (100, 2500, 2550, 1750), "C")
place(figs5, fig / "FigureS5D_modules_1_2_6_by_revised_fibroblast_state.png", (2750, 2500, 3150, 1750), "D")
place(figs5, fig / "FigureS5E_selected_module_GO_BP.png", (100, 4450, 2850, 1600), "E")
place(figs5, fig / "FigureS5F_selected_module_KEGG.png", (3050, 4450, 2850, 1600), "F")
figs5.save(root / "FigureS5_reclustered_hdwgcna_draft.tif", compression="tiff_lzw", dpi=(400, 400))
previews5 = figs5.copy()
previews5.thumbnail((2400, 2400), Image.Resampling.LANCZOS)
previews5.save(root / "FigureS5_reclustered_hdwgcna_draft_preview.png")
