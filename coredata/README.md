# Core processed data

`bulk`: complete 172-sample ComBat matrix and the 171-sample WGCNA subset, with external expression and sample annotation.

`metadata`: original global and final fibroblast cell metadata; source extraction comes from the manuscript table preparation task. Embeddings are exported separately if available.

`spatial`: each GEO section has its own RCTD weights and summary files. The NC annotation object has 4,992 spots, while its saved RCTD subset has 1,088. Do not merge these without barcode matching.

`model`, `fibroblast`, `hdwgcna`, `cellchat`: saved processed evidence associated with the archived figure workflows. Some sensitivity outputs intentionally preserve alternative results.
