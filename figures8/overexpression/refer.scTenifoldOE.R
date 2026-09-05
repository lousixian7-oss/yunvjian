#################### 虚拟过表达分析函数 ####################
# 参数说明：
#   countMatrix : 基因表达矩阵（行=基因，列=细胞）
#   gOE         : 目标过表达基因名称
#   oe_factor   : 过表达倍数（默认10）
#   ...         : 其他参数
# 返回值：列表，包含 $tensorNetworks, $manifoldAlignment, $diffRegulation

scTenifoldOE <- function(countMatrix,
                         gOE,
                         oe_factor = 10,
                         qc = TRUE,
                         qc_mtThreshold = 0.1,
                         qc_minLSize = 500,
                         qc_minCells = 25,
                         nc_lambda = 0,
                         nc_nNet = 10,
                         nc_nCells = 500,
                         nc_nComp = 3,
                         nc_scaleScores = TRUE,
                         nc_symmetric = FALSE,
                         nc_q = 0.9,
                         td_K = 3,
                         td_maxIter = 1000,
                         td_maxError = 1e-5,
                         td_nDecimal = 3,
                         ma_nDim = 2,
                         nCores = parallel::detectCores()) {
  
  cli::cli_h1("scTenifoldOE")
  cli::cli_alert_info("Simulating {gOE} gene overexpression (factor = {oe_factor})")
  
  # 检查目标基因是否存在
  if (!gOE %in% rownames(countMatrix)) {
    cli::cli_alert_danger("{gOE} is not present in the count matrix")
    cli::cli_abort("{gOE} not found")
  }
  
  # 可选的质量控制
  if (isTRUE(qc)) {
	countMatrix=as.matrix(countMatrix)
    countMatrix <- scTenifoldNet::scQC(countMatrix, 
                                       maxMTratio = qc_mtThreshold,
								 minPCT = 0,
                                       minLibSize = qc_minLSize)
    #countMatrix <- countMatrix[rowSums(countMatrix != 0) >= qc_minCells, ]
    #cli::cli_alert_success("QC applied: retained {nrow(countMatrix)} genes, {ncol(countMatrix)} cells")
    if (!gOE %in% rownames(countMatrix)) {
      cli::cli_abort("{gOE} lost after QC")
    }
  }
  
  # Step 1: 构建多个子网络
  cli::cli_alert_info("Building {nc_nNet} networks with {nc_nCells} cells each...")
  WT_networks <- scTenifoldNet::makeNetworks(
    X = countMatrix,
    q = nc_q,
    nNet = nc_nNet,
    nCells = nc_nCells,
    scaleScores = nc_scaleScores,
    symmetric = nc_symmetric,
    nComp = nc_nComp,
    nCores = nCores
  )
  
  # Step 2: 张量分解得到共识网络
  cli::cli_alert_info("Tensor decomposition (K = {td_K})...")
  WT_tensor <- scTenifoldNet::tensorDecomposition(
    xList = WT_networks,
    K = td_K,
    maxError = td_maxError,
    maxIter = td_maxIter,
    nDecimal = td_nDecimal
  )
  WT <- WT_tensor$X
  
  # Step 3: 官方流程中的后处理（strictDirection, diag, transpose）
  WT <- scTenifoldKnk:::strictDirection(WT, lambda = nc_lambda)
  WT <- as.matrix(WT)
  diag(WT) <- 0
  WT <- t(WT)
  cli::cli_alert_success("Prepared WT adjacency matrix")
  
  # Step 4: 模拟过表达
  cli::cli_alert_info("Simulating overexpression of {gOE}")
  OE <- WT
  OE[gOE, ] <- OE[gOE, ] * oe_factor

  # Step 5: 流形对齐
  cli::cli_alert_info("Manifold alignment (d = {ma_nDim})...")
  MA <- scTenifoldNet::manifoldAlignment(WT, OE, d = ma_nDim, nCores = nCores)
  
  # Step 6: 识别差异调控基因（使用官方函数）
  DR <- scTenifoldKnk:::dRegulation(MA, gKO = gOE)
  
  # 组装输出
  outputList <- list()
  outputList$tensorNetworks$WT <- Matrix::Matrix(WT)
  outputList$tensorNetworks$OE <- Matrix::Matrix(OE)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  
  cli::cli_alert_success("Finished scTenifoldOE for {gOE}")
  return(outputList)
}



#################### 方法2: 表达矩阵出发 ####################
scTenifoldOE2 <- function(countMatrix, gOE, nCores = 5) {
  
  # 1. 过滤：保留在 ≥5% 细胞中表达的基因
  min.cells <- ncol(countMatrix) * 0.05
  countMatrix <- countMatrix[rowSums(countMatrix > 0) >= min.cells, ]
  
  # 2. 随机下采样：最多 10,000 个细胞
  if(ncol(countMatrix) > 10000){
    set.seed(123)
    countMatrix <- countMatrix[, sample(ncol(countMatrix), 10000)]
  }
  
  # 3. 检查目标基因是否存在
  if(!gOE %in% rownames(countMatrix)){
    stop(paste(gOE, "not found or filtered out (expressed in <5% of cells)"))
  }
  
  # 4. 模拟过表达：additive scheme → 所有细胞 +1（文章严格写法）
  oeMatrix <- countMatrix
  oeMatrix[gOE, ] <- oeMatrix[gOE, ] + 1  
  
  # 5. 构建野生型与过表达网络（参数完全匹配文章）
  res <- scTenifoldNet::scTenifoldNet(
    X = countMatrix,        # WT
    Y = oeMatrix,           # OE
    qc = FALSE,             # 已手动过滤
    nc_nNet = 10,           # 文章：nNet = 10
    nc_nCells = 200,        # 文章：nCells = 200
    nc_nComp = 3,           # 文章：nComp = 3
    td_K = 3,               # 文章：K = 3
    nCores = nCores
  )
  
  message(paste0("Virtual overexpression completed for: ", gOE))
  return(res)
}


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


