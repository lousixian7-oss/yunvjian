$ErrorActionPreference = 'Stop'

$analysisRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$rscript = 'C:\Program Files\R\R-4.5.2\bin\x64\Rscript.exe'
$localLibrary = Join-Path $analysisRoot 'R_library'
$inputRds = Join-Path $analysisRoot 'results\objects\yjsl_batch_integration_revision.rds'
$analysisScript = Join-Path $analysisRoot 'scripts\02_cluster_stability_and_kbet.R'

if (-not (Test-Path -LiteralPath $rscript)) {
    throw "Rscript not found: $rscript"
}
if (-not (Test-Path -LiteralPath $inputRds)) {
    throw "Input RDS not found: $inputRds"
}

& $rscript $analysisScript $inputRds $analysisRoot $localLibrary
if ($LASTEXITCODE -ne 0) {
    throw "Stability analysis failed with exit code $LASTEXITCODE"
}

Write-Output "Stability analysis completed: $analysisRoot"
