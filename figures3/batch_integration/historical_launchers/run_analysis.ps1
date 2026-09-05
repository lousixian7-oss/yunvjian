$ErrorActionPreference = 'Stop'

$analysisRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$rscript = 'C:\Program Files\R\R-4.5.2\bin\x64\Rscript.exe'
$localLibrary = Join-Path $analysisRoot 'R_library'
$inputRds = Join-Path (Split-Path -Parent $analysisRoot) 'yjsl_clean_allGSE.rds'
$installScript = Join-Path $analysisRoot 'scripts\00_install_packages.R'
$analysisScript = Join-Path $analysisRoot 'scripts\01_batch_integration_revision.R'

if (-not (Test-Path -LiteralPath $rscript)) {
    throw "Rscript not found: $rscript"
}
if (-not (Test-Path -LiteralPath $inputRds)) {
    throw "Input RDS not found: $inputRds"
}

New-Item -ItemType Directory -Force -Path $localLibrary | Out-Null

& $rscript $installScript $localLibrary
if ($LASTEXITCODE -ne 0) {
    throw "R package installation failed with exit code $LASTEXITCODE"
}

& $rscript $analysisScript $inputRds $analysisRoot $localLibrary
if ($LASTEXITCODE -ne 0) {
    throw "Batch integration analysis failed with exit code $LASTEXITCODE"
}

Write-Output "Analysis completed: $analysisRoot"
