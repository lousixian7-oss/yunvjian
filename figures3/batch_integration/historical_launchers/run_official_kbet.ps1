$ErrorActionPreference = 'Stop'

$analysisRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$rscript = 'C:\Program Files\R\R-4.5.2\bin\x64\Rscript.exe'
$localLibrary = Join-Path $analysisRoot 'R_library'
$inputRds = Join-Path $analysisRoot 'results\objects\yjsl_batch_integration_revision.rds'
$analysisScript = Join-Path $analysisRoot 'scripts\03_official_kbet.R'

& $rscript $analysisScript $inputRds $analysisRoot $localLibrary
if ($LASTEXITCODE -ne 0) {
    throw "Official kBET analysis failed with exit code $LASTEXITCODE"
}

Write-Output "Official kBET analysis completed: $analysisRoot"
