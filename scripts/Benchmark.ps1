[CmdletBinding()]
param([switch]$SkipMatlab,[int]$Repetitions=3)
$ErrorActionPreference = 'Stop'
if ($Repetitions -lt 1) { throw 'Repetitions must be positive.' }
$projectPath = Split-Path -Parent $PSScriptRoot
Push-Location -LiteralPath $projectPath
try {
    if (-not $SkipMatlab) {
        if (-not (Test-Path -LiteralPath tmp)) { New-Item -ItemType Directory tmp | Out-Null }
        & matlab -batch "addpath(fullfile(pwd,'benchmarks')); runMatlabBaseline('box',$Repetitions); runMatlabBaseline('refill',$Repetitions); runMatlabBaseline('buddha',$Repetitions); runMatlabBaseline('buddha_gravity',$Repetitions);" -logfile tmp/cpp_matlab_benchmark.log
        if ($LASTEXITCODE -ne 0) { throw 'MATLAB baseline failed.' }
    }
    & python benchmarks/compare.py --scenario box refill buddha buddha_gravity --repetitions $Repetitions
    if ($LASTEXITCODE -ne 0) { throw 'C++ replay or timing comparison failed.' }
} finally { Pop-Location }
