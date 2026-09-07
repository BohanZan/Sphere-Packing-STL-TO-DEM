[CmdletBinding()]
param()
$ErrorActionPreference = 'Stop'
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio/Installer/vswhere.exe'
if (-not (Test-Path -LiteralPath $vswhere)) { throw 'Visual Studio C++ Build Tools (vswhere.exe) are required.' }
$vsPath = & $vswhere -latest -prerelease -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vsPath) { throw 'No Visual Studio installation with x64 C++ tools was found.' }
$vcvars = Join-Path $vsPath 'VC/Auxiliary/Build/vcvars64.bat'
$intelRoot = Join-Path ${env:ProgramFiles(x86)} 'Intel/oneAPI'
$intelVars = Join-Path $intelRoot 'compiler/latest/env/vars.bat'
if (-not (Test-Path -LiteralPath $intelVars)) { throw "Intel oneAPI compiler is required: $intelVars" }
# Prefer the supported v143 STL when a newer preview Visual Studio is also installed.
$toolsets = Get-ChildItem -LiteralPath (Join-Path $vsPath 'VC/Tools/MSVC') -Directory | Where-Object Name -Like '14.4*' | Sort-Object Name -Descending
$toolsetArg = if ($toolsets) { '-vcvars_ver=' + (($toolsets[0].Name -split '\.')[0..1] -join '.') } else { '' }
$envCommand = 'call "' + $vcvars + '" ' + $toolsetArg + ' >nul && call "' + $intelVars + '" intel64 >nul && set'
$environmentLines = & $env:ComSpec /d /s /c $envCommand
if ($LASTEXITCODE -ne 0) { throw 'Could not initialize x64 Visual Studio and Intel oneAPI.' }
foreach ($line in $environmentLines) {
    if ($line -match '^([^=]+)=(.*)$') { [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process') }
}
# Make sanitizer runtime available to both CTest and interactive runs.
$compilerDir = Join-Path $intelRoot 'compiler/latest'
$clangLib = Join-Path $compilerDir 'lib/clang'
if (Test-Path -LiteralPath $clangLib) {
    $asanDir = Get-ChildItem -LiteralPath $clangLib -Filter 'clang_rt.asan_dynamic-x86_64.dll' -Recurse | Select-Object -First 1 -ExpandProperty DirectoryName
    if ($asanDir) { $env:PATH = "$asanDir;$env:PATH" }
}
Write-Host "oneAPI x64 ready: $((Get-Command icx-cl.exe).Source)"
