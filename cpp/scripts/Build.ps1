[CmdletBinding()]
param([switch]$Sanitize, [switch]$NoTest)
$ErrorActionPreference = 'Stop'
. "$PSScriptRoot/Enter-OneAPI.ps1"
$projectPath = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$buildPath = Join-Path $projectPath $(if ($Sanitize) { 'build/asan' } else { 'build/release' })
$asan = if ($Sanitize) { 'ON' } else { 'OFF' }
$compilerPath = Join-Path ${env:ProgramFiles(x86)} 'Intel/oneAPI/compiler/latest/bin/icx-cl.exe'
& cmake --fresh -S $projectPath -B $buildPath -G Ninja "-DCMAKE_CXX_COMPILER=$compilerPath" '-DCMAKE_BUILD_TYPE=Release' "-DSP_SANITIZE=$asan"
if ($LASTEXITCODE -ne 0) { throw 'CMake configuration failed.' }
& cmake --build $buildPath --parallel
if ($LASTEXITCODE -ne 0) { throw 'C++ build failed.' }
# Use the real compiler directory: oneAPI's aggregate 2025.2 directory can
# contain zero-length redirection files that Windows cannot load as DLLs.
$compilerRuntimeRoot = Join-Path ${env:ProgramFiles(x86)} 'Intel/oneAPI/compiler/latest'
$runtimeNames = @('libmmd.dll', 'svml_dispmd.dll', 'libirc.dll')
if ($Sanitize) { $runtimeNames += 'clang_rt.asan_dynamic-x86_64.dll' }
foreach ($runtimeName in $runtimeNames) {
    $runtime = Get-ChildItem -LiteralPath $compilerRuntimeRoot -Recurse -Filter $runtimeName -File |
        Where-Object { $_.Length -gt 0 -and $_.FullName -notmatch 'ia32|x86_32' } | Select-Object -First 1
    if ($runtime) { Copy-Item -LiteralPath $runtime.FullName -Destination $buildPath -Force }
    elseif ($runtimeName -like 'clang_rt.asan*') { throw "Required AddressSanitizer runtime missing: $runtimeName" }
}
if (-not $NoTest) {
    & ctest --test-dir $buildPath --output-on-failure
    if ($LASTEXITCODE -ne 0) { throw 'C++ regression tests failed.' }
}
Write-Host "Executable: $buildPath/sphere_packing.exe"
