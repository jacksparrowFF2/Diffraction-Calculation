param(
    [string]$PackageName = 'COEReflection_handoff_20261112'
)
$ErrorActionPreference = 'Stop'
if ($PackageName -notmatch '^[a-zA-Z0-9_-]+$') { throw 'PackageName must be a simple file name.' }
$moduleDir = Split-Path -Parent $PSScriptRoot
$repositoryDir = Split-Path -Parent $moduleDir
$stageDir = Join-Path ([IO.Path]::GetTempPath()) ('coe_handoff_' + [guid]::NewGuid().ToString('N'))
$bundleDir = Join-Path $stageDir $PackageName
$copiedModuleDir = Join-Path $bundleDir 'topic_coe_reflection'
New-Item -ItemType Directory -Path $copiedModuleDir -Force | Out-Null
foreach ($part in @('src', 'Input', 'docs')) {
    Copy-Item -LiteralPath (Join-Path $moduleDir $part) -Destination (Join-Path $copiedModuleDir $part) -Recurse
}
New-Item -ItemType Directory -Path (Join-Path $copiedModuleDir 'Output') -Force | Out-Null
$copiedDistDir = Join-Path $copiedModuleDir 'dist'
New-Item -ItemType Directory -Path $copiedDistDir -Force | Out-Null
Copy-Item -LiteralPath $PSCommandPath -Destination (Join-Path $copiedDistDir 'BuildHandoff.ps1')
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'validate_handoff.m') -Destination (Join-Path $copiedDistDir 'validate_handoff.m')
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'validate_handoff.m') -Destination (Join-Path $bundleDir 'validate_handoff.m')
Copy-Item -LiteralPath (Join-Path $moduleDir 'docs\程序使用说明.md') -Destination (Join-Path $bundleDir 'README.md')
$testsDir = Join-Path $bundleDir 'tests'
New-Item -ItemType Directory -Path $testsDir -Force | Out-Null
foreach ($testFile in @('COESpectralTest.m', 'COESpotsTest.m')) {
    Copy-Item -LiteralPath (Join-Path $repositoryDir ('tests\' + $testFile)) -Destination (Join-Path $testsDir $testFile)
}
$referenceDir = Join-Path $bundleDir 'reference_outputs'
New-Item -ItemType Directory -Path $referenceDir -Force | Out-Null
Get-ChildItem -LiteralPath (Join-Path $moduleDir 'Output') -File |
    Where-Object { $_.Name -like 'coe_reflection_*' -and $_.Extension -in @('.csv', '.png') } |
    ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $referenceDir $_.Name) }
$manifest = Get-ChildItem -LiteralPath $bundleDir -Recurse -File | Sort-Object FullName | ForEach-Object {
    [pscustomobject]@{
        RelativePath = $_.FullName.Substring($bundleDir.Length + 1)
        Bytes = $_.Length
        SHA256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash
    }
}
$manifest | Export-Csv -LiteralPath (Join-Path $bundleDir 'MANIFEST_SHA256.csv') -NoTypeInformation -Encoding UTF8
$archivePath = Join-Path $PSScriptRoot ($PackageName + '.zip')
Compress-Archive -LiteralPath $bundleDir -DestinationPath $archivePath -Force
$archiveHash = (Get-FileHash -LiteralPath $archivePath -Algorithm SHA256).Hash
Set-Content -LiteralPath ($archivePath + '.sha256') -Value ($archiveHash + '  ' + (Split-Path -Leaf $archivePath)) -Encoding ASCII
[pscustomobject]@{
    Archive = $archivePath
    Bytes = (Get-Item -LiteralPath $archivePath).Length
    SHA256 = $archiveHash
    StagingFolder = $bundleDir
    Files = @($manifest).Count
} | ConvertTo-Json
