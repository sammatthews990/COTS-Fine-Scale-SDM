param([switch]$RefreshReview09)
$ErrorActionPreference = 'Stop'
$report = Split-Path -Parent $PSScriptRoot
$final = Split-Path -Parent $report
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$build = Join-Path $final ('cache/ccip_template/build_' + $stamp)
New-Item -ItemType Directory -Path $build -Force | Out-Null
$template = Join-Path $report 'CCIP Final Report Template - Bridging Projects.docx'
$stagedDocx = Join-Path $build 'technical_report.docx'
$stagedPdf = Join-Path $build 'technical_report.pdf'
Push-Location $report
try {
    $env:LC_ALL = 'English_Australia.utf8'
    quarto render technical_report.qmd --to docx --output technical_report_content.docx
    if ($LASTEXITCODE -ne 0) { throw 'Quarto content render failed.' }
    Move-Item -LiteralPath (Join-Path $report 'technical_report_content.docx') -Destination (Join-Path $build 'content.docx')
    python (Join-Path $PSScriptRoot 'build_ccip_report.py') --template $template --content (Join-Path $build 'content.docx') --output $stagedDocx --audit (Join-Path $build 'assembly_audit.json') --metadata (Join-Path $report 'ccip_report_metadata.json')
    if ($LASTEXITCODE -ne 0) { throw 'Template assembly failed.' }
    & (Join-Path $PSScriptRoot 'render_word_native.ps1') -InputDocx $stagedDocx -OutputPdf $stagedPdf -UpdateFields -SaveDocx -PageMapJson (Join-Path $build 'page_map.json')
    python (Join-Path $PSScriptRoot 'check_ccip_report.py') --template $template --docx $stagedDocx --source (Join-Path $report 'technical_report.qmd') --page-map (Join-Path $build 'page_map.json') --out (Join-Path $build 'checks.json')
    if ($LASTEXITCODE -ne 0) { throw 'Template fidelity/content check failed; staged files retained for inspection.' }
    $backup = Join-Path $build 'previous_exports'
    New-Item -ItemType Directory -Path $backup -Force | Out-Null
    foreach ($name in @('technical_report.docx','technical_report.pdf')) {
        $old = Join-Path $report $name
        if (Test-Path -LiteralPath $old) { Copy-Item -LiteralPath $old -Destination (Join-Path $backup $name) }
    }
    Copy-Item -LiteralPath $stagedDocx -Destination (Join-Path $report 'technical_report.docx') -Force
    Copy-Item -LiteralPath $stagedPdf -Destination (Join-Path $report 'technical_report.pdf') -Force
    if ($RefreshReview09) {
        $reviews = Join-Path $report 'reviews'
        foreach ($name in @('technical_report_review_09.docx','technical_report_review_09_baseline.docx','technical_report_review_09_baseline.qmd')) {
            $old = Join-Path $reviews $name
            if (Test-Path -LiteralPath $old) { Copy-Item -LiteralPath $old -Destination (Join-Path $backup $name) }
        }
        Copy-Item -LiteralPath $stagedDocx -Destination (Join-Path $reviews 'technical_report_review_09.docx') -Force
        Copy-Item -LiteralPath $stagedDocx -Destination (Join-Path $reviews 'technical_report_review_09_baseline.docx') -Force
        Copy-Item -LiteralPath (Join-Path $report 'technical_report.qmd') -Destination (Join-Path $reviews 'technical_report_review_09_baseline.qmd') -Force
    }
    Write-Output ('CCIP report export completed. Build evidence: ' + $build)
} finally { Pop-Location }
