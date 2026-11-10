param(
    [Parameter(Mandatory=$true)][string]$InputDocx,
    [Parameter(Mandatory=$true)][string]$OutputPdf,
    [switch]$UpdateFields,
    [switch]$SaveDocx,
    [string]$PageMapJson
)
$ErrorActionPreference = 'Stop'
$inputPath = (Resolve-Path -LiteralPath $InputDocx).Path
$pdfPath = [IO.Path]::GetFullPath($OutputPdf)
[IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($pdfPath)) | Out-Null
$word = $null
$doc = $null
try {
    # Create a separate automation instance; never attach to the user's Word.
    $word = New-Object -ComObject Word.Application
    if ($word.Documents.Count -ne 0) { throw 'Expected an empty, separate Word automation instance.' }
    $word.Visible = $false
    $word.DisplayAlerts = 0
    $word.AutomationSecurity = 3
    $doc = $word.Documents.Open($inputPath, $false, (-not $SaveDocx), $false)
    if ($UpdateFields) {
        $doc.Repaginate()
        foreach ($toc in $doc.TablesOfContents) { $toc.Update() }
        foreach ($tof in $doc.TablesOfFigures) { $tof.Update() }
        $null = $doc.Fields.Update()
        # Word can move the TOC's otherwise empty closing-field paragraph to
        # a new page. Keep that invisible terminator compact; leave every
        # visible TOC entry and its template style unchanged.
        foreach ($list in @($doc.TablesOfContents) + @($doc.TablesOfFigures)) {
            $tail = $list.Range.Paragraphs.Last
            $visible = $tail.Range.Text -replace '[\r\n\x07\x13-\x15]', ''
            if ([string]::IsNullOrWhiteSpace($visible)) {
                $tail.Range.Font.Size = 1
                $tail.Range.ParagraphFormat.SpaceBefore = 0
                $tail.Range.ParagraphFormat.SpaceAfter = 0
                $tail.Range.ParagraphFormat.LineSpacingRule = 4
                $tail.Range.ParagraphFormat.LineSpacing = 1
                $tail.Range.ParagraphFormat.KeepWithNext = 0
                $tail.Range.ParagraphFormat.PageBreakBefore = 0
            }
        }
        $doc.Repaginate()
        foreach ($toc in $doc.TablesOfContents) { $toc.UpdatePageNumbers() }
        foreach ($tof in $doc.TablesOfFigures) { $tof.UpdatePageNumbers() }
        # Hidden Word may initially cache NUMPAGES as 1. Force full layout and
        # verify the cache without unlinking the native field.
        $totalPages = [int]$doc.ComputeStatistics(2)
        foreach ($field in $doc.Fields) {
            if ($field.Code.Text.Trim() -match '^NUMPAGES\b') {
                $null = $field.Update()
                if ($field.Result.Text.Trim() -ne [string]$totalPages) {
                    $field.Result.Text = [string]$totalPages
                }
            }
        }
    }
    if ($SaveDocx) { $doc.Save() }
    if ($PageMapJson) {
        $headings = @()
        foreach ($p in $doc.Paragraphs) {
            if ($p.OutlineLevel -le 3) {
                $headings += [pscustomobject]@{
                    text = $p.Range.Text.Trim()
                    level = [int]$p.OutlineLevel
                    page = [int]$p.Range.Information(1)
                    physical_page = [int]$p.Range.Information(3)
                }
            }
        }
        [pscustomobject]@{pages=[int]$doc.ComputeStatistics(2);headings=$headings} |
            ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $PageMapJson -Encoding UTF8
    }
    $doc.ExportAsFixedFormat($pdfPath, 17)
    Write-Output ('Rendered ' + $pdfPath)
} finally {
    if ($null -ne $doc) { $doc.Close(0) }
    if ($null -ne $word) {
        if ($word.Documents.Count -eq 0) { $word.Quit() }
        [Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
    }
}
