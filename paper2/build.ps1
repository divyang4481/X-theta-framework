param(
  [string]$Main = "main.tex",
  [switch]$Fast,
  [switch]$NoClean
)

$ErrorActionPreference = "Stop"
$OutDir = "build"

# --- Prep output dir
if (!(Test-Path $OutDir)) { New-Item -ItemType Directory -Path $OutDir | Out-Null }

# --- Clean old aux (optional)
if (-not $NoClean) {
  $patterns = @(
    '*.aux','*.log','*.bbl','*.blg','*.bcf','*.run.xml','*.out','*.toc','*.pdf',
    '*.fdb_latexmk','*.fls','*.synctex.gz'
  )
  $pathsToRemove = $patterns | ForEach-Object { Join-Path $OutDir $_ }
  Remove-Item -Path $pathsToRemove -Force -ErrorAction SilentlyContinue
}

# --- Ensure figures are available in build (optional convenience)
$FigsSrc = Join-Path $PSScriptRoot "figs"
$FigsDst = Join-Path $OutDir "figs"
if (Test-Path $FigsSrc) {
  Copy-Item -Path $FigsSrc -Destination $FigsDst -Recurse -Force
}

# --- Ensure .bib files are available in build (so bibtex can find them)
$BibSrcFiles = Get-ChildItem -Path $PSScriptRoot -Filter *.bib -File -ErrorAction SilentlyContinue
if ($BibSrcFiles) {
  foreach ($bib in $BibSrcFiles) {
    Copy-Item -Path $bib.FullName -Destination $OutDir -Force -ErrorAction SilentlyContinue
  }
}

function Exists($cmd) { $null -ne (Get-Command $cmd -ErrorAction SilentlyContinue) }

# Helpers
$MainBase = [IO.Path]::GetFileNameWithoutExtension($Main)
$PdfPath  = Join-Path $OutDir ($MainBase + ".pdf")
$LogPath  = Join-Path $OutDir ($MainBase + ".log")
$AuxPath  = Join-Path $OutDir ($MainBase + ".aux")
$BblPath  = Join-Path $OutDir ($MainBase + ".bbl")
$BcfPath  = Join-Path $OutDir ($MainBase + ".bcf")  # biblatex indicator

# When running inside $OutDir, use local file names (not build/...).
$AuxName = $MainBase + ".aux"
$BcfName = $MainBase + ".bcf"

function Run-Bib {
  Push-Location $OutDir
  if ((Test-Path $BcfName) -and (Exists "biber")) {
    Write-Host ">> biber $MainBase"
    biber $MainBase
  } elseif ((Test-Path $AuxName) -and (Exists "bibtex")) {
    Write-Host ">> bibtex $MainBase"
    bibtex $MainBase
  } else {
    Write-Warning "No AUX/BCF found or no bib tool available; skipping bibliography step."
  }
  Pop-Location
}

# --- Build
if (!$Fast -and (Exists "latexmk")) {
  Write-Host ">> latexmk build"
  latexmk -pdf -interaction=nonstopmode -halt-on-error $Main
} else {
  Write-Host ">> classic build (pdflatex -> bib -> pdflatex x2)"
  pdflatex -interaction=nonstopmode -file-line-error -output-directory $OutDir $Main
  Run-Bib
  pdflatex -interaction=nonstopmode -file-line-error -output-directory $OutDir $Main
  pdflatex -interaction=nonstopmode -file-line-error -output-directory $OutDir $Main
}

# --- Diagnostics
if (Test-Path $PdfPath) {
  Write-Host "`n✅ PDF ready: $PdfPath"
} else {
  Write-Warning "⚠️ PDF not found; see $LogPath"
}

# Warn if bibliography likely missing
if (Test-Path $AuxPath) {
  $aux = Get-Content $AuxPath -ErrorAction SilentlyContinue | Out-String
  $hasCites = $aux -match '\\citation\{'
  if ($hasCites -and -not (Test-Path $BblPath)) {
    Write-Warning "Citations detected but no $($MainBase).bbl produced. Bibliography step may have failed."
  }
}

if (Test-Path $LogPath) {
  $log = Get-Content $LogPath -ErrorAction SilentlyContinue

  $undefCites = $log | Select-String "Citation `.*' on page|Warning--I didn't find a database entry" -SimpleMatch
  $undefRefs  = $log | Select-String "LaTeX Warning: Reference.*undefined" -SimpleMatch
  $missingFig = $log | Select-String "File `.*' not found: using draft setting" -SimpleMatch

  if ($undefCites) {
    Write-Warning "Undefined citations:"
    $undefCites | ForEach-Object { $_.Line }
  }
  if ($undefRefs) {
    Write-Warning "Undefined cross-refs:"
    $undefRefs | ForEach-Object { $_.Line }
  }
  if ($missingFig) {
    Write-Warning "Missing figures (LaTeX used draft placeholder):"
    $missingFig | ForEach-Object { $_.Line }
  }
}

Write-Host "`nTips:"
Write-Host "  - Use: .\build.ps1              # full build (runs bibtex/biber)"
Write-Host "  - Use: .\build.ps1 -Fast        # quick pass (no bib tool if latexmk missing)"
Write-Host "  - Use: .\build.ps1 -NoClean     # keep previous aux/bbl/log"
