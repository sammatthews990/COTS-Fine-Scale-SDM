# run_all_refits.ps1
# Runs all COTS SDM refits and validations sequentially to avoid Out of Memory errors

Write-Host "=== Starting COTS SDM Refit Pipeline ==="
Write-Host "1. Refitting 2026 Clean Model..."
Rscript DataProcessing/refit_fy2526_clean.R

Write-Host "2. Refitting 2026 Reef Guide Model..."
Rscript DataProcessing/refit_fy2526_rg.R

Write-Host "3. Refitting Agnostic Models..."
Rscript DataProcessing/refit_agnostic_fy2526.R

Write-Host "4. Running Comprehensive Validation..."
Rscript DataProcessing/validate_model_comprehensive.R

Write-Host "5. Generating Leaflet Maps..."
Rscript DataProcessing/scratch/prep_lizard_leaflet.R
Rscript -e "region_name <- 'lizard'; source('DataProcessing/scratch/build_lizard_html_map.R')"
Rscript -e "region_name <- 'cairns'; source('DataProcessing/scratch/build_lizard_html_map.R')"
Rscript -e "region_name <- 'townsville'; source('DataProcessing/scratch/build_lizard_html_map.R')"
Rscript -e "region_name <- 'whitsunday'; source('DataProcessing/scratch/build_lizard_html_map.R')"

Write-Host "=== Pipeline Complete ==="
