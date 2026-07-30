param(
    [ValidateSet("smoke", "full")]
    [string]$Profile = "smoke",
    [string]$Steps = "01,02,03,04,05,06,07,08,09,10,11,12"
)

$ErrorActionPreference = "Stop"
$projectPath = (Get-Location).Path
$image = "ghcr.io/storopoli/cmdstanr:latest"

docker info | Out-Null
docker pull $image
docker run --rm `
    --mount "type=bind,source=$projectPath,target=/work" `
    --workdir /work `
    --env "MEMOIRE_PROFILE=$Profile" `
    $image `
    bash -lc "Rscript scripts/00_setup.R && Rscript scripts/00_run_pipeline.R --profile=$Profile --steps=$Steps"
