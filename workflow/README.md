# workflow


`Snakefile` runs the four analysis steps in order (scripts 01 to 04) and only
reruns a step when its inputs, script or settings change.

Run it with `pixi run run-snakemake`, or the whole project (download, analysis
and report) with `pixi run reproduce`. Use `pixi run snakemake --cores 1 -n`
to see which steps would run without running them.

Settings are read from `config/config.yaml`; logs for each step are written to
`results/logs/`. 