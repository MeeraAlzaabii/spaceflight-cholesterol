# containers


Container recipes for this project's software environment, built from
`pixi.lock` so the container runs exactly the same packages as `pixi install`:

- `Dockerfile`: Docker image.
- `apptainer.def`: the same image for HPC clusters that use Apptainer.

Build the image from the project folder:

    docker build --platform linux/amd64 -f containers/Dockerfile -t spaceflight-cholesterol:1.0.0 .

Run the whole analysis inside it (download, Snakemake, report):

    docker run --rm -v "$PWD":/work -w /work spaceflight-cholesterol:1.0.0 \
      bash -c "python scripts/download_data.py && Rscript scripts/fetch_kegg_gene_sets.R && snakemake --cores 1 && quarto render reports/report.qmd"

Results appear in `results/` and `reports/report.html`, as with `pixi run reproduce`.