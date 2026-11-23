To report a problem or suggest a change:

1. Open an issue describing the problem or idea before larger changes.
2. Fork the repository, make a branch, and check that `pixi run reproduce`
   still works on your computer.
3. Open a pull request. The *Reproduce* check must pass before merging.



Please keep raw data out of Git (`data/raw/` is downloaded by `scripts/download_data.py`), keep analysis settings in `config/config.yaml`, commit `pixi.toml` and `pixi.lock` together, and record analysis decisions in `docs/decisions.md`.


Everyone taking part is expected to follow the [Code of Conduct](CODE_OF_CONDUCT.md).