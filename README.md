# pangraph

> a bioinformatic toolkit to align large sets of genome assemblies into a graph data structure

### install requirements

- python >= 3.8
- pipenv
- mash
- minimap2
- GNUmake (optional for build)

#### manual
run `pipenv sync` to pull all dependencies into a virtual environment
run `./bin/setup_pangraph` to download mash && minimap2 binaries [optional-linux only]
append `./bin` to your system path

### build

the repo contains a Makefile that will generate and analyze synthetic data.
the input is the targets file, each simulation run is given by a line entry.
the synthetic data can be generated (e.g. with 4 cores) by running `make -j4`

### license

[MIT License](LICENSE)
