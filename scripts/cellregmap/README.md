# CellRegMap for OneK1K
This repository contains the code to run [CellRegMap](https://www.biorxiv.org/content/10.1101/2021.09.01.458524v1) on the [OneK1K](https://www.garvan.org.au/research/garvan-weizmann/research) / [TOB-WGS](https://github.com/populationgenomics/tob-wgs) dataset.

To run a first test, use conda/mamba to install the analysis-runner, then execute the following command:

```sh
analysis-runner \
   --dataset tob-wgs --access-level test \
   --description 'first run of CellRegMap using small test inputs' \
   --output-dir '2021-12-20-acuomo-testing/ar-test' \
   python cellregmap_basic_usage.py
 ```

## Software
CellRegMap software is available at https://github.com/limix/CellRegMap/.
