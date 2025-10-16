# Tissue of Tumor Origin 

This repository contains the code for the manuscript "5-hydroxymethylcytosine Analysis Reveals Stable Epigenomic Changes in Tumor Tissue that Enable Cancer Detection in Cell-free DNA" (2025), Xue, Ning et al. 

A link to the paper will be provided once published.

## Clone the repository

```bash
git clone https://github.com/ClearNoteHealth-OpenAccess/ClearNoteHealth_TOTO_nature_comm_bio_analysis_code
cd ClearNoteHealth_TOTO_nature_comm_bio_analysis_code
```

## Data

Data is available at this repository under `data/` directory.

## Usage

To build the docker image, run the following command:

```bash
docker build -t cnh_toto .
```

Our analysis results and figures are generated using a combination of python and R scripts.

### Python

To reproduce the python analysis results, run the following commands:

```bash
#### Main Figures
docker run -it -v $PWD/:/cnh/ --rm cnh_toto python python/manuscripts/Main_Figures.py

#### Supplemental Figures
docker run -it -v $PWD/:/cnh/ --rm cnh_toto python python/manuscripts/Supplemental_Figures.py
```

### R

To reproduce the R analysis results, run the following commands:

```bash
#### Main Figures
docker run -it -v $PWD/:/cnh/ --rm cnh_toto Rscript R/manuscripts/Main_Figures.R

#### Supplemental Figures
docker run -it -v $PWD/:/cnh/ --rm cnh_toto Rscript R/manuscripts/Supplemental_Figures.R
```

## Citation

If you use the data or code in this repository, please cite:

Xue, Y., Ning, Y., Friedl, V., Haan, D., Bergamaschi, A., Guler, G. D., Hazen, K., Scott, A., Phillips, T., McCarthy, E., Ellison, C. K., Malta, R., Nguyen, A., Lopez, V., Cavet, R., Peters, M., Sheard, J. Gibb, W., Chowdhury, S., Volkmuth, W., & Levy, S. (2025). 5-hydroxymethylcytosine Analysis Reveals Stable Epigenomic Changes in Tumor Tissue that Enable Cancer Detection in Cell-free DNA (Version 1.0.0) [Computer software]. https://github.com/ClearNoteHealth-OpenAccess/ClearNoteHealth_TOTO_nature_comm_bio_analysis_code
