# Comparative Transcriptomic Analysis of the Effect of an AKT inhibitor in Ductal Adenocarcinoma Cell Lines revealed Pharmacoresistance Mechanisms

## Authors

Rebeca González-Ortiz¹,  
Leticia Vega-Alvarado²,  
Verónica Jiménez-Jacinto¹,  
Magali Espinosa³,  
Georgina Hernández-Montes¹⁴,  
Jorge Melendez-Zajgla³,  
Alejandro Sanchez-Flores¹*

**Affiliations**

1. Unidad Universitaria de Secuenciación Masiva y Bioinformática, Instituto de Biotecnología, UNAM, Cuernavaca, Morelos, Mexico  
2. Instituto de Ciencias Aplicadas y Tecnología, Universidad Nacional Autónoma de México, Ciudad de México, Mexico  
3. Laboratorio de Genómica Funcional del Cáncer, Instituto Nacional de Medicina Genómica (INMEGEN), Ciudad de México, Mexico  
4. Red de Apoyo a la Investigación, UNAM – Instituto Nacional de Ciencias Médicas y Nutrición Salvador Zubirán, Ciudad de México, Mexico  

---

## Repository Overview

This repository contains the data, analysis outputs, and figures associated with the manuscript:

**“Comparative Transcriptomic Analysis of the Effect of an AKT inhibitor in Ductal Adenocarcinoma Cell Lines revealed Pharmacoresistance Mechanisms.”**

The study analyzes RNA-Seq data from two ductal adenocarcinoma cell lines (LNCaP and BxPC-3) treated with the AKT inhibitor AZD5363, with the aim of identifying differential expression patterns and potential molecular mechanisms of drug resistance.

---

## Contents

### Quality Control Reports

FastQC quality control outputs for all sequenced samples:

- `BX_AZ1_1_fastqc.html`
- `BX_AZ1_2_fastqc.html`
- `BX_C1_1_fastqc.html`
- `LNCAP_AZ1_1_fastqc.html`
- `LNCAP_C1_1_fastqc.html`
- *(and additional reports for all replicates)*

### Differential Expression Results

Processed result files for functional enrichment and differential expression analyses:

- `BxPC-3_UP_Project_wg_result1750714681.zip`
- `BxPC-3_DOWN_Project_wg_result1750714705.zip`
- `LNCaP_UP_Project_wg_result1750714241.zip`
- `LNCaP_DOWN_Project_wg_result1750714262.zip`

### Visualizations

- `heatmap_diferencial_BxPC3_LNCaP.pdf` – final comparative heatmap
- `heat_map.txt` – data matrix used for heatmap generation

### Supporting Data

- `BxPC-3_TDC.txt`
- `LNCAP_TDC.txt`
- `Cada dato del trabajo de doctorado.xlsx`
- `code.md` – notes and commands used during analysis

---

## Objective

The purpose of this repository is to provide:

- Transparency and reproducibility of RNA-Seq analyses  
- Access to processed data and figures  
- Supporting material for peer review and publication  

---

## How to Use

The repository can be used to:

1. Inspect sequencing quality (FastQC reports)  
2. Review differential expression outputs  
3. Reproduce final figures and tables  
4. Validate results presented in the manuscript  

Scripts for full reproducibility will be added in future updates.

---

## License

This project is licensed under the Apache 2.0 License.

---

## Citation

If you use the data or materials in this repository, please cite the associated manuscript and this repository using the information provided in `CITATION.cff`.
