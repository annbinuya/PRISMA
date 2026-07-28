# External validation of BOADICEA in the PRISMA screening cohort

This repository contains the R code used to describe the study population and externally validate BOADICEA version 7.3.2 in the PRISMA cohort.

PRISMA is a population-based cohort of approximately 82,200 women participating in the Dutch breast cancer screening program. The validation used a nested case-control sample comprising 318 women diagnosed with invasive breast cancer and 771 controls.

Five-year breast cancer risks were estimated using models that progressively included:

* Age
* Established clinical and lifestyle risk factors, family history, and available pathogenic variant information
* Volpara-based mammographic breast density
* The 313-variant polygenic risk score
* Combinations of breast density and the polygenic risk score

Model performance was assessed using discrimination, calibration, risk reclassification, and decision curve analysis. The main decision context was the use of five-year risk estimates to inform consideration of risk-reducing medication.

## Code / Syntax files
| File                   | Description             |
| :----                  | :----                   |
| `1_Descriptive tables.Rmd` | Prepares the analysis datasets and generates descriptive tables and distributions for participant characteristics, outcomes, breast density, the 313-variant polygenic risk score, and predicted BOADICEA risks. The selected analysis dataset can be changed to run the main or sensitivity analyses.                |
| `2_Validation.Rmd`         | Evaluates five-year BOADICEA performance, including the C-index, observed-to-expected ratio, calibration slope, calibration plots, net reclassification improvement, and net benefit at clinically relevant risk thresholds. The selected analysis dataset can be changed to run the main or sensitivity analyses.                |

## Data availability

The individual-level PRISMA data, cleaning, and intermediate R workspaces are not included because they contain restricted research data.

## Contact

Mary Ann E. Binuya <br/>
Netherlands Cancer Institute <br/>
[m.binuya@nki.nl](m.binuya@nki.nl)


## Co-authors

* Mary Ann E. Binuya
* Jim Peters
* Renée Verdiesen
* Maartje Schreurs
* Daniëlle van der Waal
* Marjanka K. Schmidt
* Mireille J.M. Broeders
