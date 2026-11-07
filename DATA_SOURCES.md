# Data sources and version boundaries

## PBC

The extension reads `survival::pbc` from `survival` 3.8-3 in the recorded environment. It contains 418 records: 312 randomized participants plus 106 additional followed participants. Endpoint codes are **0 = censored, 1 = transplantation, 2 = death**. Age is measured in **years**, follow-up in **days**. These units differ from some Kaggle encodings.

Source documentation: [survival PBC data](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pbc.html). Bibliographic source: Therneau and Grambsch (2000), *Modeling Survival Data: Extending the Cox Model*. The primary biliary cholangitis terminology is current; the original course title used cirrhosis.

The course HTML read a local `cirrhosis.csv` that is missing from the supplied directory. The documented PBC data are a related source version, not a verified exact copy of that CSV. The original C-index 0.8574 belongs to the course split and is not a new validated result.

