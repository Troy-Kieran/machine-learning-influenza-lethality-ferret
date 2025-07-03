## Overview

R code used for the machine learning portion of: 

Kieran TJ, Sun X, Maines TR, Belser JA. Predictive models of influenza A virus lethal disease yield insights from ferret respiratory tract and brain tissues. (Accepted June 25: In Press, Scientific Reports, 2025). 

This project also includes CSV files for summary data & metrics used to create figures from the manuscript in the "ResultsSummary_Rinputs" directory. 

## References and Resources

Study makes use of previously published data:

Kieran TJ, Sun X, Creager HM, Tumpey TM, Maines TR, Belser JA. An aggregated dataset of serial morbidity and titer measurements from influenza A virus-infected ferrets. Sci Data 11, 510 (2024). https://doi.org/10.1038/s41597-024-03256-6

CDC. An aggregated dataset of serially collected influenza A virus morbidity and titer measurements from virus-infected ferrets, <https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-serially-collected-influe/cr56-k9wj/about_data> (2025).

Kieran TJ, Sun X, Tumpey TM, Maines TR, Belser JA. Spatial variation of infectious virus load in aggregated day 3 post-inoculation respiratory tract tissues from influenza A virus-infected ferrets. Under peer review.

CDC. An aggregated dataset of day 3 post-inoculation viral titer measurements from influenza A virus-infected ferret tissues, <https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-day-3-post-inoculation-vi/d9u6-mdu6/about_data> (2025).

Previous related publications:

Kieran TJ, Sun X, Maines TR, Belser JA. Machine learning approaches for influenza A virus risk assessment identifies predictive correlates using ferret model in vivo data. Commun Biol 7, 927 (2024). https://doi.org/10.1038/s42003-024-06629-0

https://github.com/CDCgov/machine-learning-influenza-ferret-model/tree/main

Kieran TJ, Sun X, Maines TR, Beauchemin CAA, Belser JA. Exploring associations between viral titer measurements and disease outcomes in ferrets inoculated with 125 contemporary influenza A viruses. J Virol 98:e01661-23. (2024). (PMID 38240592)
([https://doi.org/10.1128/jvi.01661-23](https://doi.org/10.1128/jvi.01661-23))

## Scripts
lethality_TJK.R script is the machine learning code for lethality analysis.

morbidity_TJK.R script is the machine learning code for morbidity analysis.

## Note
R script uses predicted sialic acid (receptor) binding preference (RBS) and predicted polymerase activity (PBS) markers (based on molecular sequence), along with selected HA and PB2 gene markers, that are not included in the CSV file from Sci Data, but may be cross-referenced as reported in Kieran et al (PMID 38240592).

## Manuscript Abstract
Collection of systemic tissues from influenza A virus (IAV)-infected ferrets at a fixed timepoint post-inoculation represents a frequent component of risk assessment activities to assess the capacity of IAV to replicate systemically. However, few studies have evaluated how the frequency and magnitude of IAV replication at discrete tissues contribute to within-host phenotypic outcomes, limiting our ability to fully contextualize results from scheduled necropsy into risk assessment settings. Employing aggregated data from ferrets inoculated with >100 unique IAV (both human- and avian-origin viruses, spanning H1, H2, H3, H5, H7, and H9 subtypes), we examined relationships between infectious virus detection in four discrete tissue types (nasal turbinate, lung, brain, and olfactory bulb [BnOB]) to clinical outcomes of IAV-inoculated ferrets, and the utility of including these discrete tissue data as features in machine learning (ML) models. We found that addition of viral tissue titer data maintained high performance metrics of a predictive lethality classification ML model with or without inclusion of serially-collected virological and clinical data. Interestingly, infectious virus in BnOB was detected at higher frequency and magnitude among IAV associated with high pathogenicity phenotypes in ferrets, more so than tissues from the respiratory tract; in agreement, BnOB was the highest relative ranked individual tissue specimen in predictive classification models. This study highlights the potential role of BnOB viral titers in assessing IAV pathogenicity in ferrets, and highlights the role ML approaches can contribute towards understanding the predictive benefit of in vivo-generated data in the context of pandemic risk assessment.

##
##
##
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).
