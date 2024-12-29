# Gene panel sequencing - Variant Discovery

This repository contains a complete code pipeline developed to call variants from gene panel sequencing data.

1. First please follow this [Prerequisites Preparation](./1.Preparation ) to setup necessary tools and databases for this practice.

2. Variant calling and filtering workflow can be found [Workflow](./2.Workflow).

3. **Note**
- Before you run this code, pleas change the path in main_code.sh file
- The basic syntax to run

```bash
bash /path/to/main_code.sh \
    -i /path/to/sample_sheet.csv \
    -o /path/to/output \
    -a /path/to/adapters_file \
    -r /path/to/references \
    -w /path/to/Gene_Panel_for_covered.bed
```
