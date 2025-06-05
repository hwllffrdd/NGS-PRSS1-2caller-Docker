# NGS.PRSS1-2caller Docker Pipeline

A complete Docker-based pipeline for PRSS1-PRSS2 genetic variant analysis, from FASTQ files to annotated variants with plots.

## Features

- **Complete pipeline**: FASTQ → Alignment → Variant Calling → Annotation → Plotting
- **Proper alignment**: Uses GRCh38+ALT contig to avoid pseudogene mapping issues
- **Enhanced R plotting**: Includes variant position plots with R 4.0+
- **SnpEff annotations**: Biological consequence annotations for all variants

## Requirements

- Docker
- Input FASTQ files (paired-end, named `*_R1*.fastq*` and `*_R2*.fastq*`)

## Quick Start

1. **Build the Docker image:**
   ```bash
   docker build -t prss1-pipeline .
   ```

2. **Prepare your data:**
   - Create a directory for your analysis (e.g., `/path/to/your/analysis/`)
   - Create a `fastq/` subdirectory
   - Place your FASTQ files in the `fastq/` directory

3. **Run the analysis:**
   ```bash
   docker run -it --rm -v "/path/to/your/analysis:/data" prss1-pipeline bash run_complete_pipeline.sh
   ```
   or
   ```bash
   docker run -it --rm -v "$(pwd):/data" prss1-pipeline bash run_complete_pipeline.sh
   ```

## Input Structure

```
your-analysis-directory/
├── fastq/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
```

## Output Structure

```
your-analysis-directory/
├── fastq/                              # Your input files
├── bams/                               # Aligned BAM files
├── references/                         # Downloaded reference genomes
└── prss1_analysis/                     # Main results directory
    ├── prss1_analysis_PRSS.vcf         # Variant calls
    ├── prss1_analysis_PRSS_snpEff_ann.vcf  # Annotated variants
    ├── prss1_analysis_PRSS_snpEff_ann.txt  # Annotation summary
    └── prss1_analysis_PRSS_snpEff_ann.plot.pdf  # Variant plot
```

## Example

```bash
# Build the image
docker build -t prss1-pipeline .

# Run analysis on data in /home/user/my_prss_analysis/
docker run -it --rm -v "$(pwd):/data" prss1-pipeline bash run_complete_pipeline.sh
```

## About

This pipeline implements the NGS.PRSS1-2caller tool in a complete Docker environment, adding FASTQ alignment capabilities and enhanced R plotting. It properly handles the complex PRSS1-PRSS2 locus with its pseudogenes to provide accurate variant calling for pancreatitis-related genetic analysis.

## Citation

If you use this pipeline, please cite the original NGS.PRSS1-2caller paper:

Lou H, Xie B, Wang Y, et al. Improved NGS variant calling tool for the PRSS1–PRSS2 locus. Gut 2022. doi: 10.1136/gutjnl-2022-327203
