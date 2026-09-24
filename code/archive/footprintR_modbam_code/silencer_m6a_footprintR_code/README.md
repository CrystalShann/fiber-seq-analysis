# Silencer m6A modbam footprintR pipeline

This directory contains a 3-part footprintR modbam workflow for silencer-centered analysis using m6A only.

## Directory

/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/silencer_m6a_footprintR_code

## Files

- silencer_modbam_summary_part1.R: Build silencer windows and nucleosome-aligned summaries.
- silencer_modbam_summary_part2.R: Run chunked modbam extraction and aggregate m6A accessibility summaries and metaprofiles.
- silencer_modbam_summary_part3.R: Plot silencer-centered m6A metaprofiles (optionally with nucleosome structure panel).
- submit_silencer_part1.sh: Slurm launcher for part 1.
- submit_silencer_part2.sh: Slurm launcher for part 2.
- submit_silencer_part3.sh: Slurm launcher for part 3.
- submit_silencer_all.sh: Submit parts 1 to 3 with dependencies for one silencer class.
- submit_silencer_4classes.sh: Submit full part1->part2->part3 chain for 4 requested silencer classes.

## Silencer classes (4)

- REST-Enhancers
- REST-Silencers
- STARR-Silencers.Robust
- STARR-Silencers.Stringent

## Current output root path

Outputs are now split by scope:

- all chromosome mode: /project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/m6a_modbam/all_chr
- single chromosome mode: /project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/m6a_modbam/sample_chr

Each scope stores results as:

- threshold_0.9/<SAMPLE>/<CLASS_TAG>/...

## Chromosome behavior

- CHROM=ALL in submit_silencer_all.sh submits one part1->part2->part3 chain per chromosome.
- Default CHROM in submit_silencer_all.sh and submit_silencer_4classes.sh is ALL.
- Chromosome list for ALL mode defaults to chr1..chr22,chrX,chrY and can be overridden using ALL_CHRS.

## Legacy note

Earlier versions wrote under a single root:

/project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/m6a_modbam

## Slurm resources (current)

- part1: 24:00:00, 300G
- part2: 24:00:00, 300G
- part3: 24:00:00, 300G

## Submission examples

Submit one class chain:

SAMPLE=AL10_bc2178_19130 CHROM=chr1 SILENCER_CLASS=REST-Silencers bash submit_silencer_all.sh

Submit all 4 classes:

SAMPLE=AL10_bc2178_19130 CHROM=chr1 bash submit_silencer_4classes.sh

Submit all 4 classes across all chromosomes (default behavior):

SAMPLE=AL10_bc2178_19130 CHROM=ALL bash submit_silencer_4classes.sh

## Change log

- 2026-04-23
  - Created silencer m6A-only pipeline scripts (part1, part2, part3) and Slurm wrappers.
  - Added class-aware submission wrapper for all 4 requested silencer classes.
  - Updated Slurm resources for all parts to time=24:00:00 and mem=300G.
  - Updated results output root to /project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/m6a_modbam.
  - Added CHROM=ALL submit fan-out and output scope split:
    - all_chr path when CHROM=ALL
    - sample_chr path when CHROM is a single chromosome
