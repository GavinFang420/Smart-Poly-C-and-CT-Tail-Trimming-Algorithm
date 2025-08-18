## Evaluation Framework
_Inspired by the methods used in fastp (**Chen et al., Bioinformatics, 2018**) to systematically assess our smart polyC-tail trimming approach (comparing to traditonal hard trimming methods)._

<br>

**1. Data Retention & Sequence Integrity**
- Read length distribution:
  - Compare histograms of read lengths after trimming.
- Measure total_bases_after_trimming / total_bases_before.
- Quality score profiles:
  - Plot mean Phred score per cycle (FastQC-style)
  - Check that trimming does not disproportionately remove high-quality bases in R1.
> Smart polyC trimming reduces artificial sequence while preserving more R1 info than hard clipping.

<br>

**2. Base Composition**
- Per-base nucleotide composition:
  - Before trimming, R2 starts with artificially high %C (polyC).
  - After trimming, composition should flatten, aligning closer to expected background (~0.4% C in WGBS).
- K-mer bias analysis:
  - Count overrepresented k-mers (e.g., “CCCCCCCC”) before/after trimming.
  - Smart polyC trimming should remove polyC k-mers more efficiently than fixed trimming.

<br>

**3. Mapping Quality & Alignment Efficiency**
- Mapping rate:
  - Run Bismark or Bowtie2, measure % mapped reads.
  - Expect Smart polyC trimming > fixed trimming, especially for non-overlapping R1/R2.
- Soft-clipped bases in alignments:
  - Count total soft-clip bases in BAM files.
  - Hard trimming may mask issues, while Smart polyC trimming should reduce unnecessary soft-clips by removing only artificial tails.
 
<br>

**4. Duplication rate**
- Typically by looking at read start/end positions
  - Hard trimming may falsely merge distinct reads (e.g., 148bp+2C vs 145bp+5C).
- Smart polyC trimming should better preserve true uniqueness.

<br>

**5. Speed & Efficiency**
- Runtime & memory footprint:
  - Measure trimming throughput (MB/s).
  - Smart polyC trimming method should remain linear-time like fastp (O(n)).
