# Building a GenBank Manuscript Summary Generator with Python and Tkinter

## Introduction

As a bioinformatician working with genomic data, I often find myself needing to summarize GenBank files for manuscripts and reports. To streamline this process, I created a simple desktop application using Python, Tkinter, and Biopython that automatically generates manuscript-style summaries of GenBank files.

## Key Features

1. **User-Friendly Interface**: A clean, simple GUI built with Tkinter that allows users to:
   - Load GenBank files through a file dialog
   - View generated summaries in a scrollable text area
   - Clear the output with a single click

2. **Smart Gene Detection**: The tool automatically identifies and categorizes:
   - Protein-coding sequences (CDS)
   - Transfer RNA (tRNA) genes
   - Ribosomal RNA (rRNA) genes

3. **Manuscript-Ready Output**: Generates a three-paragraph summary including:
   - Sequence statistics (length and GC content)
   - Gene count breakdown
   - Gene synteny description

## Implementation Highlights

### Gene Information Extraction

The `extract_gene_info` function processes GenBank features and smartly handles different gene types:

```python
def extract_gene_info(record):
    counts = {'CDS':0, 'tRNA':0, 'rRNA':0}
    genes = []
    for feat in record.features:
        ft = feat.type
        if ft not in counts:
            continue
        counts[ft] += 1
        # Smart name extraction with fallback options
        name = (feat.qualifiers.get('gene')
                or feat.qualifiers.get('product')
                or [ft])[0]
