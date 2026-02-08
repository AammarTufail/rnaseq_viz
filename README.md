# 🧬 RNA-seq DESeq Volcano Plot Analyzer

An interactive Streamlit application for visualizing differential gene expression results from DESeq analysis (READemption pipeline output).

## Features

- 📊 Interactive volcano plot visualization using Plotly
- 🎯 Customizable filtering criteria (padj and log2FoldChange thresholds)
- 🎨 Customizable color scheme for data points
- 🔍 Hover information showing gene names and locus tags
- 💾 Export filtered data as CSV
- 📸 Export plots as PNG and SVG (publication-ready)
- 📈 Real-time statistics of differentially expressed genes
- 📋 Tabular view of upregulated and downregulated genes

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Setup Steps

1. **Navigate to the application directory:**
   ```bash
   cd /media/codanics/ext_ssd/01_rnaseq_paired_end/application_rnaseq
   ```

2. **Run the setup script:**
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```

3. **Activate the virtual environment:**
   ```bash
   source venv/bin/activate
   ```

### Manual Setup (Alternative)

If the setup script doesn't work, follow these steps:

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate  # On Linux/Mac
# OR
venv\Scripts\activate  # On Windows

# Install dependencies
pip install -r requirements.txt
```

## Usage

1. **Start the application:**
   ```bash
   streamlit run app.py
   ```

2. **Upload your data:**
   - Click "Browse files" in the sidebar
   - Upload your DESeq output file (CSV/TSV format)
   - Select the appropriate file separator

3. **Configure filters:**
   - Set adjusted p-value threshold (default: 0.05)
   - Set log2FoldChange threshold for upregulated genes (default: 1.0)
   - Set log2FoldChange threshold for downregulated genes (default: -1.0)

4. **Customize visualization:**
   - Choose colors for upregulated, downregulated, and non-significant genes
   - Hover over data points to see gene details

5. **Export results:**
   - Download filtered data as CSV
   - Download volcano plot as PNG or SVG

## Input File Format

Your DESeq output file should contain the following columns:

### Required Columns:
- `log2FoldChange`: Log2 fold change values
- `padj`: Adjusted p-values (Benjamini-Hochberg correction)

### Optional Columns:
- `Attributes`: Gene annotations containing:
  - `gene=<gene_name>`
  - `locus_tag=<locus_tag>`

### Example:

```csv
gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,Attributes
gene_001,245.3,2.5,0.3,8.3,1.2e-16,3.4e-15,gene=gyrA;locus_tag=PA_0001;product=DNA gyrase
gene_002,156.7,-1.8,0.25,-7.2,6.5e-13,9.1e-12,gene=rpoB;locus_tag=PA_0002;product=RNA polymerase
```

## Features Explained

### Volcano Plot
- **X-axis**: log2 Fold Change
- **Y-axis**: -log10(adjusted p-value)
- **Colors**: 
  - Red: Upregulated genes (default)
  - Blue: Downregulated genes (default)
  - Gray: Non-significant genes (default)
- **Dashed lines**: Threshold cutoffs
- **Interactive**: Hover to see gene details

### Statistics Dashboard
- Total number of genes analyzed
- Count of upregulated genes
- Count of downregulated genes
- Count of non-significant genes

### Data Tables
- Separate tabs for upregulated and downregulated genes
- Sorted by fold change magnitude
- Display gene names, locus tags, fold changes, and adjusted p-values

## Suggested Future Features

### 1. **Enhanced Statistical Analysis**
   - MA plot (log2FC vs mean expression)
   - P-value distribution histogram
   - Fold change distribution plot

### 2. **Gene Set Enrichment**
   - GO term enrichment analysis
   - KEGG pathway analysis
   - Custom gene set analysis

### 3. **Multiple Comparison Support**
   - Compare multiple DESeq results
   - Venn diagrams for overlapping genes
   - Heatmap of selected genes

### 4. **Advanced Filtering**
   - Filter by gene name/locus tag
   - Filter by base mean expression
   - Multiple threshold combinations

### 5. **Data Integration**
   - Link to genome browsers
   - Link to gene databases (NCBI, UniProt)
   - Integration with STRING for protein-protein interactions

### 6. **Reporting**
   - Generate PDF reports
   - Automated statistical summaries
   - Custom annotations for genes

### 7. **Plot Customization**
   - Adjust marker sizes
   - Add gene labels for top DEGs
   - Multiple plot layouts
   - Custom axis ranges

### 8. **Batch Processing**
   - Process multiple files at once
   - Automated report generation
   - Comparison across experiments

### 9. **Data Quality Control**
   - PCA plots
   - Sample correlation heatmaps
   - Dispersion plots

### 10. **Export Options**
   - Export to R-compatible formats
   - Generate code for reproducibility
   - Save session settings

## Troubleshooting

### Issue: "Module not found" error
**Solution**: Make sure you've activated the virtual environment and installed all requirements.

### Issue: Plot not displaying
**Solution**: Ensure your browser allows JavaScript and try refreshing the page.

### Issue: File upload fails
**Solution**: Check that your file is in CSV or TSV format and contains the required columns.

### Issue: Kaleido error when downloading plots
**Solution**: Reinstall kaleido: `pip install kaleido --force-reinstall`

## Dependencies

- streamlit: Web application framework
- pandas: Data manipulation and analysis
- plotly: Interactive plotting library
- kaleido: Static image export for plotly

## License

This application is provided as-is for academic and research purposes.

## Support

For issues or questions, please check:
1. This README file
2. Streamlit documentation: https://docs.streamlit.io
3. Plotly documentation: https://plotly.com/python/

## Citation

If you use this tool in your research, please cite the relevant packages:
- DESeq2 (Love et al., 2014)
- READemption pipeline
- Plotly
- Streamlit

---

**Version**: 1.0.0  
**Last Updated**: 2024
