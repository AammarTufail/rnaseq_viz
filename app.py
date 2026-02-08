import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import re
from io import BytesIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Page configuration with theme toggle
st.set_page_config(
    page_title="RNA-seq DESeq Volcano Plot Analyzer", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Theme toggle in sidebar
with st.sidebar:
    st.header("⚙️ Settings")
    theme_mode = st.toggle("🌙 Dark Mode", value=False)
    
    # Apply custom CSS based on theme
    if theme_mode:
        st.markdown("""
        <style>
        .stApp {
            background-color: #0e1117;
            color: #fafafa;
        }
        </style>
        """, unsafe_allow_html=True)

# Helper functions
def extract_gene_info(attributes_str):
    """Extract gene name and locus tag from Attributes column"""
    gene_name = "Unknown"
    locus_tag = "Unknown"
    
    if pd.isna(attributes_str):
        return gene_name, locus_tag
    
    gene_match = re.search(r'gene=([^;]+)', str(attributes_str))
    if gene_match:
        gene_name = gene_match.group(1)
    
    locus_match = re.search(r'locus_tag=([^;]+)', str(attributes_str))
    if locus_match:
        locus_tag = locus_match.group(1)
    
    if gene_name == "Unknown":
        name_match = re.search(r'Name=([^;]+)', str(attributes_str))
        if name_match:
            gene_name = name_match.group(1)
        else:
            product_match = re.search(r'product=([^;]+)', str(attributes_str))
            if product_match:
                gene_name = product_match.group(1)[:30]
    
    return gene_name, locus_tag

def filter_deseq_data(df, padj_threshold, log2fc_up, log2fc_down):
    """Filter and categorize DESeq data based on thresholds"""
    df_copy = df.copy()
    
    conditions = [
        (df_copy['padj'] < padj_threshold) & (df_copy['log2FoldChange'] >= log2fc_up),
        (df_copy['padj'] < padj_threshold) & (df_copy['log2FoldChange'] <= log2fc_down)
    ]
    choices = ['Upregulated', 'Downregulated']
    
    df_copy['Significance'] = np.select(conditions, choices, default='Not Significant')
    df_copy['-log10(padj)'] = -np.log10(df_copy['padj'].replace(0, 1e-300))
    
    return df_copy

def create_volcano_plot(df, color_up, color_down, color_ns, padj_threshold, log2fc_up, log2fc_down, width, height):
    """Create an interactive volcano plot"""
    
    df['hover_text'] = (
        '<b>' + df['gene_name'] + '</b><br>' +
        'Locus: ' + df['locus_tag'] + '<br>' +
        'log2FC: ' + df['log2FoldChange'].round(3).astype(str) + '<br>' +
        'padj: ' + df['padj'].apply(lambda x: f'{x:.2e}') + '<br>' +
        'Significance: ' + df['Significance']
    )
    
    fig = px.scatter(
        df,
        x='log2FoldChange',
        y='-log10(padj)',
        color='Significance',
        color_discrete_map={
            'Upregulated': color_up,
            'Downregulated': color_down,
            'Not Significant': color_ns
        },
        hover_data={'hover_text': True, 'log2FoldChange': False, '-log10(padj)': False, 'Significance': False},
        title='Volcano Plot: Differential Gene Expression'
    )
    
    fig.update_traces(hovertemplate='%{customdata[0]}<extra></extra>')
    
    fig.add_hline(
        y=-np.log10(padj_threshold),
        line_dash="dash",
        line_color="black",
        opacity=0.5,
        annotation_text=f"padj = {padj_threshold}",
        annotation_position="right"
    )
    
    fig.add_vline(x=log2fc_up, line_dash="dash", line_color="black", opacity=0.5)
    fig.add_vline(x=log2fc_down, line_dash="dash", line_color="black", opacity=0.5)
    
    fig.update_layout(
        xaxis_title='log2(Fold Change)',
        yaxis_title='-log10(adjusted p-value)',
        height=height,
        width=width,
        template='plotly_white',
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    return fig

@st.cache_data
def load_data(uploaded_file):
    """Load and cache the uploaded data"""
    df = pd.read_csv(
        uploaded_file, 
        sep='\t',
        comment='#',
        na_values=['NA', 'NaN', '']
    )
    return df

st.title("🧬 RNA-seq DESeq Volcano Plot Analyzer")
st.markdown("---")

# Sidebar
with st.sidebar:
    uploaded_file = st.file_uploader("Upload DESeq Output File", type=['csv', 'tsv', 'txt'])
    
    if uploaded_file:
        st.success("File uploaded successfully!")
    
    st.markdown("---")
    st.subheader("📊 Filter Criteria")
    padj_threshold = st.number_input("Adjusted p-value threshold", 
                                      min_value=0.0, max_value=1.0, 
                                      value=0.05, step=0.01, format="%.3f")
    
    log2fc_up = st.number_input("log2FC threshold (Upregulated)", 
                                 min_value=0.0, value=1.0, step=0.1)
    
    log2fc_down = st.number_input("log2FC threshold (Downregulated)", 
                                   max_value=0.0, value=-1.0, step=0.1)
    
    st.markdown("---")
    st.subheader("🎨 Color Settings")
    color_up = st.color_picker("Upregulated", "#FF0000")
    color_down = st.color_picker("Downregulated", "#0000FF")
    color_ns = st.color_picker("Not Significant", "#808080")
    
    st.markdown("---")
    st.subheader("📐 Plot Dimensions")
    plot_width = st.slider("Plot Width (px)", 600, 1600, 1200, 50)
    plot_height = st.slider("Plot Height (px)", 400, 1200, 700, 50)
    
    st.markdown("---")
    st.subheader("💾 Export Settings")
    dpi_value = st.selectbox("DPI for Image Export", [150, 300, 600, 900], index=1)

# Main content
if uploaded_file is not None:
    try:
        # Load data with caching
        df = load_data(uploaded_file)
        st.success(f"✅ File loaded successfully! {len(df)} genes found.")
        
        # Extract gene information
        if 'Attributes' in df.columns:
            df[['gene_name', 'locus_tag']] = df['Attributes'].apply(
                lambda x: pd.Series(extract_gene_info(x)))
        else:
            st.warning("'Attributes' column not found. Using row indices for labels.")
            df['gene_name'] = df.index.astype(str)
            df['locus_tag'] = df.index.astype(str)
        
        # Remove rows with NA in critical columns
        df = df.dropna(subset=['log2FoldChange', 'padj'])
        
        # Filter data
        df_filtered = filter_deseq_data(df, padj_threshold, log2fc_up, log2fc_down)
        
        # Display statistics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Genes", len(df_filtered))
        with col2:
            st.metric("Upregulated", len(df_filtered[df_filtered['Significance'] == 'Upregulated']))
        with col3:
            st.metric("Downregulated", len(df_filtered[df_filtered['Significance'] == 'Downregulated']))
        with col4:
            st.metric("Not Significant", len(df_filtered[df_filtered['Significance'] == 'Not Significant']))
        
        st.markdown("---")
        
        # Tabs for different visualizations
        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
            "🌋 Volcano Plot", "📊 MA Plot", "📈 Box Plot", 
            "🔥 Heatmap", "📉 P-value Distribution", "📋 Data Table"
        ])
        
        with tab1:
            st.subheader("Volcano Plot")
            fig_volcano = create_volcano_plot(df_filtered, color_up, color_down, color_ns, 
                                            padj_threshold, log2fc_up, log2fc_down, 
                                            plot_width, plot_height)
            st.plotly_chart(fig_volcano, use_container_width=True)
            
            # Download options for volcano plot
            col1, col2, col3 = st.columns(3)
            with col1:
                img_png = fig_volcano.to_image(format="png", width=plot_width, height=plot_height, scale=dpi_value/100)
                st.download_button("📥 Download PNG", data=img_png, 
                                 file_name=f"volcano_plot_{dpi_value}dpi.png", 
                                 mime="image/png")
            with col2:
                img_svg = fig_volcano.to_image(format="svg", width=plot_width, height=plot_height)
                st.download_button("📥 Download SVG", data=img_svg, 
                                 file_name="volcano_plot.svg", 
                                 mime="image/svg+xml")
            with col3:
                html_str = fig_volcano.to_html()
                st.download_button("📥 Download HTML", data=html_str, 
                                 file_name="volcano_plot.html", 
                                 mime="text/html")
        
        with tab2:
            st.subheader("MA Plot")
            if 'baseMean' in df_filtered.columns:
                df_filtered['log10(baseMean)'] = np.log10(df_filtered['baseMean'].replace(0, 1e-10))
                
                fig_ma = px.scatter(df_filtered,
                                  x='log10(baseMean)',
                                  y='log2FoldChange',
                                  color='Significance',
                                  color_discrete_map={
                                      'Upregulated': color_up,
                                      'Downregulated': color_down,
                                      'Not Significant': color_ns
                                  },
                                  title="MA Plot: Expression vs Fold Change",
                                  width=plot_width,
                                  height=plot_height)
                
                st.plotly_chart(fig_ma, use_container_width=True)
                
                # Download MA plot
                img_png = fig_ma.to_image(format="png", width=plot_width, height=plot_height, scale=dpi_value/100)
                st.download_button("📥 Download MA Plot (PNG)", data=img_png, 
                                 file_name=f"ma_plot_{dpi_value}dpi.png", 
                                 mime="image/png")
        
        with tab3:
            st.subheader("Box Plot: Log2 Fold Change Distribution")
            
            fig_box = go.Figure()
            for sig_type in ['Upregulated', 'Not Significant', 'Downregulated']:
                data = df_filtered[df_filtered['Significance'] == sig_type]['log2FoldChange']
                fig_box.add_trace(go.Box(
                    y=data,
                    name=sig_type,
                    marker_color=color_up if sig_type == 'Upregulated' 
                               else (color_down if sig_type == 'Downregulated' else color_ns)
                ))
            
            fig_box.update_layout(
                title="Distribution of Log2 Fold Changes",
                yaxis_title="log2(Fold Change)",
                width=plot_width,
                height=plot_height
            )
            st.plotly_chart(fig_box, use_container_width=True)
            
            img_png = fig_box.to_image(format="png", width=plot_width, height=plot_height, scale=dpi_value/100)
            st.download_button("📥 Download Box Plot (PNG)", data=img_png, 
                             file_name=f"boxplot_{dpi_value}dpi.png", 
                             mime="image/png")
        
        with tab4:
            st.subheader("Heatmap: Top Differentially Expressed Genes")
            
            # Get count columns
            count_cols = [col for col in df_filtered.columns if 'countings' in col.lower()]
            
            if len(count_cols) > 0:
                # Select top genes by padj
                n_genes = st.slider("Number of genes to display", 10, 50, 20)
                top_genes = df_filtered.nsmallest(n_genes, 'padj')
                
                # Prepare data for heatmap
                count_data = top_genes[count_cols + ['gene_name']].set_index('gene_name')
                
                # Create matplotlib figure for better heatmap
                fig, ax = plt.subplots(figsize=(plot_width/100, plot_height/100))
                sns.heatmap(count_data.T, cmap='RdYlBu_r', center=0, 
                           robust=True, ax=ax, cbar_kws={'label': 'Normalized Counts'})
                ax.set_xlabel('Gene')
                ax.set_ylabel('Sample')
                ax.set_title(f'Top {n_genes} Differentially Expressed Genes')
                plt.tight_layout()
                st.pyplot(fig)
                
                # Save heatmap
                buf = BytesIO()
                fig.savefig(buf, format='png', dpi=dpi_value, bbox_inches='tight')
                buf.seek(0)
                st.download_button("📥 Download Heatmap (PNG)", data=buf, 
                                 file_name=f"heatmap_{dpi_value}dpi.png", 
                                 mime="image/png")
            else:
                st.warning("No count columns found in the dataset.")
        
        with tab5:
            st.subheader("P-value Distribution")
            
            fig_pval = go.Figure()
            fig_pval.add_trace(go.Histogram(
                x=df_filtered['padj'],
                nbinsx=50,
                name='Adjusted P-value',
                marker_color='steelblue'
            ))
            
            fig_pval.update_layout(
                title="Distribution of Adjusted P-values",
                xaxis_title="Adjusted P-value",
                yaxis_title="Frequency",
                width=plot_width,
                height=plot_height
            )
            st.plotly_chart(fig_pval, use_container_width=True)
            
            # QQ plot for p-values
            st.subheader("Q-Q Plot for P-values")
            fig_qq, ax_qq = plt.subplots(figsize=(plot_width/100, plot_height/100))
            
            # Remove NA values
            pvalues = df_filtered['padj'].dropna()
            stats.probplot(pvalues, dist="uniform", plot=ax_qq)
            ax_qq.set_title("Q-Q Plot of Adjusted P-values")
            ax_qq.set_xlabel("Theoretical Quantiles")
            ax_qq.set_ylabel("Sample Quantiles")
            st.pyplot(fig_qq)
            
            img_png = fig_pval.to_image(format="png", width=plot_width, height=plot_height, scale=dpi_value/100)
            st.download_button("📥 Download P-value Distribution (PNG)", data=img_png, 
                             file_name=f"pvalue_dist_{dpi_value}dpi.png", 
                             mime="image/png")
        
        with tab6:
            st.subheader("Filtered Data Table")
            
            # Filter options
            col1, col2 = st.columns(2)
            with col1:
                show_sig = st.selectbox("Show", ["All", "Upregulated", "Downregulated", "Not Significant"])
            with col2:
                search_term = st.text_input("Search genes", "")
            
            # Apply filters
            display_df = df_filtered.copy()
            if show_sig != "All":
                display_df = display_df[display_df['Significance'] == show_sig]
            if search_term:
                display_df = display_df[display_df['gene_name'].str.contains(search_term, case=False, na=False)]
            
            st.dataframe(display_df, use_container_width=True, height=400)
            
            # Download filtered data
            csv = display_df.to_csv(index=False)
            st.download_button("📥 Download Filtered Data (CSV)", data=csv, 
                             file_name="filtered_deseq_results.csv", 
                             mime="text/csv")
        
    except Exception as e:
        st.error(f"Error processing file: {str(e)}")
        st.exception(e)

else:
    st.info("👆 Please upload a DESeq output file to begin analysis.")
    
    with st.expander("📋 Expected File Format"):
        st.markdown("""
        Your file should:
        - Be tab-separated (TSV format)
        - Have comment lines starting with `#` (will be skipped)
        - Contain these key columns:
          - `log2FoldChange`
          - `padj` (adjusted p-value)
          - `baseMean`
          - Count columns (ending with 'countings')
          - `Attributes` (optional, for gene names)
        """)

st.markdown("---")
st.markdown("**Created with ❤️ for RNA-seq Analysis | Powered by Streamlit & Plotly**")
