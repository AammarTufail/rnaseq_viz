#!/bin/bash

echo "Setting up RNA-seq Volcano Plot Analyzer..."

# Create virtual environment
python3 -m venv .venv

# Activate virtual environment
source .venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install requirements
# pip install -r requirements.txt
pip install streamlit plotly kaleido

echo "Setup complete! Activate the environment with: source venv/bin/activate"
echo "Then run the app with: streamlit run app.py"
