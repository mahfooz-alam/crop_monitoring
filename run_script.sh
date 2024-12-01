#!/bin/bash

echo "Running Python script..."
python3 -m venv myenv
source myenv/bin/activate
pip install -r requirements.txt
python3 main.py
read -p "Press Enter to continue..."
