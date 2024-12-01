@echo off
python -m venv myenv
call myenv\Scripts\activate
pip install -r requirements.txt
python main.py
pause