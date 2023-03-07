# BlastUI
Web User Interface for blast

## Requirements
- python 3.10
- streamlit
- streamlit-aggrid 0.3.3 right now 0.3.4post2 doesn't not work well
- streamlit-extras
- streamlit-option-menu
- pandas
- bokeh 2.4.3 no dependencies

## Installation
If you have conda
```
git pull https://github.com/Alessandro201/BlastUI/
conda create -n BlastUI python=3.10
conda activate BlastUI
pip install pandas streamlit streamlit-aggrid==0.3.3 streamlit-extras streamlit-option-menu
pip install --no-deps bokeh==2.4.3

```

## Quickstart
```
python src/main.py
```
