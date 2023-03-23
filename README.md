# BlastUI
From the [blast website](https://blast.ncbi.nlm.nih.gov/Blast.cgi): 
> The Basic Local Alignment Search Tool (BLAST) finds regions of local similarity between sequences. 
The program compares nucleotide or protein sequences to sequence databases and calculates the 
statistical significance of matches. BLAST can be used to infer functional and evolutionary relationships 
between sequences as well as help identify members of gene families. 

This WebApp is a simple interface to run blast locally. 

## Features
- Blast queries against your own sequences
- Visualize the results, filter and aggregate them
- Graphical summary of the alignments
- Save the results as Excel or CSV files
- Save the hit sequences or the alignments
- Find multiple hits from the same sequence

It also has a dark mode.

<br>
<br>


![immagine](https://user-images.githubusercontent.com/61567683/227249073-3cb94f8e-e045-40be-8ff9-91de799537bb.png)

![immagine](https://user-images.githubusercontent.com/61567683/227252687-d1fb102a-72c4-47b4-91eb-17f617ef9a5e.png)

![immagine](https://user-images.githubusercontent.com/61567683/227253947-c1a8f3ec-d255-406b-848f-33985cc26c14.png)

![immagine](https://user-images.githubusercontent.com/61567683/227254938-732ed1ac-27a5-4f04-a49e-186d47fb180c.png)


The data used to produce the images is not the same in all of them. 

## Requirements
- python >= 3.10
- streamlit
- streamlit-aggrid 0.3.3 right now 0.3.4post2 does not work well
- streamlit-extras
- streamlit-option-menu
- pandas
- bokeh 2.4.3 no dependencies

## Installation
#### Conda

```
git pull https://github.com/Alessandro201/BlastUI/
conda create -n BlastUI python=3.10
conda activate BlastUI
pip install pandas streamlit streamlit-aggrid==0.3.3 streamlit-extras streamlit-option-menu
pip install --no-deps bokeh==2.4.3

```

## Quickstart
```
python -m streamlit src/Home.py
```
