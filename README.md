# BlastUI

A graphical interface for BLAST to run queries against your own local sequences.

From the [blast website](https://blast.ncbi.nlm.nih.gov/Blast.cgi):
> The Basic Local Alignment Search Tool (BLAST) finds regions of local similarity between sequences.
> The program compares nucleotide or protein sequences to sequence databases and calculates the
> statistical significance of matches. BLAST can be used to infer functional and evolutionary relationships
> between sequences as well as help identify members of gene families.

## Features

- Blast queries against your own sequences
- Visualize the results, filter and aggregate them
- Graphical summary of the alignments
- Save the results as Excel or CSV files
- Add your own scripts to analyze the data
- Save the hits sequences or the alignments
- Find multiple hits from the same sequence

## Next features

- [ ] Auto update the app when a new version is available
- [ ] Right now you can run the program on a server and all the users in the network can access it and perform analysis. I plan on implementing the ability to let users log in, have their own sequences and queue queries so that many users can use the software at the same time

<br>
<br>



https://github.com/Alessandro201/BlastUI/assets/61567683/6b409c6d-38cb-428f-86cd-4c8f79f12485



![immagine](https://user-images.githubusercontent.com/61567683/227249073-3cb94f8e-e045-40be-8ff9-91de799537bb.png)

![immagine](https://user-images.githubusercontent.com/61567683/227252687-d1fb102a-72c4-47b4-91eb-17f617ef9a5e.png)

![immagine](https://user-images.githubusercontent.com/61567683/227253947-c1a8f3ec-d255-406b-848f-33985cc26c14.png)

![immagine](https://user-images.githubusercontent.com/61567683/227254938-732ed1ac-27a5-4f04-a49e-186d47fb180c.png)

> The data used to produce the images is not the same in all of them.

## Installation

#### Conda

```
git clone https://github.com/Alessandro201/BlastUI/
conda env create -n BlastUI --file environment.yaml
conda activate BlastUI
```

Launch the app with
```
python ./run_app.py
```

## Contributing 

If you want to package the app into an executable, you need to install pyinstaller to package it:

```
pip install pyinstaller
```

And once you're done, you can package the app with:

```
pyinstaller --clean -y .\create_package.spec 
```

Check that `./.streamlit`, `./pages` and `./scripts` folders are present in the top-level directory of the output
otherwise the application won't work.
