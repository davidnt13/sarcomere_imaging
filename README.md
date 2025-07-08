# sarcomere_imaging
Sarcomere imaging analysis

# Installation
1. Create an environment and activate it.
If using conda:
```
conda create -n sarc-env python=3.11
conda activate sarc-env
```
If using venv and on Mac/Linux:
```
python -m venv sarc-venv
source sarc-venv/bin/activate
```
If using venv and on Windows:
```
python -m venv sarc-venv
sarc-venv\Scripts\activate
```
In VScode you can select the environment by using Ctrl+Shift+P, type Python Select Interpreter and selecting ```./venv/bin/python```.

2. Clone the Repository
```
git clone https://github.com/davidnt13/sarcomere_imaging.git
cd sarcomere_imaging
```

3. Install necessary packages and modules,
```
python -m pip install -e .
```
And If Git is installed:
```
python -m pip install -r requirements.txt
```
But If Git is not Installed:
i. Go to https://github.com/davidnt13/SarcGraph
ii. Click the green <> code button and select Download Zip
iii. cd to the downloaded SarcGraph directory and install in editable mode:
```
pip install -e .
```


# Image Analysis Pipeline
Once in the sarcomere_imaging directory, the main step is running the run_pipeline.py file.
- This will present a popup window asking for content such as the data (tif images) directory, the date (for organization of files), and whether plots are desired.
- Input the necessary selections, hit submit, and the code will run, producing the following:
    - A folder for each image containing:
        - Sarcomere channel of the image
        - The filtered image sarcomere channel
        - Sarcomeres.csv: SarcGraph output of the file (Used for data analysis)
        - sarcomere_mask.tif: SarcAsM's sarcomere mask of the image
        - zbands.tif: SarcAsM's z-band mask
        - {Image Name}_sg_sarcasm_output: Txt file featuring the SarcAsM's analysis for comparison purposes
    - image_analysis_{stat_title}.csv: Summary Statistics of the Image Data
    - image_analysis_{stat_title}.png: Plot of the data if desired
