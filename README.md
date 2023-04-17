# TAO23

Tutorials for Techniques in Astronomical Observation, 2023A semester, Seoul National University.

## Tutorials
* Tutorial 1: [PSF photometry](https://nbviewer.org/github/hbahk/TAO23/blob/main/tutorials/PSF_photometry.ipynb)

## 1. Preparation

### 1.1. Clone this repository

```shell
cd <PathToGitClone>
git clone https://github.com/hbahk/TAO23.git
```



### 1.2. Software installations

#### 1.2.1 Python packages

* Useful links
  
  [SNU_AOclass/00-3_Prepare_Python.md at master · ysBach/SNU_AOclass · GitHub](https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/00-3_Prepare_Python.md)
  
  [SNU_AOclass/00-1_Softwares.md at master · ysBach/SNU_AOclass · GitHub](https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/00-1_Softwares.md)
  
  
1. Create a new conda environment (optional)

```shell
conda create -n tao23 python numpy scipy pandas jupyter
conda activate tao23
```

2. Install packages via conda

```shell
conda install -c conda-forge astropy sep tqdm
conda install -c conda-forge photutils scikit-learn-intelex
conda install -c astropy ccdproc astroscrappy specutils
```

3. Install package via Github repositories

```shell
cd <PathToGitClone> 
git clone https://github.com/jrjohansson/version_information.git && git clone https://github.com/astropy/astroquery.git && git clone https://github.com/ejeschke/ginga.git && git clone https://github.com/quatrope/astroalign
```

```shell
cd version_information && pip install -e . && cd ..
cd astroquery && git pull && pip install -e . && cd ..
cd ginga && git pull && pip install -e . && cd .. 
cd astroalign && git pull && pip install -e . && cd .. 
```

4. Set the python path of your own kernels/editors with this environments

cf) For Spyder users, install `spyder-kernels` for your version of Spyder following the table in this link [Common Illnesses &#8212; Spyder 5 documentation](https://docs.spyder-ide.org/current/troubleshooting/common-illnesses.html#spyder-kernels-not-installed-incompatible). In my case, the version of my Spyder is 5.2.2, so I had to install `spyder-kernels=2.2.1`.

```shell
conda install spyder‑kernels=2.2.1
```

#### 1.2.2 Source Extractor

* **For Linux (or Window) users**, please refer to the official documentation:
  [Installing the software](https://sextractor.readthedocs.io/en/latest/Installing.html)
  For Ubuntu users, for example, simple command will do:
  ```shell
  sudo apt-get install sextractor
  ```


* **For Mac users**, using Homebrew is easiest way, I suppose.
  ```shell
  brew install sextractor
  ```
  If you don't have Homebrew intalled yet, then [give that a try!](https://brew.sh)