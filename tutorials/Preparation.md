# Preparation

## 1. Clone the repository

```shell
cd <PathToGitClone>
git clone https://github.com/hbahk/TAOTutorials.git
```



## 2. Software installations

### 2.1 Python packages

* Useful links
  
  [SNU_AOclass/00-3_Prepare_Python.md at master · ysBach/SNU_AOclass · GitHub](https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/00-3_Prepare_Python.md)
  
  [SNU_AOclass/00-1_Softwares.md at master · ysBach/SNU_AOclass · GitHub](https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/00-1_Softwares.md)
  
  
1. Create a new conda environment (optional)

```shell
conda create -n tao python numpy scipy pandas jupyter
conda activate tao
```

2. Install packages via conda

```shell
conda install -c conda-forge astropy sep tqdm
conda install -c conda-forge photutils
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

:::{admonition} For Spyder users
:class: dropdown tip

For Spyder users, install `spyder-kernels` for your version of Spyder following the table in this link [Common Illnesses &#8212; Spyder 5 documentation](https://docs.spyder-ide.org/current/troubleshooting/common-illnesses.html#spyder-kernels-not-installed-incompatible). In my case, the version of my Spyder is 5.2.2, so I had to install `spyder-kernels=2.2.1`.


```shell
conda install spyder‑kernels=2.2.1
```
:::

### 2.2 Source Extractor

* **For Linux (or Window) users**, please refer to the official documentation: [Installing the software](https://sextractor.readthedocs.io/en/latest/Installing.html).

  For Ubuntu users, for example, simple command will do:
  ```shell
  sudo apt-get install sextractor
  ```



* **For Mac users**, using Homebrew is easiest way, I suppose.
  ```shell
  brew install sextractor
  ```
  If you don't have Homebrew intalled yet, then [give that a try!](https://brew.sh)

* You can check whether SExtractor is installed or not with the command of:
  ```shell
  sex
  ```
  Or if you installed with the `apt-get` (and other possible distributors) then try
  ```shell
  sextractor
  ```
  or
  ```shell
  source-extractor
  ```

### 2.3. PSF Extractor

* **For Linux (or Window) users**, please refer to the official documentation: [Installing the software](https://psfex.readthedocs.io/en/latest/Installing.html).
  For Ubuntu users, for example, executing simple command will install this software:
  ```shell
  sudo apt-get install psfex
  ```

* **For Mac users** (Tested with Apple M1 Pro Ventura 13.2.1),
  ```shell
  cd <PathToGitClone>
  git clone https://github.com/astromatic/psfex.git
  ```

  Check whether required packages are installed: `autoconf`, `automake`, `libtool`.
  ```shell
  brew info autoconf
  brew info automake
  brew info libtool
  ```
  If you are missing one of them, then install with the command `brew install <MissingPackage>`.

  Then 
  ```shell
  cd psfex
  ./autogen.sh
  ```
  Check the result and find the pathes of `lib` and `include` directory of the `OpenBLAS` package. You can find their location with command of `brew info openblas`, if you installed `OpenBLAS` with Homebrew, which should be the case when you installed SExtractor with Homebrew. Similarly, find `lib` and `include` for `FFTW`, and then try
  ```shell
  ./configure --disable-silent-rules --enable-openblas --with-openblas-libdir="/opt/homebrew/opt/openblas/lib" --with-openblas-incdir="/opt/homebrew/opt/openblas/include" --with-fftw-libdir="/opt/homebrew/opt/fftw/lib" --with-fftw-incdir="/opt/homebrew/opt/fftw/include"
  ```
  Here you should put the `libdir` and `incdir` with your own pathes.

  Finally,
  ```shell
  sudo make install
  ```
  and now you are all set.




