# paleogeo_functions

This is a repository for helper function library based on pyGPlates python library and matplotlib for static paleogeography reconstructions. 

## Installation

Currently this package is not registered on PyPI. You can install it by adding this GitHub repository directly to your Python environment. 

### for Mac OS/Linux

Determine whether your Python environment uses .bash or .zsh as the shell. Download this GitHub repository to your local machine, for example to `/Users/username/Github/paleogeo_functions/`.

Add the following line to your `~/.bash_profile` or `~/.zshrc` file:

```
export PYTHONPATH="$PYTHONPATH:/Users/username/Github/paleogeo_functions/"
```
replacing `username` with your actual username.

Then, run the following command to apply the changes:

```
source ~/.bash_profile
```
or
```
source ~/.zshrc
```

## Usage
first install pyGPlates using conda (https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#install-using-conda):


```
conda install -c conda-forge pygplates
```

You can import the functions in your Python scripts or Jupyter notebooks as follows:

```
python
from paleogeo_functions import *
```
