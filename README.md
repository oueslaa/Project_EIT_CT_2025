## Requirements

You need Python 3.9 or later.I use VS code but you can use a virtual environment too.
## Install Python
You can install the latest Python version here : https://www.python.org/downloads/. During the installation, it is VERY IMPORTANT to Check the box “Add Python to PATH” at the bottom, then click “Install".


In VS code: Go to Extensions in the menu at left, and search for `Python` than install  the official extension.

It is posssbile that you need to install pip as well.
Open a terminal in VS Code:  Go to `Terminal > New Terminal`
    Type :
     `pip --version`
    If it shows a version number, then pip is already installed.
    If not type :      `python -m ensurepip --upgrade` and it should install pip
    Now : 
        `pip --version`    should    show a version number.

### Install dependencies:
To install the libraries to run the code, type in a VS code terminal
`pip install -r requirements.txt`

Make sure your path is `C:\Users\'%USERNAME%'\Dropbox\Project_CT` where the file requirements.txt is.


Once the libraries installed,you can choose a script or open the folder in WORKSPACE.Open a script and  with the `Start` button at the the top right corner -> Run Python file.
You can also juste type in the terminal : `python '.\src_aatae\Liver V and S vs age auto.py'` for example.


the code is still repetitive and redundant for now, so if you would like to test it for differant Organs/subjects, you can adapt it manually.(At the start of each script where there is a comment  #To adapt)


## Interactive visualization with Napari

To inspect the CT scans and  segmented volumes interactively(just like visualize_nii.m in MATLAB), you can use the Napari viewer:

`pip install napari[all]`

The installation will probably fail with an error like:(if not great)

`error: Microsoft Visual C++ 14.0 or greater is required`

You need to install Microsoft C++ Build Tools --> https://visualstudio.microsoft.com/visual-cpp-build-tools/

Only the script "3D with Napari.py" uses Napari, for now.
