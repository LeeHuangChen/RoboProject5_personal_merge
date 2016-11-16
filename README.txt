
Please follow the steps below to install MMTK, a toolkit we are using for this project.
##########################
#   INSTALLATION STEPS   #
##########################

ADD THE DOWNLOAD DIRECTORY FOR MMTK
  (Referance: http://installion.co.uk/ubuntu/vivid/universe/p/python-mmtk/install/index.html)
- goto terminal
- sudo gedit /etc/apt/sources.list
- add following line to the file if it doesn't have it already
  - deb http://us.archive.ubuntu.com/ubuntu vivid main universe
- terminal: type
  - sudo apt-get update

DOWNGRADE NUMPY TO VER:1.8
- pip install numpy==1.8
- go to python by typing 'python'
- type the following to make sure the version is now 1.8.0
  >>> import numpy
  >>> numpy.__version__

INSTALL SCIENTIFIC PYTHON:
- download scientific python 2.9.4 with this link:
  - https://sourcesup.renater.fr/frs/download.php/file/4570/ScientificPython-2.9.4.tar.gz
- extract it by 
  - tar xzf ScientificPython-2.9.4.tar.gz
- install
  - cd ScientificPython-2.9.4
  - python setup.py build
  - sudo python setup.py install

INSTAL MMTK:
- sudo apt-get install python-mmtk

TO TEST:
- goto python
  >>> from MMTK import *
- the line should not return any errors.