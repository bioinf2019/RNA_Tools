# RNA_Tools

## Introduction

In this repository you will find the python script of software called "All_In_One_shRNA_Generator".

## Technical details

* The software was developed using Python v3.5
* The software GUI was developed using PyQt5

## Requirenments before run the script
### **Windows**

Before to run the script in Windows you need to install:

* [Python 3.5.x](/download/python-3.5.0.exe)
* [PyQt5](/download/PyQt5-5.6-gpl-Py3.5-Qt5.6.0-x32-2.exe)

#### Installing Python3.5

These are the files that we have after download the previous files

![](../images/installing_python3.5.PNG)

Double click on the Python installer and follow the steps

![](/images/installing_python3.5_00.PNG)

Step 1:

![](/images/installing_python3.5_01.PNG)

Step 2:

![](/images/installing_python3.5_02.PNG)

Step 3:

![](/images/installing_python3.5_03.PNG)

#### Installing PyQt5

Double click on the PyQt5 installer and follow the steps

![](/images/installing_pyQT5.PNG)

Step 1:

![](/images/installing_pyQT5_00.PNG)

Step 2:

![](/images/installing_pyQT5_01.PNG)

Step 3:

![](/images/installing_pyQT5_02.PNG)

Step 4:

![](/images/installing_pyQT5_03.PNG)

Step 5:

![](/images/installing_pyQT5_04.PNG)

#### Running the script on windows

The final step for running the script is open the python script directory

![](/images/runnig_the_script_00.PNG)

Right click on the path and select copy

![](/images/runnig_the_script_01.PNG)

Click on the start button

![](/images/runnig_the_script_02.PNG)

Type the word **cmd** and press **ENTER**

![](/images/runnig_the_script_03.PNG)

Type **cd** and paste the copied path and press **ENTER**

![](/images/runnig_the_script_04.PNG)

Then type the word **si** that are the inital letters od the script that we want to run, and press **TAB** and then press **ENTER**

![](/images/runnig_the_script_05.PNG)

Here we have the software running

![](/images/runnig_the_script_06.PNG)

### **Linux (Ubuntu)**

Before to run the script in Ubuntu Linux we need to install:
* [Installing PIP](https://www.rosehosting.com/blog/how-to-install-pip-on-ubuntu-16-04/)
Update and upgrade the system
````
sudo apt-get update && sudo apt-get -y upgrade
````
Next install pip3 
````
sudo apt-get install python-pip3
````
Finally verify the installation

````
pip3 -V
````
* [Installing PyQt5](https://pypi.org/project/PyQt5/)
The version that its needed is for pip3

````
pip3 install PyQt5
````
#### Run the script on Ubuntu

To run the script open the directory where is placed your python file on the Terminal.

````
python3 All_In_One_shRNA_Generator_release_v1.0.py
````