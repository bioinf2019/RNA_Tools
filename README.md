# RNA_Tools

## Introduction

In this repository you will find the python script of software called "All_In_One_shRNA_Generator".

## Technical details

* The software was developed using Python v3.5
* The software GUI was developed using PyQt5

## Requirenments before run the script
### **Windows** [Instructive video](https://youtu.be/jG7qKrKMu8M)

Before to run the script in Windows you need to download:

1. [Python version 3.5.0](https://www.python.org/ftp/python/3.5.0/python-3.5.0.exe)
2. [PyQt5 version 5.6.0](https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.6/PyQt5-5.6-gpl-Py3.5-Qt5.6.0-x32-2.exe/download)
3. [SSD Software script](/download/siRNA_and_shRNA_designer_(SSD)_release_v1.0.zip)

The orther of the installations must be follows as shown.

#### Installing Python3.5

These are the files that we have after download the previous files

![](/images/installing_python3.5.png)

Double click on the Python installer and follow the steps

![](/images/installing_python3.5_00.png)

Step 1:

![](/images/installing_python3.5_01.png)

Step 2:

![](/images/installing_python3.5_02.png)

Step 3:

![](/images/installing_python3.5_03.png)

#### Installing PyQt5

Double click on the PyQt5 installer and follow the steps

![](/images/installing_pyQT5.png)

Step 1:

![](/images/installing_pyQT5_00.png)

Step 2:

![](/images/installing_pyQT5_01.png)

Step 3:

![](/images/installing_pyQT5_02.png)

Step 4:

![](/images/installing_pyQT5_03.png)

We need to verify the path where the PyQt5 will be installed its looks like as follow:

![](/images/pyqt5_path.png)

Step 5:

![](/images/installing_pyQT5_04.png)

#### Running the script on windows

Once you have downloaded the zipped file that is inside the Python script you will see like this:

![](/images/runnig_the_script_000.png)

First **right click** onthe zipped file and press the **Extract Here** option.

![](/images/runnig_the_script_001.png)

Then you will have the Python script unzipped. To runn it you just need to **Double click** on the script and it will run.

![](/images/runnig_the_script_002.png)

Finally you will appear the window of the software as shown:

![](/images/runnig_the_script_06.png)


##### NOTES:
If you cannot run with double cick you have to change the extension file from **py** to **pyw** or just use the command prompt to run it. 

[Running the SSD software using the command prompt](https://youtu.be/X0S5jYU3vnU)

### **Linux (Ubuntu)** [Instructive video](https://www.youtube.com/embed/FC1ttM7NY)

Before to run the script in Windows you need to download:

1. PIP3
2. PyQt5
3. [SSD Software script](/download/siRNA_and_shRNA_designer_(SSD)_release_v1.0.zip)

First you need to checke the version of python3 you have using the command:
````
python3 --version
````

![](/images/installation_pip3.png)

You have to see the version you have installed like this:

![](/images/installation_pip3_00.png)

Next step is to update and upgrade the system in order to install pip3:

````
sudo apt-get update
````

![](/images/installation_pip3_04.png)

````
sudo apt-get upgrade
````

![](/images/installation_pip3_05.png)

The system will show you how to install pip3 just runnin this command:
````
pip3
````

![](/images/installation_pip3_01.png)

You have to copy the command that shows you.

![](/images/installation_pip3_02.png)

Then paste on the terminal and run it.

![](/images/installation_pip3_03.png)

You can check the pip3 version using:
````
pip3 -V
````

![](/images/installation_pip3_06.png)

After that we need to install PyQt5, so it will be done with command:
````
pip3 install pyqt5
````

![](/images/installation_pip3_07.png)

#### Running the script on Linux (Ubuntu)

Open the directory where is you script target and the right click and select **Open on terminal**.
Then type the following command:
````
python3 siRNA\ and\ shRNA\ designer\ (SSD\)_release_v1.0.py
````

![](/images/installation_pip3_08.png)

And finally will appear the software window

![](/images/installation_pip3_09.png)



