# RNA_Tools

## Introduction

In this repository you will find the python script of the software tool called “siRNA and multimeric shRNA designer” (SSD).

## Technical details

* The software was developed using Python v3.5
* The software GUI was developed using PyQt5

## Requirenments to run the script
### **Windows** [Instructive video](https://youtu.be/jG7qKrKMu8M)

Before running the script on Windows you need to download:

1. [Python version 3.5.0](https://www.python.org/ftp/python/3.5.0/python-3.5.0.exe)
2. [PyQt5 version 5.6.0](https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.6/PyQt5-5.6-gpl-Py3.5-Qt5.6.0-x32-2.exe/download)
3. [SSD Software script](/download/siRNA_and_shRNA_designer_(SSD)_release_v1.0.zip)

Instructive video - [How to download the SSD repository](https://www.youtube.com/watch?v=K_eLmiqBpSo&)

The order of the installations must be followed as shown below.

#### Installing Python3.5

These are the files that we have after downloading the previous files

![](/images/installing_python3.5.png)

Double click on the Python installer and follow the steps shown on the screen

![](/images/installing_python3.5_00.png)

Step 1:

![](/images/installing_python3.5_01.png)

Step 2:

![](/images/installing_python3.5_02.png)

Step 3:

![](/images/installing_python3.5_03.png)

#### Installing PyQt5

Double click on the PyQt5 installer and follow the steps shown on the screen

![](/images/installing_pyQT5.png)

Step 1:

![](/images/installing_pyQT5_00.png)

Step 2:

![](/images/installing_pyQT5_01.png)

Step 3:

![](/images/installing_pyQT5_02.png)

Step 4:

![](/images/installing_pyQT5_03.png)

We need to verify that the path where the PyQt5 will be installed looks as follows:

![](/images/pyqt5_path.png)

Step 5:

![](/images/installing_pyQT5_04.png)

#### Running the script on windows

Once you have downloaded the zipped file that is inside the Python script you will see something like this:

![](/images/runnig_the_script_000.png)

First **right click** on the zipped file and press the **Extract Here** option.

![](/images/runnig_the_script_001.png)

Then you will have the Python script unzipped. To run it you just need to **double click** on the script and it will run.

![](/images/runnig_the_script_002.png)

Finally, the software window will appear as shown below:

![](/images/runnig_the_script_06.png)


##### NOTES:
If you cannot run the script by double-clicking on it, try changing the file extension from **.py** to **.pyw**, or just use the **prompt** command  to run it. 

#### Python installation problem: No module named 'encodings'
This is a crash of Python during a script python running caused by un issue amid the Python installation process on Windows, specifically in the PATH environment variable of Python. 

If the crash occurs in the terminal a message output like this:
````
"py_initialize unable to get the locale encoding"
````
The solution for this issue is reinstalling Python following this installation method:

* Installing the exe/msi as admin.
* During installation, select "Add Python 3.x to PATH" and "Customize installation".
* Under "Advanced Options," select "Install for current user only"



[Running the SSD software using the prompt command](https://youtu.be/X0S5jYU3vnU)

### **Linux (Ubuntu)** [Instructive video](https://youtu.be/FC1ttM7NY-0)

Before running the script on Windows, you need to download:

1. PIP3
2. PyQt5
3. [SSD Software script](/download/siRNA_and_shRNA_designer_(SSD)_release_v1.0.zip)

First, you need to check the version of python3 installed in your device by typing the command:
````
python3 --version
````

![](/images/installation_pip3.png)

Once you have typed in the command and pressed enter, the terminal will return the version of python3 that is installed in your device, like this:

![](/images/installation_pip3_00.png)

The next step is to update and upgrade the system in order to install pip3:

````
sudo apt-get update
````

![](/images/installation_pip3_04.png)

````
sudo apt-get upgrade
````

![](/images/installation_pip3_05.png)

The system will show you how to install pip3 by running this command:
````
pip3
````

![](/images/installation_pip3_01.png)

You have to copy, paste and run the command displayed on the terminal.

![](/images/installation_pip3_02.png)

![](/images/installation_pip3_03.png)

You can check the pip3 version using:
````
pip3 -V
````

![](/images/installation_pip3_06.png)

Next, we need to install PyQt5, through the command:
````
pip3 install pyqt5
````

![](/images/installation_pip3_07.png)

#### Running the script on Linux (Ubuntu)

Open the directory where your target script is and right click on it, and select **Open on terminal**.
Then type the following command:
````
python3 siRNA\ and\ shRNA\ designer\ (SSD\)_release_v1.0.py
````

![](/images/installation_pip3_08.png)

Finally, the software window will appear:

![](/images/installation_pip3_09.png)


### **Windows** [How to use SSD software - an example](https://youtu.be/7pfQ7EVX5w8)

[Example file](https://github.com/bioinf2019/RNA_Tools/blob/master/Caenorhabditis%20elegans%20Phosphatidylinositol%203-kinase%20catalytic%20subunit%20type%203%20(vps-34)%2C%20partial%20mRNA.txt) used for the SSD software instructional video. *Caenorhabditis elegans* Phosphatidylinositol 3-kinase catalytic subunit type 3 (vps-34), partial mRNA.


Once you are done installing the software and have the sequence you want to use, you open the software and manually insert the sequence, or click on the button "Open FASTA (.txt)" to import the file containing the sequence. Do not forget to choose a name for the output files and put it on the "Target name" field.

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/upload%20fasta.png)

In order to initiate the strand analysis, click on one of the buttons below according to your needs. When you get a successful strand analysis, you will see a pop-up message like the one shown below:

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/successful%20strand%20analysis.png)

Once you have generated the monovalent sequence file, open it. Go to the monovalent shRNA lid on the SSD software, copy the sense and antisense sequences and paste them on the respective fields. 

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/monovalent.png)

Then press the "Generate shRNA" button, and you will see a pop-up message like the one shown below:

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/monovalent%20success.png)

The resulting monovalent shRNA file will look like this:


![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/01_Monovalent.png)



To generate bivalent and trivalent shRNAs, you must proceed in a similar manner, by going to the bivalent or trivalent shRNA lid of the SSD software respectively, opening the respective nivalent or trivalent sequence file, copying and pasting the sequences into the respective fields, and pressing the "Generate shRNA" button.

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/bivalent.png)

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/bivalent%20success.png)


The bivalent shRNA file will look like this:

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/02_Bivalent.png)


![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/trivalent.png)

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/trivalent%20success.png)


The trivalent shRNA file will look like this (click to enlarge):

![](https://github.com/bioinf2019/RNA_Tools/blob/master/images/03_Trivalent.png)

