#!/usr/bin/env python3
import os
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QLabel, QCheckBox, QWidget, QPushButton, QFrame, QApplication
from PyQt5.QtCore import QSize
from PyQt5.QtGui import QColor
import sys
import math
import datetime
import time

flag = True
#target name
def openPastedSequence(fileName): 
   file = open(fileName, 'r')
   seq_1 = ''   
   for line in file:
      if (line[0] != '>'):
         x = 0
         for x in range(0, len(line)):
            if((line[x] != '\n')):
               seq_1 += line[x]	
   file.close()  
   return seq_1	
   
#validation functiona
def validatorOfCharacters(seq):
   value = False
   for x in range(len(seq)):
      if ((seq[x] != 'A') and (seq[x] != 'T') and (seq[x] != 'G') and (seq[x] != 'C') and (seq[x] != 'U')):
         value = True
   return value
	

def format(value):
   return "%.3f" % value


def bubblesort(list_0, list_1, list_2):
# Swap the elements to arrange in order
    for iter_num in range(len(list_2)-1,0,-1):
        for idx in range(iter_num):
            if list_2[idx]<list_2[idx+1]:
                temp_0 = list_0[idx]			
                temp_1 = list_1[idx]
                temp_2 = list_2[idx]
                list_0[idx] = list_0[idx+1]				
                list_1[idx] = list_1[idx+1]
                list_2[idx] = list_2[idx+1]
                list_0[idx+1] = temp_0				
                list_1[idx+1] = temp_1
                list_2[idx+1] = temp_2	
    return list_1, list_2

def acumGibbsOfEnergy(string):
  GibbsValue = 0.000000
  for x in range(0, len(string) - 1):
    if((string[x] == 'A' and string[x + 1] == 'A') or (string[x] == 'A' and string[x + 1] == 'U') or (string[x] == 'U' and string[x + 1] == 'U')):
      GibbsValue += 1.100000
    elif((string[x] == 'A' and string[x + 1] == 'C') or (string[x] == 'G' and string[x + 1] == 'U')):
      GibbsValue += 2.400000
    elif((string[x] == 'A' and string[x + 1] == 'G') or (string[x] == 'C' and string[x + 1] == 'U')):
      GibbsValue += 1.900000 
    elif((string[x] == 'C' and string[x + 1] == 'A') or (string[x] == 'C' and string[x + 1] == 'G') or (string[x] == 'U' and string[x + 1] == 'G')):
      GibbsValue += 2.200000	  
    elif((string[x] == 'C' and string[x + 1] == 'C') or (string[x] == 'G' and string[x + 1] == 'G')):
      GibbsValue += 3.300000	  
    elif((string[x] == 'G' and string[x + 1] == 'A')):
      GibbsValue += 2.700000	  
    elif((string[x] == 'G' and string[x + 1] == 'C')):
      GibbsValue += 3.800000	  
    elif((string[x] == 'U' and string[x + 1] == 'A')):
      GibbsValue += 1.400000	  
    elif((string[x] == 'U' and string[x + 1] == 'C')):
      GibbsValue += 2.600000	  
  return GibbsValue	
	
def deltaGibbCalculator(string):
  #Gibbs Calculator
  initialSubString = string[2:7]
  initialSubString = initialSubString
  finalSubString = string[16:21] 
  finalSubString = finalSubString
  finalSubString = invertSeq(finalSubString)
  finalSubString = complementSeqUracil(finalSubString)
  initialGibbs = acumGibbsOfEnergy(initialSubString)
  finalGibbs = acumGibbsOfEnergy(finalSubString)  
  deltaGibbs = initialGibbs - finalGibbs
  return(deltaGibbs)	
	
#Read FASTA file

def read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))	

def invertSeq(seq):
  return seq[::-1]

def moveNucleotides(seq):
  seq = seq[19] + seq[20] + seq[0:19]
  return seq
  
def complementSeq(seq):
  compSeq = ''
  for index in seq:
    if index == 'A':
      compSeq += 'T'
    elif index == 'T':
      compSeq += 'A'
    elif index == 'C':
      compSeq += 'G'
    elif index == 'G':
      compSeq += 'C'  
    elif index == '-':
      compSeq += '-'  	  
  return compSeq  
  
def complementSeqUracil(seq):
  compSeq = ''
  for index in seq:
    if index == 'A':
      compSeq += 'U'
    elif index == 'U':
      compSeq += 'A'
    elif index == 'C':
      compSeq += 'G'
    elif index == 'G':
      compSeq += 'C'  
    elif index == '-':
      compSeq += '-'  	  
  return compSeq  

def multimeric_shRNA(siRNAsSenseStrandFiveThree_21_mer, siRNAsAntienseStrandFiveThree_21_mer, minorLoop, largerLoop): 
  reverseComplementOfTheLastTwoNucleotides = siRNAsAntienseStrandFiveThree_21_mer[-2] + siRNAsAntienseStrandFiveThree_21_mer[-1]
  reverseComplementOfTheLastTwoNucleotides = complementSeqUracil(reverseComplementOfTheLastTwoNucleotides)
  sequence_23_mer = invertSeq(reverseComplementOfTheLastTwoNucleotides) + siRNAsSenseStrandFiveThree_21_mer
  replacedNucleotidesSequence = sequence_23_mer.replace('U','T')
  complementSequence = complementSeq(replacedNucleotidesSequence) 
  invertedComplementSequence = invertSeq(complementSequence)
  completeLinearSequenceGeneratorOf_shRNA = minorLoop + replacedNucleotidesSequence + largerLoop + invertedComplementSequence
  return completeLinearSequenceGeneratorOf_shRNA
#MessageBox
class MessageBox(QMessageBox):
  changedValue = pyqtSignal(int)
  
  def __ini__(self2):
    QMessageBox.__init__(self)
    self2.setText("This is a MessageBox, typically used to convey short messages to the user.")
    self2.setInformativeText("Informative text provides more space to explain the message purpose.")
    self2.setIcon(QMessageBox.Information)
    self2.setStandardButtons(QMessageBox.Close)
	
#Main window 
class Window(QWidget):
   def __init__(self, parent=None):
      super(Window, self).__init__(parent)      
      self.setMinimumSize(1360, 600)
      self.setMaximumSize(1360, 600)
	  
      QWidget.__init__(self)
      layout = QGridLayout()
      self.setLayout(layout)
      self.tab4 = QWidget()
      self.tab1 = QWidget()
      self.tab2 = QWidget()
      self.tab3 = QWidget()
      self.tab5 = QWidget()
      tabwidget = QTabWidget()
      tabwidget.addTab(self.tab4, "siRNA designer")
      tabwidget.addTab(self.tab1, "Monovalent shRNA")
      tabwidget.addTab(self.tab2, "Bivalent shRNA")
      tabwidget.addTab(self.tab3, "Trivalent shRNA")
      tabwidget.addTab(self.tab5, "General Instructions")	  
      layout.addWidget(tabwidget, 100, 100)
      self.tab4UI()
      self.tab1UI()
      self.tab2UI()
      self.tab3UI()
      self.tab5UI()	  
      self.setWindowTitle("siRNA and shRNA designer (SSD)V1.0")

   def tab1UI(self):
      layout=QFormLayout()
      self.nameLine_001 = QLineEdit()	  
      self.nameLine_1 = QLineEdit()
      self.nameLine_2 = QLineEdit()
      self.nameLine_3 = QLineEdit()
      self.nameLine_4 = QLineEdit()
      layout.addRow("Target name:",self.nameLine_001)	  
      label_1 = QLabel("Insert Sequences of One siRNA (21 nt; 5' to 3'):")
      layout.addRow(label_1)	  
      layout.addRow("siRNA's sense strand:",self.nameLine_1)
      layout.addRow("siRNA's antisense strand:",self.nameLine_2)
      label_2 = QLabel("OR\n\nOPEN FASTA FILE (the same sequences, in the same order, as above.):")
      layout.addRow(label_2)	  
      self.openFileButton_1 = QPushButton("Open FASTA file")
      self.generateButton_1 = QPushButton("Generate shRNA")	  
      layout.addRow(self.openFileButton_1, self.nameLine_3)
      layout.addRow(self.generateButton_1)
      label_4 = QLabel("Monovalent Multmeric shRNA sequence:")
      layout.addRow(label_4)	  
      layout.addRow(self.nameLine_4)	  
      label_3 = QLabel("RESULT FILES ARE GENERATED IN SOFTWARE'S ORIGINAL FOLDER")
      layout.addRow(label_3)
	  #Input validation
      self.generateButton_1.setEnabled(False)
      self.nameLine_1.textChanged[str].connect(self.your_function_1)
	  #Conditions When buttons are pressed
      self.openFileButton_1.clicked.connect(self.open)
      self.generateButton_1.clicked.connect(self.generateMultimeric)  
      self.tab1.setLayout(layout)

   def your_function_1(self,text):
      if not text:
         self.generateButton_1.setEnabled(False)
      else:
         self.generateButton_1.setEnabled(True)	  
	  
   def open(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_3.setText(fileName[0])
      if not os.path.isfile(fileName[0]):
         messageButton = QMessageBox.about(self, "Error opening file", "You did not select a file, please select a file")
      else:	  
         with open(fileName[0]) as fp:
            nameList = []
            seqList = []
            for name, seq in read_fasta(fp):
               nameList.append(name)
               seqList.append(seq)
         sense = seqList[0]
         antisense = seqList[1]	   
         self.nameLine_1.setText(sense)
         self.nameLine_2.setText(antisense)	   
	  
   def generateMultimeric(self):     
      siRNAsSenseStrandFiveThree_21_mer = self.nameLine_1.text()
      siRNAsAntienseStrandFiveThree_21_mer = self.nameLine_2.text()
	  #validation of characters
      value_1 = validatorOfCharacters(siRNAsSenseStrandFiveThree_21_mer)
      value_2 = validatorOfCharacters(siRNAsAntienseStrandFiveThree_21_mer)
      now = datetime.datetime.now()
      outputName = self.nameLine_001.text()
      if (outputName == ''):
         QMessageBox.about(self, "Warning", "Please enter a Target name")
      else:	  
         if ((siRNAsAntienseStrandFiveThree_21_mer == '') or siRNAsSenseStrandFiveThree_21_mer == ''):
            QMessageBox.about(self, "Warning", "Please complete with all the required sequences")
         else:
            if (value_1 == True or value_2 == True):
               self.nameLine_1.setText('')
               self.nameLine_2.setText('') 
               self.generateButton_1.setEnabled(False)		 
               QMessageBox.about(self, "Warning", "Error: Contains illegal characters\nYou must enter a secuence with this nucleotides 'AGTCU'")
            else:
               reverseComplementOfTheLastTwoNucleotides = siRNAsAntienseStrandFiveThree_21_mer[-2] + siRNAsAntienseStrandFiveThree_21_mer[-1]
               reverseComplementOfTheLastTwoNucleotides = complementSeqUracil(reverseComplementOfTheLastTwoNucleotides)
               sequence_23_mer = invertSeq(reverseComplementOfTheLastTwoNucleotides) + siRNAsSenseStrandFiveThree_21_mer
               replacedNucleotidesSequence = sequence_23_mer.replace('U','T')  
               complementSequence = complementSeq(replacedNucleotidesSequence) 
               invertedComplementSequence = invertSeq(complementSequence) 
               minorLoop = 'CTCTC'
               largerLoop = 'TGACAGGAAG'
               completeLinearSequenceGeneratorOf_shRNA = minorLoop + replacedNucleotidesSequence + largerLoop + invertedComplementSequence
               inicialPortion = completeLinearSequenceGeneratorOf_shRNA[:16]
               finalPortion = completeLinearSequenceGeneratorOf_shRNA[16:]
               shRNA = finalPortion + inicialPortion	 
               saveFile = open(outputName + '_Monovalent_shRNA_' + now.strftime("SSD_%Y-%m-%d_[%H h %M min].txt"), 'w')
               saveFile.write('********************************************************************************\n')
               saveFile.write('                         siRNA and shRNA designer (SSD)\n')
               saveFile.write('********************************************************************************\n')
               saveFile.write('Version 1.0 - 2018\n')
               saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
               saveFile.write('********************************************************************************\n')
               saveFile.write('>Monovalent_shRNA\n')
               saveFile.write(shRNA)
               saveFile.close()
               self.nameLine_4.setText(shRNA)
               #Success Message	  
               QMessageBox.about(self, "Process", "The Monovalent shRNA sequence was successfully generated")

   def automaticGenerateMultimeric_shRNA(self):
      rawSequence = self.nameLine_20.toPlainText()
      now = datetime.datetime.now()		 
      outputName = self.nameLine_19.text()
      if (outputName == ''):
         QMessageBox.about(self, "Warning", "Please enter a Target name")
      else:	
	     #SAve the pasted sequence in a file
         sequenceName = outputName + now.strftime("_used_sequence_SSD %Y-%m-%d_[%H h %M min].txt")	  
         saveFile = open(sequenceName, 'w')	
         saveFile.write('>' + outputName + '\n')		 
         saveFile.write(rawSequence)
         saveFile.close()  		 
         mRNA_Sequence = openPastedSequence(sequenceName)		 
         #validation of characters	  
         value_1 = validatorOfCharacters(mRNA_Sequence)
         if (value_1 == True):
            self.nameLine_20.setText('')
            self.checkbox1.setEnabled(False)
            self.checkbox2.setEnabled(False)
            self.checkbox3.setEnabled(False)		 
            QMessageBox.about(self, "Warning", "Error: Contains illegal characters\nYou must enter a secuence with this nucleotides 'AGTCU'")
         else:   
            #Generate siRNAs ordered by Gibbs free energy
            mRNA_Sequence = '--' + mRNA_Sequence
            mRNA_Sequence = mRNA_Sequence.replace('T','U')
            first_nt = 0
            last_nt = 23
            flag = True
            sequence = []
            deltaGibbs = []
            while(flag == True):
               if (last_nt <= (len(mRNA_Sequence))):
                  sequence.append(mRNA_Sequence[first_nt:last_nt])
                  deltaGibbs.append(deltaGibbCalculator(mRNA_Sequence[first_nt:last_nt]))
               else:
                  flag = False
               first_nt+=1
               last_nt+=1
	        #This part sorts the sequences by gibbs free energy 
            positionOfSeq = []		 
            for x in range(len(sequence)):
               positionOfSeq.append(x)     
            bubblesort(positionOfSeq,sequence, deltaGibbs)
            sense = []
            antisense = []
            for x in range(len(sequence)):
               sense.append(sequence[x][2:24])
               antisense.append(complementSeqUracil(invertSeq(sequence[x][0:21])))
            x = 0
	        #Send the first pair sequences
	        #To monovalent tab
            siRNAsSenseStrandFiveThree_21_mer = sense[x]
            siRNAsAntienseStrandFiveThree_21_mer = antisense[x]	  
	        #To bivalent tab
            siRNAsSenseStrandFiveThree_21_merNumberOne = sense[x]
            siRNAsAntienseStrandFiveThree_21_merNumberOne = antisense[x]	  
	        #To trivalent tab
            siRNAsSenseStrandFiveThree_21_merNumberOne = sense[x]
            siRNAsAntienseStrandFiveThree_21_merNumberOne = antisense[x]	  		 
            y = x + 1
            while(abs(positionOfSeq[x] - positionOfSeq[y]) <= 23):
               y = y + 1
	        #Send the second pair sequences
	        #To bivalent tab
            siRNAsSenseStrandFiveThree_21_merNumberTwo = sense[y]
            siRNAsAntienseStrandFiveThree_21_merNumberTwo = antisense[y]	  
	        #To trivalent tab
            siRNAsSenseStrandFiveThree_21_merNumberTwo = sense[y]
            siRNAsAntienseStrandFiveThree_21_merNumberTwo = antisense[y]	  	
            y = y + 1		 
            while(abs(positionOfSeq[x] - positionOfSeq[y]) <= 23):
               y = y + 1
            #Send the third pair sequences
	        #To trivalent tab
            siRNAsSenseStrandFiveThree_21_merNumberThree = sense[y]
            siRNAsAntienseStrandFiveThree_21_merNumberThree = antisense[y]	  	      
            #Monovalent shRNA generation
            reverseComplementOfTheLastTwoNucleotides = siRNAsAntienseStrandFiveThree_21_mer[-2] + siRNAsAntienseStrandFiveThree_21_mer[-1]
            reverseComplementOfTheLastTwoNucleotides = complementSeqUracil(reverseComplementOfTheLastTwoNucleotides)
            sequence_23_mer = invertSeq(reverseComplementOfTheLastTwoNucleotides) + siRNAsSenseStrandFiveThree_21_mer
            replacedNucleotidesSequence = sequence_23_mer.replace('U','T')  
            complementSequence = complementSeq(replacedNucleotidesSequence) 
            invertedComplementSequence = invertSeq(complementSequence) 
            minorLoop = 'CTCTC'
            largerLoop = 'TGACAGGAAG'
            completeLinearSequenceGeneratorOf_shRNA = minorLoop + replacedNucleotidesSequence + largerLoop + invertedComplementSequence
            inicialPortion = completeLinearSequenceGeneratorOf_shRNA[:16]
            finalPortion = completeLinearSequenceGeneratorOf_shRNA[16:]
            shRNA = finalPortion + inicialPortion
            saveFile = open(outputName + '_Monovalent_shRNA_' + now.strftime("SSD_%Y-%m-%d_[%H h %M min].txt"), 'w')
            saveFile.write('********************************************************************************\n')
            saveFile.write('                         siRNA and shRNA designer (SSD)\n')
            saveFile.write('********************************************************************************\n')
            saveFile.write('Version 1.0 - 2018\n')
            saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
            saveFile.write('********************************************************************************\n')
            saveFile.write('>Monovalent_shRNA\n')
            saveFile.write(shRNA)
            saveFile.close()
            #Bivalent shRNA generation	
            minorLoopNumberOne = 'CTCTC'
            minorLoopNumberTwo = 'CTCTC'
            largerLoopNumberAll = 'TGACAGGAAG'
            completeLinearSequence = multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberOne, siRNAsAntienseStrandFiveThree_21_merNumberOne, minorLoopNumberOne, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberTwo, siRNAsAntienseStrandFiveThree_21_merNumberTwo, minorLoopNumberTwo, largerLoopNumberAll)
            inicialPortion = completeLinearSequence[:16]
            finalPortion = completeLinearSequence[16:]
            BivalentMultimeric_shRNA = finalPortion + inicialPortion	  
            saveFile = open(outputName + '_Bivalent_shRNA_' + now.strftime("SSD_%Y-%m-%d_[%H h %M min].txt"), 'w')
            saveFile.write('********************************************************************************\n')
            saveFile.write('                         siRNA and shRNA designer (SSD)\n')
            saveFile.write('********************************************************************************\n')
            saveFile.write('Version 1.0 - 2018\n')
            saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
            saveFile.write('********************************************************************************\n')
            saveFile.write('>Bivalent_shRNA\n')
            saveFile.write(BivalentMultimeric_shRNA)
            saveFile.close() 
            #Trivalent shRNA generation  	  
            minorLoopNumberOne = 'CTCTC'
            minorLoopNumberTwo = 'CTCTC'
            minorLoopNumberThree = 'CTCTC'
            largerLoopNumberAll = 'TGACAGGAAG'
            completeLinearSequence = multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberOne, siRNAsAntienseStrandFiveThree_21_merNumberOne, minorLoopNumberOne, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberTwo, siRNAsAntienseStrandFiveThree_21_merNumberTwo, minorLoopNumberTwo, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberThree, siRNAsAntienseStrandFiveThree_21_merNumberThree, minorLoopNumberThree, largerLoopNumberAll)
            inicialPortion = completeLinearSequence[:16]
            finalPortion = completeLinearSequence[16:]	  
            TrivalentMultimeric_shRNA = finalPortion + inicialPortion
            saveFile = open(outputName + '_Trivalent_shRNA_' + now.strftime("SSD_%Y-%m-%d_[%H h %M min].txt"), 'w')
            saveFile.write('********************************************************************************\n')
            saveFile.write('                         siRNA and shRNA designer (SSD)\n')
            saveFile.write('********************************************************************************\n')
            saveFile.write('Version 1.0 - 2018\n')
            saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
            saveFile.write('********************************************************************************\n')
            saveFile.write('>Trivalent_shRNA\n')	  
            saveFile.write(TrivalentMultimeric_shRNA)
            saveFile.close()   
            #Success Message		  
            QMessageBox.about(self, "Process(automatic generation of shRNAs)", "The generation of Monovalent, Bivalent and Trivalent shRNAs was successfully done")		  	  
	  
   def tab2UI(self):
      layout=QFormLayout()
      self.nameLine_002 = QLineEdit()	  
      self.nameLine_5 = QLineEdit()
      self.nameLine_6 = QLineEdit()
      self.nameLine_7 = QLineEdit()
      self.nameLine_8 = QLineEdit()
      self.nameLine_9 = QLineEdit()
      self.nameLine_10 = QLineEdit()	
      layout.addRow("Target name:",self.nameLine_002)	  
      label_5 = QLabel("Insert Sequences of Two siRNAs (21 nt; 5' to 3'):")
      layout.addRow(label_5)
      layout.addRow("First siRNA's sense strand:",self.nameLine_5)
      layout.addRow("First siRNA's antsense strand:",self.nameLine_6)
      layout.addRow("Second siRNA's sense strand:",self.nameLine_7)
      layout.addRow("Second siRNA's antsense strand:",self.nameLine_8)	  
      label_5 = QLabel("OR\n\nOPEN FASTA FILE (the same sequences, in the same order, as above.):")
      layout.addRow(label_5)
      self.openFileButton_2 = QPushButton("Open FASTA file")
      self.generateButton_2 = QPushButton("Generate shRNA")	  
      layout.addRow(self.openFileButton_2, self.nameLine_9)
      layout.addRow(self.generateButton_2)
      label_6 = QLabel("Bivalent Multmeric shRNA sequence:")
      layout.addRow(label_6)	
      layout.addRow(self.nameLine_10)  
      label_7 = QLabel("RESULT FILES ARE GENERATED IN SOFTWARE'S ORIGINAL FOLDER")
      layout.addRow(label_7)
	  #Validation
      self.generateButton_2.setEnabled(False)
      self.nameLine_5.textChanged[str].connect(self.your_function_2)
	  #Conditions When buttons are pressed
      self.openFileButton_2.clicked.connect(self.open2)
      self.generateButton_2.clicked.connect(self.generateMultimeric2)	  
      self.tab2.setLayout(layout)	  

   def your_function_2(self,text):
      if not text:
         self.generateButton_2.setEnabled(False)
      else:
         self.generateButton_2.setEnabled(True)
	  
   def open2(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_9.setText(fileName[0])
      if not os.path.isfile(fileName[0]):
         messageButton = QMessageBox.about(self, "Error opening file", "You did not select a file, please select a file")
      else:			
         with open(fileName[0]) as fp:
            nameList = []
            seqList = []
            for name, seq in read_fasta(fp):
               nameList.append(name)
               seqList.append(seq)
         sense_1 = seqList[0]
         antisense_1 = seqList[1]	   
         sense_2 = seqList[2]
         antisense_2 = seqList[3]
         self.nameLine_5.setText(sense_1)
         self.nameLine_6.setText(antisense_1)
         self.nameLine_7.setText(sense_2)
         self.nameLine_8.setText(antisense_2)	  

   def generateMultimeric2(self):
      siRNAsSenseStrandFiveThree_21_merNumberOne = self.nameLine_5.text()
      siRNAsAntienseStrandFiveThree_21_merNumberOne = self.nameLine_6.text()
      siRNAsSenseStrandFiveThree_21_merNumberTwo = self.nameLine_7.text()
      siRNAsAntienseStrandFiveThree_21_merNumberTwo = self.nameLine_8.text()
      now = datetime.datetime.now()	
      outputName = self.nameLine_002.text()	
      if (outputName == ''):
         QMessageBox.about(self, "Warning", "Please enter a Target name")
      else:
         if ((siRNAsAntienseStrandFiveThree_21_merNumberOne == '') or (siRNAsSenseStrandFiveThree_21_merNumberOne == '') or (siRNAsAntienseStrandFiveThree_21_merNumberTwo == '') or (siRNAsSenseStrandFiveThree_21_merNumberTwo == '')):
            QMessageBox.about(self, "Warning", "Please complete with all the required sequences")
         else:   	  
	        #validation of characters
            value_1 = validatorOfCharacters(siRNAsSenseStrandFiveThree_21_merNumberOne)
            value_2 = validatorOfCharacters(siRNAsAntienseStrandFiveThree_21_merNumberOne)
            value_3 = validatorOfCharacters(siRNAsSenseStrandFiveThree_21_merNumberTwo)
            value_4 = validatorOfCharacters(siRNAsAntienseStrandFiveThree_21_merNumberTwo)	  
            if (value_1 == True or value_2 == True or value_3 == True or value_4 == True):
               self.nameLine_5.setText('')
               self.nameLine_6.setText('')
               self.nameLine_7.setText('')
               self.nameLine_8.setText('')
               self.generateButton_2.setEnabled(False)		 
               QMessageBox.about(self, "Warning", "Error: Contains illegal characters\nYou must enter a secuence with this nucleotides 'AGTCU'")
            else:	  
               minorLoopNumberOne = 'CTCTC'
               minorLoopNumberTwo = 'CTCTC'
               largerLoopNumberAll = 'TGACAGGAAG'
               completeLinearSequence = multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberOne, siRNAsAntienseStrandFiveThree_21_merNumberOne, minorLoopNumberOne, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberTwo, siRNAsAntienseStrandFiveThree_21_merNumberTwo, minorLoopNumberTwo, largerLoopNumberAll)
               inicialPortion = completeLinearSequence[:16]
               finalPortion = completeLinearSequence[16:]
               BivalentMultimeric_shRNA = finalPortion + inicialPortion	
               saveFile = open(outputName + '_Bivalent_shRNA_' + now.strftime("SSD_%Y-%m-%d_[%H h %M min].txt"), 'w')
               saveFile.write('********************************************************************************\n')
               saveFile.write('                         siRNA and shRNA designer (SSD)\n')
               saveFile.write('********************************************************************************\n')
               saveFile.write('Version 1.0 - 2018\n')
               saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
               saveFile.write('********************************************************************************\n')
               saveFile.write('>Bivalent_shRNA\n')
               saveFile.write(BivalentMultimeric_shRNA)
               saveFile.close()
               self.nameLine_10.setText(BivalentMultimeric_shRNA)  	  
	           #Success message
               QMessageBox.about(self, "Process", "The Bivalent shRNA sequence was successfully generated")
	  
   def tab3UI(self):
      layout=QFormLayout()
      self.nameLine_003 = QLineEdit()
      self.nameLine_11 = QLineEdit()
      self.nameLine_12 = QLineEdit()
      self.nameLine_13 = QLineEdit()
      self.nameLine_14 = QLineEdit()
      self.nameLine_15 = QLineEdit()
      self.nameLine_16 = QLineEdit()	  
      self.nameLine_17 = QLineEdit()
      self.nameLine_18 = QLineEdit()	
      layout.addRow("Target name:",self.nameLine_003)	  
      label_8 = QLabel("Insert Sequences of Three siRNAs (21 nt; 5' to 3'):")
      layout.addRow(label_8)
      layout.addRow("First siRNA's sense strand:",self.nameLine_11)
      layout.addRow("First siRNA's antsense strand:",self.nameLine_12)
      layout.addRow("Second siRNA's sense strand:",self.nameLine_13)
      layout.addRow("Second siRNA's antsense strand:",self.nameLine_14)	  
      layout.addRow("Third siRNA's sense strand:",self.nameLine_15)
      layout.addRow("Third siRNA's antsense strand:",self.nameLine_16)	
      label_9 = QLabel("OR\n\nOPEN FASTA FILE (the same sequences, in the same order, as above.):")
      layout.addRow(label_9)	  
      self.openFileButton_3 = QPushButton("Open FASTA file")
      self.generateButton_3 = QPushButton("Generate shRNA")	  
      layout.addRow(self.openFileButton_3, self.nameLine_17)
      layout.addRow(self.generateButton_3)
      label_10 = QLabel("Trivalent Multmeric shRNA sequence:")
      layout.addRow(label_10)
      layout.addRow(self.nameLine_18)	  
      label_11 = QLabel("RESULT FILES ARE GENERATED IN SOFTWARE'S ORIGINAL FOLDER")
      layout.addRow(label_11)
	  #Validation
      self.generateButton_3.setEnabled(False)
      self.nameLine_11.textChanged[str].connect(self.your_function_3)	  
	  #Conditions When buttons are pressed
      self.openFileButton_3.clicked.connect(self.open3)
      self.generateButton_3.clicked.connect(self.generateMultimeric3)
      self.tab3.setLayout(layout)

   def your_function_3(self,text):
      if not text:
         self.generateButton_3.setEnabled(False)
      else:
         self.generateButton_3.setEnabled(True)	  
	  
   def open3(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_17.setText(fileName[0])
      if not os.path.isfile(fileName[0]):
         messageButton = QMessageBox.about(self, "Error opening file", "You did not select a file, please select a file")
      else:	  
         with open(fileName[0]) as fp:
            nameList = []
            seqList = []
            for name, seq in read_fasta(fp):
               nameList.append(name)
               seqList.append(seq)
         sense_1 = seqList[0]
         antisense_1 = seqList[1]	   
         sense_2 = seqList[2]
         antisense_2 = seqList[3]
         sense_3 = seqList[4]
         antisense_3 = seqList[5]	  
         self.nameLine_11.setText(sense_1)
         self.nameLine_12.setText(antisense_1)
         self.nameLine_13.setText(sense_2)
         self.nameLine_14.setText(antisense_2)
         self.nameLine_15.setText(sense_3)
         self.nameLine_16.setText(antisense_3)

   def generateMultimeric3(self):
      siRNAsSenseStrandFiveThree_21_merNumberOne = self.nameLine_11.text()
      siRNAsAntienseStrandFiveThree_21_merNumberOne = self.nameLine_12.text()
      siRNAsSenseStrandFiveThree_21_merNumberTwo = self.nameLine_13.text()
      siRNAsAntienseStrandFiveThree_21_merNumberTwo = self.nameLine_14.text()
      siRNAsSenseStrandFiveThree_21_merNumberThree = self.nameLine_15.text()
      siRNAsAntienseStrandFiveThree_21_merNumberThree = self.nameLine_16.text()
      outputName = self.nameLine_003.text()
      now = datetime.datetime.now()
      if (outputName == ''):
         QMessageBox.about(self, "Warning", "Please enter a Target name")
      else:
         if ((siRNAsAntienseStrandFiveThree_21_merNumberOne == '') or (siRNAsSenseStrandFiveThree_21_merNumberOne == '') or (siRNAsAntienseStrandFiveThree_21_merNumberTwo == '') or (siRNAsSenseStrandFiveThree_21_merNumberTwo == '') or (siRNAsAntienseStrandFiveThree_21_merNumberThree == '') or (siRNAsSenseStrandFiveThree_21_merNumberThree == '')):
            QMessageBox.about(self, "Warning", "Please complete with all the required sequences")
         else: 	  
	        #validation of characters
            value_1 = validatorOfCharacters(siRNAsSenseStrandFiveThree_21_merNumberOne)
            value_2 = validatorOfCharacters(siRNAsAntienseStrandFiveThree_21_merNumberOne)
            value_3 = validatorOfCharacters(siRNAsSenseStrandFiveThree_21_merNumberTwo)
            value_4 = validatorOfCharacters(siRNAsAntienseStrandFiveThree_21_merNumberTwo)
            value_5 = validatorOfCharacters(siRNAsSenseStrandFiveThree_21_merNumberThree)
            value_6 = validatorOfCharacters(siRNAsAntienseStrandFiveThree_21_merNumberThree)	  
            if (value_1 == True or value_2 == True or value_3 == True or value_4 == True or value_5 == True or value_6 == True):
               self.nameLine_11.setText('')
               self.nameLine_12.setText('')
               self.nameLine_13.setText('')
               self.nameLine_14.setText('')	
               self.nameLine_15.setText('')
               self.nameLine_16.setText('')	
               self.generateButton_3.setEnabled(False)		 
               QMessageBox.about(self, "Warning", "Error: Contains illegal characters\nYou must enter a secuence with this nucleotides 'AGTCU'")
            else:	 	  
               minorLoopNumberOne = 'CTCTC'
               minorLoopNumberTwo = 'CTCTC'
               minorLoopNumberThree = 'CTCTC'
               largerLoopNumberAll = 'TGACAGGAAG'
               completeLinearSequence = multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberOne, siRNAsAntienseStrandFiveThree_21_merNumberOne, minorLoopNumberOne, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberTwo, siRNAsAntienseStrandFiveThree_21_merNumberTwo, minorLoopNumberTwo, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberThree, siRNAsAntienseStrandFiveThree_21_merNumberThree, minorLoopNumberThree, largerLoopNumberAll)
               inicialPortion = completeLinearSequence[:16]
               finalPortion = completeLinearSequence[16:]	  
               TrivalentMultimeric_shRNA = finalPortion + inicialPortion         		 
               saveFile = open(outputName + '_Trivalent_shRNA_' + now.strftime("SSD_%Y-%m-%d_[%H h %M min].txt"), 'w')
               saveFile.write('********************************************************************************\n')
               saveFile.write('                         siRNA and shRNA designer (SSD)\n')
               saveFile.write('********************************************************************************\n')
               saveFile.write('Version 1.0 - 2018\n')
               saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
               saveFile.write('********************************************************************************\n')
               saveFile.write('>Trivalent_shRNA\n')	  
               saveFile.write(TrivalentMultimeric_shRNA)
               saveFile.close()
               self.nameLine_18.setText(TrivalentMultimeric_shRNA)  	  
               #Success Message	  
               QMessageBox.about(self, "Process", "The Trivalent shRNA sequence was successfully generated")	  

   def tab4UI(self):
      layout=QFormLayout()
      self.nameLine_19 = QLineEdit()	
      self.nameLine_20 = QTextEdit()
      self.nameLine_21 = QLineEdit()	  
      self.checkbox1 = QPushButton("Generate siRNAs Ordered by Positon within mRNA")
      self.checkbox2 = QPushButton("Generate siRNAs Ordered by Gibbs free energy value")  
      self.checkbox3 = QPushButton("Generate shRNAs automatically from mRNA")	  
      layout.addRow("Target name:",self.nameLine_19)
      layout.addRow("Insert full or partial mRNA (FASTA FORMAT):",self.nameLine_20)
      self.openFileButton_4 = QPushButton("Or Open FASTA (.txt)")
#      self.generateButton_4 = QPushButton("Generate siRNAs")	  
      layout.addRow(self.openFileButton_4, self.nameLine_21)
      label2 = QLabel("======================\nGENERATION OF siRNAs BY ORDER OF PRESENTATION")
      layout.addRow(label2)	  
      layout.addRow(self.checkbox1)	
      layout.addRow(self.checkbox2)	
      label02 = QLabel("======================\nAUTOMATIC GENERATION OF shRNAs FROM mRNA")
      layout.addRow(label02)	 	  
      layout.addRow(self.checkbox3)	  
      label3 = QLabel("======================")	  
      layout.addRow(label3)  	  
#      layout.addRow(self.generateButton_4)
      label1 = QLabel("RESULT FILES ARE GENERATED IN SOFTWARE'S ORIGINAL FOLDER")
      layout.addRow(label1)
	  #Validation
      self.checkbox1.setEnabled(False)
      self.checkbox2.setEnabled(False)
      self.checkbox3.setEnabled(False)
      self.nameLine_20.textChanged.connect(self.your_function_4)	  
	  #Conditions When buttons are pressed
      self.openFileButton_4.clicked.connect(self.open4)
      self.checkbox1.clicked.connect(self.generate_siRNAs_Order_by_position)
      self.checkbox2.clicked.connect(self.generate_siRNAs_Order_by_Gibbs_free_energy)	
      self.checkbox3.clicked.connect(self.automaticGenerateMultimeric_shRNA)	  
      self.tab4.setLayout(layout)

   def your_function_4(self):
      self.checkbox1.setEnabled(True)
      self.checkbox2.setEnabled(True)
      self.checkbox3.setEnabled(True)		 
	  
   def open4(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_21.setText(fileName[0])
      if not os.path.isfile(fileName[0]):
         messageButton = QMessageBox.about(self, "Error opening file", "You did not select a file, please select a file")
      else:	  
         with open(fileName[0]) as fp:
            nameList = []
            seqList = []
            for name, seq in read_fasta(fp):
               nameList.append(name)
               seqList.append(seq)
         seq_1 = seqList[0]	  
         self.nameLine_20.setText(seq_1)
		 
 

   def generate_siRNAs_Order_by_position(self):
      rawSequence = self.nameLine_20.toPlainText()
      now = datetime.datetime.now()		 
      outputName = self.nameLine_19.text()
      if (outputName == ''):
         QMessageBox.about(self, "Warning", "Please enter a Target name")
      else:	
	     #SAve the pasted sequence in a file
         sequenceName = outputName + now.strftime("_used_sequence_SSD %Y-%m-%d_[%H h %M min].txt")	  
         saveFile = open(sequenceName, 'w')	
         saveFile.write('>' + outputName + '\n')		 
         saveFile.write(rawSequence)
         saveFile.close()  		 
         mRNA_Sequence = openPastedSequence(sequenceName)  		 
         #validation of characters	  
         value_1 = validatorOfCharacters(mRNA_Sequence)
         if (value_1 == True):
            self.nameLine_19.setText('')		 
            self.nameLine_20.setText('') 
            self.checkbox1.setEnabled(False)
            self.checkbox2.setEnabled(False)
            self.checkbox3.setEnabled(False)		 
            QMessageBox.about(self, "Warning", "Error: Contains illegal characters\nYou must enter a secuence with this nucleotides 'AGTCU'")
         else:   
            mRNA_Sequence = '--' + mRNA_Sequence
            mRNA_Sequence = mRNA_Sequence.replace('T','U')
            first_nt = 0
            last_nt = 23
            flag = True
            sequence = []
            deltaGibbs = []
            while(flag == True):
               if (last_nt <= (len(mRNA_Sequence))):
                  sequence.append(mRNA_Sequence[first_nt:last_nt])
                  deltaGibbs.append(deltaGibbCalculator(mRNA_Sequence[first_nt:last_nt]))
               else:
                  flag = False
               first_nt+=1
               last_nt+=1
            sense = []
            antisense = []
            for x in range(len(sequence)):
               sense.append(sequence[x][2:24])
               antisense.append(complementSeqUracil(invertSeq(sequence[x][0:21])))		 	 
            saveFile = open(outputName + now.strftime("_siRNA_position_SSD %Y-%m-%d_[%H h %M min].txt"), 'w')
            saveFile.write('********************************************************************************\n')
            saveFile.write('                         siRNA and shRNA designer (SSD)\n')
            saveFile.write('********************************************************************************\n')
            saveFile.write('Version 1.0 - 2018\n')
            saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
            saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
            saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
            saveFile.write('siRNAs below are ordered by position.\n')
            saveFile.write("All sequences are oriented 5' to 3'.\n")
            saveFile.write('********************************************************************************\n')
            saveFile.write('Gene name: '+outputName+'\n')
            saveFile.write('********************************************************************************\n')	  
            for x in range(len(sense)):
               saveFile.write('POSITION:'+ str(x + 1) +'\n')
               saveFile.write('\nSense RNA:	')
               saveFile.write(sense[x])
               saveFile.write('\nAntisense RNA:	')
               saveFile.write(antisense[x])
               saveFile.write('\n\n')
               saveFile.write('RESULTS:\n')
               saveFile.write('\nDelta Gibbs: ')
               formatted = format(deltaGibbs[x])
               saveFile.write(str(formatted))
               if(deltaGibbs[x] >= 0):
                  saveFile.write(' => Functional siRNA\n')
               else:
                  saveFile.write(' => Non-Functional siRNA\n')
               saveFile.write('-----------------------------------------------------------------------\n\n')			
            saveFile.close()		  
            #Success Message  
            QMessageBox.about(self, "Process(Ordered by position)", "The Strand Analysis was successfully done")

   def generate_siRNAs_Order_by_Gibbs_free_energy(self):
      rawSequence = self.nameLine_20.toPlainText()
      now = datetime.datetime.now()		 
      outputName = self.nameLine_19.text()
      if (outputName == ''):
         QMessageBox.about(self, "Warning", "Please enter a Target name")
      else:	
	     #SAve the pasted sequence in a file
         sequenceName = outputName + now.strftime("_used_sequence_SSD %Y-%m-%d_[%H h %M min].txt")	  
         saveFile = open(sequenceName, 'w')	
         saveFile.write('>' + outputName + '\n')		 
         saveFile.write(rawSequence)
         saveFile.close()  		 
         mRNA_Sequence = openPastedSequence(sequenceName)   
         #validation of characters	   	  
         value_1 = validatorOfCharacters(mRNA_Sequence)
         if (value_1 == True):
            self.nameLine_20.setText('') 
            self.checkbox1.setEnabled(False)
            self.checkbox2.setEnabled(False)
            self.checkbox3.setEnabled(False)		 
            QMessageBox.about(self, "Warning", "Error: Contains illegal characters\nYou must enter a secuence with this nucleotides 'AGTCU'")
         else:   
            mRNA_Sequence = '--' + mRNA_Sequence
            mRNA_Sequence = mRNA_Sequence.replace('T','U')
            first_nt = 0
            last_nt = 23
            flag = True
            sequence = []
            deltaGibbs = []
            while(flag == True):
               if (last_nt <= (len(mRNA_Sequence))):
                  sequence.append(mRNA_Sequence[first_nt:last_nt])
                  deltaGibbs.append(deltaGibbCalculator(mRNA_Sequence[first_nt:last_nt]))
               else:
                  flag = False
               first_nt+=1
               last_nt+=1
	        #This part sorts the sequences by gibbs free energy 
            positionOfSeq = []		 
            for x in range(len(sequence)):
               positionOfSeq.append(x)     
            bubblesort(positionOfSeq,sequence, deltaGibbs)   
            sense = []
            antisense = []
            for x in range(len(sequence)):
               sense.append(sequence[x][2:24])
               antisense.append(complementSeqUracil(invertSeq(sequence[x][0:21])))
            now = datetime.datetime.now()		 
            outputName = self.nameLine_19.text()
            if (outputName == ''):
               QMessageBox.about(self, "Warning", "Please enter a Target name")
            else:		 
               saveFile = open(outputName + now.strftime("_siRNA_Gibbs_SSD %Y-%m-%d [%H h %M min].txt"), 'w')
               saveFile.write('********************************************************************************\n')
               saveFile.write('                         siRNA and shRNA designer (SSD)\n')
               saveFile.write('********************************************************************************\n')
               saveFile.write('Version 1.0 - 2018\n')
               saveFile.write('Authors: Gabriel Jose de Carli (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Abdon Troche Rotela (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Danilo Fernandez Rios (National University of Asuncion - Paraguay)\n')
               saveFile.write('         Greice Lubini (University of Sao Paulo - Brazil)\n')
               saveFile.write('         Tiago Campos Pereira  (University of Sao Paulo - Brazil)\n\n')
               saveFile.write('siRNAs below are ordered from the most functional to the least one\n')
               saveFile.write("All sequences are oriented 5' to 3'.\n")
               saveFile.write('********************************************************************************\n')
               saveFile.write('Gene name: '+outputName+'\n')
               saveFile.write('********************************************************************************\n')	  
               for x in range(len(sense)):
                  saveFile.write('POSITION:'+ str(positionOfSeq[x] + 1) +'\n')
                  saveFile.write('\nSense RNA:	')
                  saveFile.write(sense[x])
                  saveFile.write('\nAntisense RNA:	')
                  saveFile.write(antisense[x])
                  saveFile.write('\n\n')
                  saveFile.write('RESULTS:\n')
                  saveFile.write('\nDelta Gibbs: ')
                  formatted = format(deltaGibbs[x])
                  saveFile.write(str(formatted))
                  if(deltaGibbs[x] >= 0):
                     saveFile.write(' => Functional siRNA\n')
                  else:
                     saveFile.write(' => Non-Functional siRNA\n')
                  saveFile.write('-----------------------------------------------------------------------\n\n')			
               saveFile.close()	 
               #Success Message 
               QMessageBox.about(self, "Process(Ordered by Gibbs free energy)", "The Strand Analysis was successfully done")
	    
	  
   def tab5UI(self):
      layout=QFormLayout()
      label_0 = QLabel("Always use CAPITAL LETTERS to insert nucleotde sequences (5' to 3'). Never use sequences containing spaces or numbers. ")
      label_1 = QLabel("The software's input file format is FASTA(.txt).\n\tExample: 'input_file_name.txt' ")
      label_2 = QLabel("The FASTA fle content for shRNA design is:\n----------------------------------------------------------\n\t>sense1\n\tGAGGCAUUCUUACGGUAUCAU\n\t>antisense1\n\tGAUACCGUAAGAAUGCCUCCU\n\t>sense2\n\tCUGCAGUCUUCGUUCUUAUCA\n\t>antisense2\n\tAUAAGAACGAAGACUGCCGUG")
      label_3 = QLabel("----------------------------------------------------------")
      label_01 = QLabel("The FASTA fle content for siRNA design is:\n------------------------------------------------------------------------------------------------------------------------------------\n\t>Gene name 1")
      label_03 = QLabel("\tGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGTTGGCATGGT")
      label_04 = QLabel("\tAGCTGGGATGTTAGGGCTCAGGGAAGAAAAGTCAGAAGACCAGGACCTCCAGGGCCTCAAGGACAAACCC")
      label_05 = QLabel("\tCTCAAGTTTAAAAAGGTGAAGAAAGATAAGAAAGAAGAGAAAGAGGGCAAGCATGAGCCCTCTCAGCCAT")
      label_06 = QLabel("\tCAGCCTCTCACTCTGCTGAGCCCGCAGAGGCAGGCAAAGCAGAGACATCAGAAGGGTCAGGCTCCGCCCC")
      label_07 = QLabel("\tGGCTGTGCCGGAAGCTTCTGCCTCCCCCAAACAGCGGCGCTCCATCATCCGTGACCGGGGACCCATGTAT")
      label_08 = QLabel("\tGATGACCCTCTCCTGCCTGAAGGCTGGACACGGAAGCTTAAGCAAAGGAAATCTGGCCGCTCTGCTGGGA")
      label_09 = QLabel("\tAGTATGATGTGTATTTGATCAATCCCCAGGGAAAAGCCTTTCGCTCTAAAGTGGAGTTGATTGCGTACTT")
      label_010 = QLabel("\tCGAAAAGGTAGGCGACACATCCCTGGACCCTAATGATTTTGACTTCACGGTAACTGGGAGAGGGAGCCCC")
      label_011 = QLabel("\tTCCCGGCGAGAGCAGAAACTCTCTAAGAAGCCCAAATCTCCCAAAGCTCCAGGAACTGGCAGAGGCCGGG")
      label_012 = QLabel("\tGACGCCCCAAAGGGAGCGGCACCACGAGACCCAAGGCGGCCACGTCAGAGGGTGTGCAGGTGAAAAGGGT")
      label_013 = QLabel("\tCCTGGAGAAAAGTCCTGGGAAGCTCCTTGTCAAGATGCCTTTTCAAACTTCGCCAGGGGGCAAGGCTGAG")
      label_014 = QLabel("\tGGGGGTGGGGCTCTCACATCTCTCCAGGTCATGGTGATCAAACGCCCCGGCAGGAAGCGAAAAGCTGAGG")
      label_015 = QLabel("\tCCGACCCTCAGGCCATTCCCAAGAAACGGGGCCGAAAGCCGGGGAGTGTGGTGGCAGCCGCTGCCGCCGA")
      label_02 = QLabel("------------------------------------------------------------------------------------------------------------------------------------\nIt is very important to follow the input structure showed above to avoid any errors")
      label_4 = QLabel("The output file will be a TXT format and will be created in the software's folder.\n\tExample: 'output_file_name_siRNA_Gibbs_SSD_Date_[Hour].txt' for the siRNA designer\n\t\t'output_file_name_Trivalent_shRNA_Date_[Hour].txt' for the shRNA designer")	  
      layout.addRow(label_0)
      layout.addRow(label_1)
      layout.addRow(label_2)
      layout.addRow(label_3)
      layout.addRow(label_01)
      layout.addRow(label_03)
      layout.addRow(label_04)
      layout.addRow(label_05)
      layout.addRow(label_06)
      layout.addRow(label_07)
      layout.addRow(label_08)
      layout.addRow(label_09)
      layout.addRow(label_010)
      layout.addRow(label_011)
      layout.addRow(label_012)
      layout.addRow(label_013)
      layout.addRow(label_014)
      layout.addRow(label_015)	  
      layout.addRow(label_02)	  
      layout.addRow(label_4)	  
      self.tab5.setLayout(layout)
	     
app = QApplication(sys.argv)
screen = Window()
screen.showMaximized()
sys.exit(app.exec_())