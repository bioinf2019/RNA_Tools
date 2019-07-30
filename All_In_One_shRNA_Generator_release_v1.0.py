#!/usr/bin/env python3
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import sys

def format(value):
    return "%.3f" % value


def bubblesort(list_1, list_2):
# Swap the elements to arrange in order
    for iter_num in range(len(list_2)-1,0,-1):
        for idx in range(iter_num):
            if list_2[idx]<list_2[idx+1]:
                temp_1 = list_1[idx]
                temp_2 = list_2[idx]
                list_1[idx] = list_1[idx+1]
                list_2[idx] = list_2[idx+1]
                list_1[idx+1] = temp_1
                list_2[idx+1] = temp_2	
    return list_1, list_2

def deltaGibbCalculator(string):
  #Gibbs Calculator
  initialSubString = string[2:6]
  initialSubString = initialSubString + initialSubString[3]
  finalSubString = string[16:20] 
  finalSubString = finalSubString + finalSubString[3]

  initialGibbs = 0.0
  finalGibbs = 0.0
  deltaGibbs = 0.0

  first_nt = 0
  last_nt = 2
  flag = True
  while(flag == True):
    if (last_nt <= (len(initialSubString))):
      if((initialSubString[first_nt:last_nt] == 'AA') or (initialSubString[first_nt:last_nt] == 'AU') or (initialSubString[first_nt:last_nt] == 'UU')):
        initialGibbs += 1.1
      elif((initialSubString[first_nt:last_nt] == 'AC') or (initialSubString[first_nt:last_nt] == 'GU')):
        initialGibbs += 2.4
      elif((initialSubString[first_nt:last_nt] == 'AG') or (initialSubString[first_nt:last_nt] == 'CU')):
        initialGibbs += 1.9
      elif((initialSubString[first_nt:last_nt] == 'CA') or (initialSubString[first_nt:last_nt] == 'CG') or (initialSubString[first_nt:last_nt] == 'UG')):
        initialGibbs += 2.2
      elif((initialSubString[first_nt:last_nt] == 'CC') or (initialSubString[first_nt:last_nt] == 'GG')):
        initialGibbs += 3.3
      elif((initialSubString[first_nt:last_nt] == 'GA')):
        initialGibbs += 2.7
      elif((initialSubString[first_nt:last_nt] == 'GC')):
        initialGibbs += 3.8
      elif((initialSubString[first_nt:last_nt] == 'UA')):
        initialGibbs += 1.4
      elif((initialSubString[first_nt:last_nt] == 'UC')):
        initialGibbs += 1.6
    else:
      flag = False
    first_nt+=1
    last_nt+=1

  first_nt = 0
  last_nt = 2
  flag = True
  while(flag == True):
    if (last_nt <= (len(finalSubString))):
      if((finalSubString[first_nt:last_nt] == 'AA') or (finalSubString[first_nt:last_nt] == 'AU') or (finalSubString[first_nt:last_nt] == 'UU')):
        finalGibbs += 1.1
      elif((finalSubString[first_nt:last_nt] == 'AC') or (finalSubString[first_nt:last_nt] == 'GU')):
        finalGibbs += 2.4
      elif((finalSubString[first_nt:last_nt] == 'AG') or (finalSubString[first_nt:last_nt] == 'CU')):
        finalGibbs += 1.9
      elif((finalSubString[first_nt:last_nt] == 'CA') or (finalSubString[first_nt:last_nt] == 'CG') or (finalSubString[first_nt:last_nt] == 'UG')):
        finalGibbs += 2.2
      elif((finalSubString[first_nt:last_nt] == 'CC') or (finalSubString[first_nt:last_nt] == 'GG')):
        finalGibbs += 3.3
      elif((finalSubString[first_nt:last_nt] == 'GA')):
        finalGibbs += 2.7
      elif((finalSubString[first_nt:last_nt] == 'GC')):
        finalGibbs += 3.8
      elif((finalSubString[first_nt:last_nt] == 'UA')):
        finalGibbs += 1.4
      elif((finalSubString[first_nt:last_nt] == 'UC')):
        finalGibbs += 1.6
    else:
      flag = False
    first_nt+=1
    last_nt+=1

  deltaGibbs = finalGibbs - initialGibbs
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
  return compSeq  

def multimeric_shRNA(siRNAsSenseStrandFiveThree_21_mer, siRNAsAntisenseStrandFiveThree_21_mer, minorLoop, largerLoop): 
  reverseComplementOfTheLastTwoNucleotides = siRNAsAntisenseStrandFiveThree_21_mer[-2] + siRNAsAntisenseStrandFiveThree_21_mer[-1]
  reverseComplementOfTheLastTwoNucleotides = complementSeqUracil(reverseComplementOfTheLastTwoNucleotides)
  sequence_23_mer = invertSeq(reverseComplementOfTheLastTwoNucleotides) + siRNAsSenseStrandFiveThree_21_mer
  replacedNucleotidesSequence = sequence_23_mer.replace('U','T')
  complementSequence = complementSeq(replacedNucleotidesSequence) 
  invertedComplementSequence = invertSeq(complementSequence)
  completeLinearSequenceGeneratorOf_shRNA = minorLoop + replacedNucleotidesSequence + largerLoop + invertedComplementSequence
  return completeLinearSequenceGeneratorOf_shRNA
  
class Window(QWidget):
   def __init__(self, parent=None):
      super(Window, self).__init__(parent)      
      self.setMinimumSize(1360, 600)
      self.setMaximumSize(1360, 600)
	  
      QWidget.__init__(self)
      layout = QGridLayout()
      self.setLayout(layout)
      self.tab1 = QWidget()
      self.tab2 = QWidget()
      self.tab3 = QWidget()
      self.tab4 = QWidget()
      tabwidget = QTabWidget()
      tabwidget.addTab(self.tab1, "shRNA Monovalent")
      tabwidget.addTab(self.tab2, "shRNA Bivalent")
      tabwidget.addTab(self.tab3, "shRNA Trivalent")
      tabwidget.addTab(self.tab4, "Strand Analysis")
      layout.addWidget(tabwidget, 100, 100)
      self.tab1UI()
      self.tab2UI()
      self.tab3UI()
      self.tab4UI()
      self.setWindowTitle("shRNA Multivalent Generator")
   def tab1UI(self):
      layout=QFormLayout()
      self.nameLine_1 = QLineEdit()
      self.nameLine_2 = QLineEdit()
      self.nameLine_3 = QLineEdit()
      self.nameLine_4 = QLineEdit()
      layout.addRow("siRNA's sense strand:",self.nameLine_1)
      layout.addRow("siRNA's antisense strand:",self.nameLine_2)
      layout.addRow("shRNA Multimeric:",self.nameLine_4)
      self.openFileButton_1 = QPushButton("Open FASTA file")
      self.generateButton_1 = QPushButton("Generate shRNA")	  
      layout.addRow(self.openFileButton_1, self.nameLine_3)
      layout.addRow(self.generateButton_1)
	  #Conditions When buttons are pressed
      self.openFileButton_1.clicked.connect(self.open)
      self.generateButton_1.clicked.connect(self.generateMultimeric)
	  
      self.tab1.setLayout(layout)
	  
   def open(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_3.setText(fileName[0])
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
      siRNAsAntisenseStrandFiveThree_21_mer = self.nameLine_2.text()
      reverseComplementOfTheLastTwoNucleotides = siRNAsAntisenseStrandFiveThree_21_mer[-2] + siRNAsAntisenseStrandFiveThree_21_mer[-1]
      reverseComplementOfTheLastTwoNucleotides = complementSeqUracil(reverseComplementOfTheLastTwoNucleotides)
      sequence_23_mer = invertSeq(reverseComplementOfTheLastTwoNucleotides) + siRNAsSenseStrandFiveThree_21_mer
      replacedNucleotidesSequence = sequence_23_mer.replace('U','T')  
      complementSequence = complementSeq(replacedNucleotidesSequence) 
      invertedComplementSequence = invertSeq(complementSequence) 
      minorLoop = 'TTGTT'
      largerLoop = 'TGACAGGAAG'
      completeLinearSequenceGeneratorOf_shRNA = minorLoop + replacedNucleotidesSequence + largerLoop + invertedComplementSequence
      inicialPortion = completeLinearSequenceGeneratorOf_shRNA[:16]
      finalPortion = completeLinearSequenceGeneratorOf_shRNA[16:]
      shRNA = finalPortion + inicialPortion
      saveFile = open('Multimeric_shRNA.fasta', 'w')
      saveFile.write('>Multimeric_shRNA.fasta of GUI\n')
      saveFile.write(shRNA)
      saveFile.close()
      self.nameLine_4.setText(shRNA)

   def tab2UI(self):
      layout=QFormLayout()
      self.nameLine_5 = QLineEdit()
      self.nameLine_6 = QLineEdit()
      self.nameLine_7 = QLineEdit()
      self.nameLine_8 = QLineEdit()
      self.nameLine_9 = QLineEdit()
      self.nameLine_10 = QLineEdit()	  
      layout.addRow("siRNA's sense strand one:",self.nameLine_5)
      layout.addRow("siRNA's antisense strand one:",self.nameLine_6)
      layout.addRow("siRNA's sense strand two:",self.nameLine_7)
      layout.addRow("siRNA's antisense strand two:",self.nameLine_8)	  
      layout.addRow("shRNA Multimeric:",self.nameLine_10)
      self.openFileButton_2 = QPushButton("Open FASTA file")
      self.generateButton_2 = QPushButton("Generate shRNA")	  
      layout.addRow(self.openFileButton_2, self.nameLine_9)
      layout.addRow(self.generateButton_2)
	  #Conditions When buttons are pressed
      self.openFileButton_2.clicked.connect(self.open2)
      self.generateButton_2.clicked.connect(self.generateMultimeric2)
	  
      self.tab2.setLayout(layout)	  
	  
   def open2(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_9.setText(fileName[0])
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
      siRNAsAntisenseStrandFiveThree_21_merNumberOne = self.nameLine_6.text()
      siRNAsSenseStrandFiveThree_21_merNumberTwo = self.nameLine_7.text()
      siRNAsAntisenseStrandFiveThree_21_merNumberTwo = self.nameLine_8.text()
      minorLoopNumberOne = 'TTGTT'
      minorLoopNumberTwo = 'CCACC'
      largerLoopNumberAll = 'TGACAGGAAG'
      completeLinearSequence = multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberOne, siRNAsAntisenseStrandFiveThree_21_merNumberOne, minorLoopNumberOne, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberTwo, siRNAsAntisenseStrandFiveThree_21_merNumberTwo, minorLoopNumberTwo, largerLoopNumberAll)
      inicialPortion = completeLinearSequence[:16]
      finalPortion = completeLinearSequence[16:]
      BivalentMultimeric_shRNA = finalPortion + inicialPortion
      saveFile = open('BivalentMultimeric_shRNA_of GUI.fasta', 'w')
      saveFile.write('>BivalentMultimeric_shRNA_of GUI\n')
      saveFile.write(BivalentMultimeric_shRNA)
      saveFile.close()
      self.nameLine_10.setText(BivalentMultimeric_shRNA)
	  
   def tab3UI(self):
      layout=QFormLayout()
      self.nameLine_11 = QLineEdit()
      self.nameLine_12 = QLineEdit()
      self.nameLine_13 = QLineEdit()
      self.nameLine_14 = QLineEdit()
      self.nameLine_15 = QLineEdit()
      self.nameLine_16 = QLineEdit()	  
      self.nameLine_17 = QLineEdit()
      self.nameLine_18 = QLineEdit()	  
      layout.addRow("siRNA's sense strand one:",self.nameLine_11)
      layout.addRow("siRNA's antisense strand one:",self.nameLine_12)
      layout.addRow("siRNA's sense strand two:",self.nameLine_13)
      layout.addRow("siRNA's antisense strand two:",self.nameLine_14)	  
      layout.addRow("siRNA's sense strand three:",self.nameLine_15)
      layout.addRow("siRNA's antisense strand three:",self.nameLine_16)	  
      layout.addRow("shRNA Multimeric:",self.nameLine_18)
      self.openFileButton_3 = QPushButton("Open FASTA file")
      self.generateButton_3 = QPushButton("Generate shRNA")	  
      layout.addRow(self.openFileButton_3, self.nameLine_17)
      layout.addRow(self.generateButton_3)
	  #Conditions When buttons are pressed
      self.openFileButton_3.clicked.connect(self.open3)
      self.generateButton_3.clicked.connect(self.generateMultimeric3)

      self.tab3.setLayout(layout)

   def open3(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_17.setText(fileName[0])
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
      siRNAsAntisenseStrandFiveThree_21_merNumberOne = self.nameLine_12.text()
      siRNAsSenseStrandFiveThree_21_merNumberTwo = self.nameLine_13.text()
      siRNAsAntisenseStrandFiveThree_21_merNumberTwo = self.nameLine_14.text()
      siRNAsSenseStrandFiveThree_21_merNumberThree = self.nameLine_15.text()
      siRNAsAntisenseStrandFiveThree_21_merNumberThree = self.nameLine_16.text()
		
      minorLoopNumberOne = 'TTGTT'
      minorLoopNumberTwo = 'CCACC'
      minorLoopNumberThree = 'CGTGC'
      largerLoopNumberAll = 'TGACAGGAAG'
      completeLinearSequence = multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberOne, siRNAsAntisenseStrandFiveThree_21_merNumberOne, minorLoopNumberOne, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberTwo, siRNAsAntisenseStrandFiveThree_21_merNumberTwo, minorLoopNumberTwo, largerLoopNumberAll) + multimeric_shRNA(siRNAsSenseStrandFiveThree_21_merNumberThree, siRNAsAntisenseStrandFiveThree_21_merNumberThree, minorLoopNumberThree, largerLoopNumberAll)
      inicialPortion = completeLinearSequence[:16]
      finalPortion = completeLinearSequence[16:]
      TrivalentMultimeric_shRNA = finalPortion + inicialPortion
      saveFile = open('TrivalentMultimeric_shRNA_of GUI.fasta', 'w')
      saveFile.write('>TrivalentMultimeric_shRNA_of GUI\n')	  
      saveFile.write(TrivalentMultimeric_shRNA)
      saveFile.close()
      self.nameLine_18.setText(TrivalentMultimeric_shRNA)

   def tab4UI(self):
      layout=QFormLayout()
      self.nameLine_19 = QLineEdit()	
      self.nameLine_20 = QLineEdit()
      self.nameLine_21 = QLineEdit()	  
      layout.addRow("Gene name:",self.nameLine_19)
      layout.addRow("Gene:",self.nameLine_20)	  
      self.openFileButton_4 = QPushButton("Open FASTA file")
      self.generateButton_4 = QPushButton("Generate siRNAs")	  
      layout.addRow(self.openFileButton_4, self.nameLine_21)
      layout.addRow(self.generateButton_4)
	  #Conditions When buttons are pressed
      self.openFileButton_4.clicked.connect(self.open4)
      self.generateButton_4.clicked.connect(self.generate_siRNAs)

      self.tab4.setLayout(layout)

   def open4(self):
      fileName = QFileDialog.getOpenFileName(self, 'OpenFile')
      self.nameLine_21.setText(fileName[0])
      with open(fileName[0]) as fp:
         nameList = []
         seqList = []
         for name, seq in read_fasta(fp):
            nameList.append(name)
            seqList.append(seq)
      seq_1 = seqList[0]	  
      self.nameLine_20.setText(seq_1)	  

   def generate_siRNAs(self):
      miRNA_Sequence = self.nameLine_20.text()
      miRNA_Sequence = miRNA_Sequence.replace('T','U')
      first_nt = 0
      last_nt = 23
      flag = True
      sequence = []
      deltaGibbs = []
      while(flag == True):
         if (last_nt <= (len(miRNA_Sequence))):
            sequence.append(miRNA_Sequence[first_nt:last_nt])
            deltaGibbs.append(deltaGibbCalculator(miRNA_Sequence[first_nt:last_nt]))
         else:
            flag = False
         first_nt+=1
         last_nt+=1

      bubblesort(sequence, deltaGibbs)

      sense = []
      antisense = []

      for x in range(len(sequence)):
         sense.append(sequence[x][2:24])
         antisense.append(complementSeqUracil(invertSeq(sequence[x][0:21])))	
      outputName = self.nameLine_19.text()
      saveFile = open(outputName + '.txt', 'w')
      saveFile.write('Gene name:'+outputName+'\n')
      for x in range(len(sense)):
         saveFile.write('siRNA:\n')
         saveFile.write('Caracteristic:\n')
         if(deltaGibbs[x] >= 0):
            saveFile.write('Functional siRNA:\n')
         else:
            saveFile.write('Non-Functional siRNA:\n')
         if((sense[x].find('AAAA') != -1) or (sense[x].find('UUUU') != -1) or (sense[x].find('GGGG') != -1) or (sense[x].find('CCCC') != -1)):
            saveFile.write('4 or more idential nuileotdes in a row - inadequate for use:\n')
         saveFile.write('Sense strand: ')
         saveFile.write(sense[x])
         saveFile.write('\nAntisense strand: ')
         saveFile.write(antisense[x])
         saveFile.write('\nTarget: ')
         saveFile.write(sequence[x])
         saveFile.write('\nDelta Gibbs: ')
         formatted = format(deltaGibbs[x])
         saveFile.write(str(formatted))
         saveFile.write('\n\n')
      saveFile.close()  
	  
app = QApplication(sys.argv)
screen = Window()
screen.show()
sys.exit(app.exec_())