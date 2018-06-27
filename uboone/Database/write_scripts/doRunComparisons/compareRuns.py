import numpy as np
import csv
import sys

def declareDataStruct():
  dataStruct = {}
  for i in range(8256):
    dataStruct[i] = {}
  return dataStruct

def extractToDataStruct(inFile, runNumber, nColumns, dataStruct):
  for i, line in enumerate(inFile):
    if i<4: continue
    elements = line.rstrip("\n").split(',')
    dataStruct[int(elements[0])][runNumber] = elements

def getHeaders(inFile, hStruct):
  for i, line in enumerate(inFile):
    if i==2:
      elements = line.rstrip("\n").split(',')
      for i in range(0, len(elements)):
        hStruct.append(elements[i])


with open('electronicsDbFiles.list', 'r') as elecFileList, open('chanStatDbFiles.list', 'r') as chanStatFileList: 

  elecDataStruct = declareDataStruct()
  chanStatDataStruct = declareDataStruct()
  
  elecHeaders = []
  chanStatHeaders = []
  runNumbers = []

  list = [[elecDataStruct,elecFileList,6,elecHeaders],[chanStatDataStruct,chanStatFileList,2,chanStatHeaders]]
  

  j=0
  for struct,fileList,nColumns,headerStruct in list:
    for i,file in enumerate(fileList):
      inFile = open(file.rstrip("\n"))
      if i == 0:
        getHeaders(inFile, headerStruct)
      inFile = open(file.rstrip("\n"))
      runNumber = int(file.split('_')[1].rstrip(".csv\n"))
      if j == 0:
        runNumbers.append(runNumber)
      extractToDataStruct(inFile, runNumber, nColumns, struct)
    j+=1
  
  origStdOut = sys.stdout
  outputFile=open('./dbStatusChanges.txt', 'w+')
  sys.stdout = outputFile
  
  m = 0
  for index,i in enumerate(runNumbers):
    if m == 0:
      m+=1
      continue
    print "------- Run ", i,"------- " 
    for j in range (0, 8256):
      for k in range(1, len(elecHeaders)):
        if elecDataStruct[j][runNumbers[index-1]][k] != elecDataStruct[j][runNumbers[index]][k]:
          print "|--CHANNEL ",j,": ",elecHeaders[k],": ",elecDataStruct[j][runNumbers[index-1]][k],">>",elecDataStruct[j][runNumbers[index]][k]
      for l in range(1,len(chanStatHeaders)):
        if chanStatDataStruct[j][runNumbers[index-1]][l] != chanStatDataStruct[j][runNumbers[index]][l]:
          print "|--CHANNEL ",j,": ",chanStatHeaders[l],": ", chanStatDataStruct[j][runNumbers[index-1]][l],">>", chanStatDataStruct[j][runNumbers[index]][l]
  sys.stdout = origStdOut
  outputFile.close()
