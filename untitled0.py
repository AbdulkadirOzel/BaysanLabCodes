# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 00:16:03 2022

@author: Abdulkadir
"""

import  jpype
import  asposecells     
jpype.startJVM() 
from asposecells.api import Workbook
workbook = Workbook("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastp//SAA11A2//fastp.json")
workbook.Save("C://Users//Abdulkadir//Desktop//Genome_Project_Documents//haplo//fastp//SAA11A2//Output.txt")
jpype.shutdownJVM()