# Load Packages

import numpy
import csv

from class_definition import *

#####################################################################################################################
def initializer(self, x):
    self.x = x
#####################################################################################################################
def create_classes(volume_='',classDict_={}):
    """
    Add the classes described in the dictionary 'classDict_' to the object 'volume_'.
    
    Arguments:
    ----------
    
    Returns:
    --------

    """
    print('aqui')
    classDict_['__init__'] = initializer
    volume_ = type('field', (object,), classDict_)
    print("clases creadas...")
    return volume_
#####################################################################################################################
