'''
Created on Jan 19, 2016
Link : http://manta.bmrb.wisc.edu/api/a/chemshifts.php?idlist=1234,5678,90
idlist can be list of ids
@author: kumaran
'''


import pynmrstar 

class PyBMRB(object):
    '''
    Plotly tools for BMRB
    
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def fetch_entry_data(self,entryid):
        try:
            dat = pynmrstar.Entry.from_database(entryid)
            print dat
        except IOError as e:
            print e


   

if __name__=="__main__":
    p = PyBMRB()
    p.fetch_entry_data('15060x')