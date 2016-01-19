'''
Created on Jan 19, 2016
Link : http://manta.bmrb.wisc.edu/api/a/chemshifts.php?idlist=1234,5678,90
idlist can be list of ids
@author: kumaran
'''

import requests,re

class PyBMRB(object):
    '''
    API to fetch chemical shift data from BMRB database
    
    '''


    def __init__(self, bmrbids):
        '''
        Constructor
        '''
        self.api_url="http://manta.bmrb.wisc.edu/api/a/chemshifts.php?idlist="+bmrbids
        


    def get_data(self):
        self.resp=requests.get(self.api_url)
        
        if self.resp.status_code != 200:
            # This means something went wrong
            raise ApiError(self.resp.status_code)
        
        self.data=self.parse_data(self.resp.content)
        
        
    def parse_data(self,data):
        tmp_dat=re.findall(r'\[(\S+)\]', data)
        tmp_dat2=tmp_dat[0].split("],[")
        tmp_dat.remove(tmp_dat[0])
        tmp_dat.reverse()
        tmp_dat.append(tmp_dat2[1])
        tmp_dat.append(tmp_dat2[0])
        tmp_dat.reverse()
        return tmp_dat
        


class ApiError(Exception):
    
    def __init__(self, errvalue):
        self.errvalue=errvalue
        if self.errvalue==404:
            self.message="Server Not Found"
        else:
            self.message="Unknown Error"
    def __str__(self):
        return repr(self.message)

if __name__=="__main__":
    b=PyBMRB('1234,90')
    b.get_data()
    for i in b.data:
        print i