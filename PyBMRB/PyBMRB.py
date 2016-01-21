'''
Created on Jan 19, 2016
Link : http://manta.bmrb.wisc.edu/api/a/chemshifts.php?idlist=1234,5678,90
idlist can be list of ids
@author: kumaran
'''

import requests,re
from operator import itemgetter
from string import atoi
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
        tmp_dat3=[i.split(",") for i in tmp_dat]
        for i in tmp_dat3:
            for j in range(len(i)):
                i[j]=i[j].replace("\"","")
                try:
                    i[j]=atoi(i[j])
                except ValueError:
                    pass
        tmp_dat4=sorted(tmp_dat3, key=itemgetter(0,10,3,1))
        tmp_dat5=[[str(j) for j in i] for i in tmp_dat4]
        return tmp_dat5
    
    def header(self):
        return self.data[-1]
    
    def value(self):
        return self.data[:-1]
    
    def writeFile(self,path):
        bmrbid=self.value()[0][0].replace("\"","")
        f=open('%s/bmrb%s.txt'%(path,bmrbid),'w')
        f.write("#%s\n"%("\t".join(self.header())))
        for i in self.value():
            if i[0].replace("\"","")==bmrbid:
                f.write("%s\n"%("\t".join(i)))
            else:
                f.close()
                bmrbid=i[0].replace("\"","")
                f=open('%s/bmrb%s.txt'%(path,bmrbid),'w')
                f.write("#%s\n"%("\t".join(self.header())))
                f.write("%s\n"%("\t".join(i)))
        f.close()
        



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
    b=PyBMRB('90,18857')
    b.get_data()
    b.writeFile('/home/kumaran/Desktop')