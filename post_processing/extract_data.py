from numpy import *
from scipy import *


#This is a python module used to load data from file
#reject_el represent all marked character to reject a line (default to none).
#caution: this module can only load numeric data !

def read_file(filename,reject_el="void",separator="space"):
    
    dummy=open(filename,'r')
    
    #count line
    if reject_el=="void":
        reject_list=0
    else:
        reject_symbol=reject_el.split()
        reject_list=len(reject_symbol)
    
    #print reject_symbol
    row_count=0 #count datalength
    marker=0
    str_array=array([])
    
    for line in dummy:
        #print line.find(reject_symbol[0],0,100)
        for i in range(0,reject_list):
            if (line.find(reject_symbol[i],0,100)!=-1):
                marker=1
                break
        
        if (marker==0):
            row_count=row_count+1
            str_array=append(str_array,line)
        else:
            marker=0
            
    if (separator=="space"):
        col_count=len(str_array[0].split())  
    else:
        col_count=len(str_array[0].split(separator))
        
    data_matrix=ndarray(shape=(row_count,col_count))
    
    for i in range(0,row_count):
        
        if (separator=="space"):
            dummy_data=str_array[i].split()  
        else:
            dummy_data=str_array[i].split(separator)
        
        for k in range(0,col_count):
						#print (dummy_data)
						data_matrix[i][k]=float(dummy_data[k])
    
    dummy.close()
    
    return data_matrix
        
def save_ndarray(output,data,sep=",",form='%d'):
    
    nrow=data.shape[0]
    ncol=data.shape[1]
    dummy=open(output,'w')
    
    for i in range(0,nrow):
        for j in range(0,ncol):
            if (j!=(ncol-1)):
                formatting=form+','
                dummy.write(formatting % (data[i,j]))
            else:
                formatting=form+'\n'
                dummy.write(formatting % (data[i,j]))
    
    dummy.close()
        
          
        