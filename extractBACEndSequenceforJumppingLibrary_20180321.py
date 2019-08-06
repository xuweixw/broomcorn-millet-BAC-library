
# coding: utf-8

# ## 提取pacbio数据双端序列
# 
# #### 利用三代长读段高通量获取BAC BES序列的技术
# #### BES sequencing with Pacbio long-read technology
# 
# constructing sequence == ConstructingSeq

# In[1]:


import pandas as pd


# ## 创建ReadAnnotation
# #### 传入单条读段的blast注释数据框和读段长度
# > DataFrame 应仅包含3列，且列名分别为 ['tag','start','end']

# In[2]:




class ReadAnnotation:
    def __init__(self,anno=None,length=None):
        self.raw_data = anno
        self.read_length = length
        self.annotated_data = None
        
        '''
        A real read must is such: 
        len(Amp) ----- (1,2,3,4...n)
        len(R1) or len(R2) ----- (0,1)
        '''
        self.isException = False
        
        self.annotate()
            
    def annotate(self):
        a = self.raw_data.sort_values(by='start')
        amp = a.loc[a.tag=='Amp',:].sort_values(by='start')
        R1 = a.loc[a.tag=='R1',:]
        R2 = a.loc[a.tag=='R2',:]
        
        if len(amp) > 0 :     # check amp
            BES_R_end = int(amp[0:1].start) - 1
            BES_F_start = int(amp[-1:].end) + 1
            
            if len(R1) == 0 and len(R2) == 0:    #  R1 and R2 have not
                BES_R_start = 1
                BES_F_end = self.read_length
            elif len(R1) == 1 and len(R2) == 0:  # R1 has, R2 not
                if int(R1.start) < BES_R_end:         # R1 in left
                    BES_R_start = int(R1.end) +1
                    BES_F_end = self.read_length
                elif int(R1.end) > BES_F_start:       # R1 in right
                    BES_R_start = 1
                    BES_F_end = int(R1.start) -1
                else:
                    self.isException = True
            elif len(R1) == 0 and len(R2) == 1:  # R1 not, R2 have
                if int(R2.start) < BES_R_end:         # R2 in left
                    BES_R_start = int(R2.end) + 1
                    BES_F_end = self.read_length
                elif int(R2.end) > BES_F_start:       # R2 in right
                    BES_R_start = 1
                    BES_F_end = int(R2.start) -1
                else:
                    self.isException = True
            elif len(R1) == 1 and len(R2) == 1:  # R1 and R2 have
                if int(R1.start) < BES_R_end:         # R1 in the front of R2
                    BES_R_start = int(R1.end) + 1
                    BES_F_end = int(R2.start) - 1
                elif int(R2.start) < BES_R_end:       # R2 in thr front of R1
                    BES_R_start = int(R2.end) + 1
                    BES_F_end = int(R1.start) - 1
                else:
                    self.isException = True
            else:
                self.isException = True
        else:
            self.isException = True
        
        if not self.isException:
            BES_R = {'tag':'BES_R','start':BES_R_start,'end':BES_R_end}
            BES_F = {'tag':'BES_F','start':BES_F_start,'end':BES_F_end}
            b = a.append([BES_R,BES_F])
            self.annotated_data = b.sort_values(by='start')


# ## 创建ConstractingSeq类
# 
# ### 类的调用
# > 构建类实例时需传入seq参数，seq必须是SeqRecord对象
# 
# >随后调用ConstractingSeq.annotate()方法，传入参数，内部调用ReadAnnotation类进行注释，结果记录在ConstractingSeq.annotation属性中。

# In[3]:



class ConstructingSeq:
    def __init__(self,seq,min_BES_len=50 ): # seq is a SeqRecord class
        self.id = None
        self.len = len(seq)
        self.seq = seq
        self.min_BES_len = min_BES_len
        self.annotation = None
        self.BES_F = None
        self.BES_R = None
        self.isException = False
        
    def check_type(self):
        if ( isinstance(self.len,int) and 
             isinstance(self.annotation, pd.core.frame.DataFrame) and 
             isinstance(self.BES_F, list) and 
             isinstance(self.BES_R, list) and 
             isinstance(self.isException, bool) ) :
            return True
        else:
            return False
        
    def addSeq(self,seq):
        self.seq = seq
        
    def annotate(self,anno_DataFrame):
        RA = ReadAnnotation(anno_DataFrame.loc[:,['tag','start','end']],len(self.seq))
        self.annotation = RA.annotated_data
        if isinstance(self.annotation,pd.core.frame.DataFrame):
            self.BES_F = [int(self.annotation.loc[self.annotation.tag=='BES_F','start']),
                          int(self.annotation.loc[self.annotation.tag=='BES_F','end'])]
            self.BES_R = [int(self.annotation.loc[self.annotation.tag=='BES_R','start']),
                          int(self.annotation.loc[self.annotation.tag=='BES_R','end'])]
        self.isException = RA.isException
                
        
    def extractPairedSeq(self):
        '''
        BES_F_record, BES_R_record
        '''
        
        if  self.check_type() and self.isPaired() :
            BES_F_record = self.seq[self.BES_F[0]-1:self.BES_F[-1]]
            BES_R_record = self.seq[self.BES_R[0]-1:self.BES_R[-1]]
            return (BES_F_record, BES_R_record)
        else:
            return ('', '')
    
    def extractSingleSeq(self):
        if self.isSingle() :
            if self.BES_F[-1] - self.BES_F[0] > self.min_BES_len :
                return self.seq[self.BES_F[0]-1:self.BES_F[-1]]
            elif self.BES_F[-1] - self.BES_F[0] < self.min_BES_len :
                return self.seq[self.BES_R[0]-1:self.BES_R[-1]]
            else:
                return ''
        else:
            return ''
    
    def isPaired(self):
        if self.check_type() :
            if self.BES_F[-1] - self.BES_F[0] > self.min_BES_len and self.BES_R[-1] -self.BES_R[0] > self.min_BES_len:
                return True
            else:
                return False
        else:
            False
    
    def isSingle(self):
        if self.check_type() :
            if self.BES_F[-1] - self.BES_F[0] < self.min_BES_len and self.BES_R[-1] -self.BES_R[0] > self.min_BES_len:
                return True
            elif self.BES_F[-1] - self.BES_F[0] > self.min_BES_len and self.BES_R[-1] -self.BES_R[0] < self.min_BES_len:
                return True
            else:
                return False
        else:
            return False
        


# ## 编写主程序，调用上面定义的类

# In[4]:


from Bio.SeqIO import parse,write


fasta_path = '/home/wxu/workspace/DZZ/data/m54191_171219_101619.ccs.fasta'
blastn_annotation_path = '/home/wxu/workspace/DZZ/data/result.blast'

BES_F_records = []
BES_R_records = []

BES_Single_records = []

unknow_records = []

## reads counts
unknown_count = 0
pairedBES_count =0
singleBES_count = 0

# read loc file, convert to dict class
loc_dict = {}
loc = pd.read_table(blastn_annotation_path,names=['query_id','query_len','tag','tag_len','start','end','identity'])
for i in set(loc.query_id):
    loc_dict[i] = loc.loc[loc.query_id == i,]
    
    
# read sequence file
for seq_record in parse(fasta_path,'fasta'):
    seq_id = seq_record.id
    
    try:    # if sequence id is not in dict, function will skip this id record.
        seq_anno = loc_dict[seq_id]
        CS = ConstructingSeq(seq_record)
        CS.annotate(seq_anno)
        if not CS.isException:
            if CS.isPaired():     # this sequence include BES two sequences. 
                BES_F_record,BES_R_record = CS.extractPairedSeq()
                if not isinstance(BES_F_record,str): 
                    BES_F_record.description = 'BES_R1'
                    BES_F_records.append(BES_F_record)

                    BES_R_record.description = 'BES_R2'
                    BES_R_records.append(BES_R_record)
                    pairedBES_count += 1
            elif CS.isSingle():   # this sequence only include one BES sequence.
                BES_Single_record = CS.extractSingleSeq()
                if not isinstance(BES_Single_record,str):
                    BES_Single_record.description = 'BES_Single'
                    BES_Single_records.append(BES_Single_record)
                    singleBES_count += 1
            else:
                unknow_records.append(seq_record)
                unknown_count += 1
        else:
            unknow_records.append(seq_record)
            unknown_count += 1
            
    except KeyError:
        unknow_records.append(seq_record)
        unknown_count += 1
        next
        
write(BES_F_records,'BES_R1.fasta','fasta')
write(BES_R_records,'BES_R2.fasta','fasta')
write(BES_Single_records,'BES_Single.fasta','fasta')
write(unknow_records,'unknown_pacbio_reads.fasta','fasta')


with open('run-log','w') as log:
    log.write('Paired-end BES sequences have ' + str(pairedBES_count) + ' pairs.\n')
    log.write('Single-end BES sequences have ' + str(singleBES_count) + ' single.\n')
    log.write(str(unknown_count) + ' reads do not be ananlysed.')
    


# ## Run exception
# 
# ```
# [wxu@comput11 extractBES]$ more nohup.out 
# Traceback (most recent call last):
#   File "./extractBACEndSequence_20180321.py", line 203, in <module>
#     CS.annotate(seq_anno)
#   File "./extractBACEndSequence_20180321.py", line 118, in annotate
#     RA = ReadAnnotation(anno_DataFrame.loc[:,['tag','start','end']],len(self.seq))
#   File "./extractBACEndSequence_20180321.py", line 39, in __init__
#     self.annotate()
#   File "./extractBACEndSequence_20180321.py", line 86, in annotate
#     BES_R = {'tag':'BES_R','start':BES_R_start,'end':BES_R_end}
# UnboundLocalError: local variable 'BES_R_start' referenced before assignment
# [wxu@comput11 extractBES]$ vi extractBACEndSequence_20180321.py 
# 
# ```
