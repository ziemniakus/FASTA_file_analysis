from Bio import SeqIO
import pandas as pd
import numpy as np

try:
    f = open("dna2.fasta")
except IOError:
    print("File dna.example.fasta does not exist!!")
    
seqs={}
for line in f:
# let's discard the newline at the end (if any)
    line=line.rstrip()
# distinguish header from sequence
    if line[0]=='>': # or line.startswith('>')
        words=line.split()
        name=words[0][1:]
        seqs[name]=''
    else : # sequence, not header
        seqs[name] = seqs[name] + line

count=0
lengths_sequences=[]
dictionary_lengths={}

for name,seq in seqs.items():
    #print(name)
    # (1) How many records are in the file?
    count = count + 1  
    # (2) What are the lengths of the sequences in the file? What is the longest sequence and what is the shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers? 
    length=len(seq)
    lengths_sequences.append(length)
    dictionary_lengths[name]=length
    
max_length = max(lengths_sequences)
min_length = min(lengths_sequences)

#print(dictionary_lengths.values())
    
print ("The sequences number is: %d" %count) 
print ("The length of all sequences are:")
print(lengths_sequences)

value_max = [key for key, value in dictionary_lengths.items() if value == max_length]
print("Max value %d for key:" %max_length)
print(value_max)

value_min = [key for key, value in dictionary_lengths.items() if value == min_length]
print("Min value %d for key:" %min_length)
print(value_min)

# (3) In molecular biology, a reading frame is a way of dividing the DNA sequence of nucleotides into a set of consecutive, non-overlapping triplets (or codons). Depending on where we start, there are six possible reading frames: three in the forward (5' to 3') direction and three in the reverse (3' to 5'). For instance, the three possible forward reading frames for the sequence AGGTGACACCGCAAGCCTTATATTAGC are: 

#AGG TGA CAC CGC AAG CCT TAT ATT AGC

#A GGT GAC ACC GCA AGC CTT ATA TTA GC

#AG GTG ACA CCG CAA GCC TTA TAT TAG C 

#These are called reading frames 1, 2, and 3 respectively. An open reading frame (ORF) is the part of a reading frame that has the potential to encode a protein. It starts with a start codon (ATG), and ends with a stop codon (TAA, TAG or TGA). For instance, ATGAAATAG is an ORF of length 9.

#Given an input reading frame on the forward strand (1, 2, or 3) your program should be able to identify all ORFs present in each sequence of the FASTA file, and answer the following questions: what is the length of the longest ORF in the file? What is the identifier of the sequence containing the longest ORF? For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? What is the starting position of the longest ORF in the sequence that contains it? The position should indicate the character number in the sequence. For instance, the following ORF in reading frame 1:

#>sequence1

#ATGCCCTAG

#starts at position 1.

#Note that because the following sequence:

#>sequence2

#ATGAAAAAA

#does not have any stop codon in reading frame 1, we do not consider it to be an ORF in reading frame 1. 

def orfFINDER(dna,frame):
    
    stop_codons = ['tga', 'tag', 'taa']
    start_codon = ['atg']
    start_positions = []
    stop_positions = []
    num_starts=0
    num_stops=0
    
    for i in range(frame,len(dna),3):
        codon=dna[i:i+3].lower()
        if codon in start_codon:
            start_positions += str(i+1).splitlines()
        if codon in stop_codons:
            stop_positions += str(i+1).splitlines()
    
    for line in stop_positions:
        num_stops += 1
    
    for line in start_positions:
        num_starts += 1
    
    orffound = {}
    
    if num_stops >=1 and num_starts >=1: #first statment: the number of stop codons and start condos are greater than or equal to 1;
        
        orfs = True
        stop_before = 0
        start_before = 0
        
        if num_starts > num_stops:
            num_runs = num_starts
        if num_stops > num_starts:
            num_runs = num_stops
        if num_starts == num_stops:
            num_runs = num_starts
            
        position_stop_previous = 0
        position_start_previous = 0
        counter = 0
        
        for position_stop in stop_positions:
            
            position_stop = int(position_stop.rstrip()) + 2
            
            for position_start in start_positions:
                
                position_start = position_start.rstrip()
                
                if int(position_start) < int(position_stop) and int(position_stop) > int(position_stop_previous) and int(position_start) > int(position_stop_previous):
                    
                    counter += 1
                    nameorf = "orf"+str(counter)
                    position_stop_previous += int(position_stop) - int(position_stop_previous)
                    position_start_previous += int(position_start) - int(position_start_previous)
                    sizeorf = int(position_stop) - int(position_start) + 1
                    
                    orffound[nameorf] = position_start,position_stop,sizeorf,frame
                        
                else:
                    pass
    else:  
        orfs = False
    return orffound
    
def add_values_in_dict(sample_dict, key, list_of_values):
    ''' Append multiple values to a key in 
        the given dictionary '''
    if key not in sample_dict:
        sample_dict[key] = list()
    sample_dict[key].append(list_of_values)
    return sample_dict

#DEFINE FRAME TO FIND ORF
#if frame = 0, start from the first position in the sequence
frame=0
len_list=[]
num_list=[]
start_list=[]
stop_list=[]
frame_list=[]
header_list=[]
#EXECUTE THE ORFFINDER FUNCTION
for i in seqs.items():
    header= i[0]
    seq = i[1]
    orf = orfFINDER(seq,frame)
     
    for i in orf.items():
        
        numorf=i[0]
        startorf=orf[numorf][0]
        stoporf=orf[numorf][1]
        lengthorf=orf[numorf][2]
        frameorf=orf[numorf][3]
        #print(header,numorf,"start",startorf,"stop",stoporf,"length",lengthorf,"frame",frameorf)
        header_list.append(header)
        num_list.append(numorf)
        start_list.append(startorf)
        stop_list.append(stoporf)
        len_list.append(lengthorf)
        frame_list.append(frameorf)

#zipped = list(zip(header_list,num_list, start_list, stop_list, len_list, frame_list))
#df = pd.DataFrame(zipped, columns=['Sequence','Number ORF','Start ORF','Stop ORF','Length ORF','Frame ORF'])

frame=1
#EXECUTE THE ORFFINDER FUNCTION
for i in seqs.items():
    header= i[0]
    seq = i[1]
    orf = orfFINDER(seq,frame)
     
    for i in orf.items():
        
        numorf=i[0]
        startorf=orf[numorf][0]
        stoporf=orf[numorf][1]
        lengthorf=orf[numorf][2]
        frameorf=orf[numorf][3]
        #print(header,numorf,"start",startorf,"stop",stoporf,"length",lengthorf,"frame",frameorf)
        header_list.append(header)
        num_list.append(numorf)
        start_list.append(startorf)
        stop_list.append(stoporf)
        len_list.append(lengthorf)
        frame_list.append(frameorf)

#zipped = list(zip(header_list,num_list, start_list, stop_list, len_list, frame_list))
#df1 = pd.DataFrame(zipped, columns=['Sequence','Number ORF','Start ORF','Stop ORF','Length ORF','Frame ORF'])

frame=2

#EXECUTE THE ORFFINDER FUNCTION
for i in seqs.items():
    header= i[0]
    seq = i[1]
    orf = orfFINDER(seq,frame)
     
    for i in orf.items():
        
        numorf=i[0]
        startorf=orf[numorf][0]
        stoporf=orf[numorf][1]
        lengthorf=orf[numorf][2]
        frameorf=orf[numorf][3]
        #print(header,numorf,"start",startorf,"stop",stoporf,"length",lengthorf,"frame",frameorf)
        header_list.append(header)
        num_list.append(numorf)
        start_list.append(startorf)
        stop_list.append(stoporf)
        len_list.append(lengthorf)
        frame_list.append(frameorf)

zipped = list(zip(header_list,num_list, start_list, stop_list, len_list, frame_list))
df = pd.DataFrame(zipped, columns=['Sequence','Number ORF','Start ORF','Stop ORF','Length ORF','Frame ORF'])

#print(df['Length ORF'].tolist())
max_length = max(df['Length ORF'].tolist())
print("Max value %d for sequence:" %max_length)

# MAX LENGTH OF ALL FRAMES AND SEQUENCES:
print(df.loc[df['Length ORF'] == max(df['Length ORF'].tolist())])

# MAX LENGTH OF ORF FOR EACH SEQUENCE IDENTIFIER:
print("MAX LENGTH OF ORF FOR EACH SEQUENCE IDENTIFIER:")
df2=df.groupby('Sequence')['Length ORF'].max().reset_index()
df3 = pd.merge(df.drop_duplicates(), df2.drop_duplicates(), on=['Sequence','Length ORF'],how='inner')
print(df3)

# (4) A repeat is a substring of a DNA sequence that occurs in multiple copies (more than one) somewhere in the sequence. Although repeats can occur on both the forward and reverse strands of the DNA sequence, we will only consider repeats on the forward strand here. Also we will allow repeats to overlap themselves. For example, the sequence ACACA contains two copies of the sequence ACA - once at position 1 (index 0 in Python), and once at position 3. Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA file. Your program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.
		
# This code is contributed by Mohit kumar 29	

def repeated_substring(dna,n):
    prefix_array=[]
    start_position=[]
    stop_position=[]
    dic={}
    for i in range(len(dna)):
        for start in range(len(dna)):
            if dna[start:i]!='' and len(dna[start:i])==n:
                prefix_array.append(dna[start:i])
    for j in range(len(prefix_array)):
        count=0
        if prefix_array[j] in dna:
            substring = prefix_array[j]
        for position in range(len(dna)):
            if dna[position:].startswith(prefix_array[j]):
                count += 1 
                
                #stop_position = start_position+n
                #print(start_position) 
                
        dic[substring] = [count]
            
    return dic

#EXECUTE THE ORFFINDER FUNCTION
substring_list=[]
count_list=[]
header2_list=[]
seq_list=[]
n=7

for j in seqs.items():
    header2= j[0]
    #print(header2)
    seq2 = j[1]
    #print(seq2)
    repeats = repeated_substring(seq2,n)

    for j in repeats.items():
        
        substring=j[0]
        #startorf=repeats[substring][0]
        #stoporf=repeats[substring][1]
        count=repeats[substring][0]
        #print(header2, substring,"count",count)
        header2_list.append(header2)
        #seq_list.append(seq2)
        substring_list.append(substring)
        count_list.append(count)
    
n_list=[]
for k in range(len(substring)):
    n_list.append(k)

zipped_4 = list(zip(header2_list,substring_list,count_list))
df_4 = pd.DataFrame(zipped_4, columns=['Sequence','Substring','Count'])
#print(df_4)

#print(df['Length ORF'].tolist())
max_length_4 = max(df_4['Count'].tolist())
print("Max value %d for sequence:" %max_length_4)

# MAX LENGTH OF ALL FRAMES AND SEQUENCES:
print(df_4.loc[df_4['Count'] == max(df_4['Count'].tolist())])

from collections import Counter
 
def most_frequent(List):
    occurence_count = Counter(List)
    return occurence_count.most_common(1)[0][0]

#The most frequent repeated substring:
print('The most frequent repeated substring:')
print(most_frequent(df_4['Substring'].tolist()))
#The sequences with the most repeated substring:
df_5 = df_4[df_4['Substring'] == most_frequent(df_4['Substring'].tolist())]
print('The sequences with the most repeated substring:') 
print(df_5)
print(df_5['Count'].sum())

