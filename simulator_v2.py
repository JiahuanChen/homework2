#! /usr/bin/python
"""
Created on Mon Oct  3 21:30:10 2016

@author: Jiahuan Chen
"""
import sys
import math
import string
import random

length = 70

def checkinput(arg):
    if len(arg) == 1:
        return 1;
    elif arg[1] == "--calc":
        if len(arg) == 3 or arg[3] == "--save":
            return 2
        else:
            return -1
    elif arg[1] == "--load":
        if len(arg) == 3 or arg[3] == "--save":
            return 3
        else:
            return -1

#input: list of ditribution(dictionary) and sequence_length(int)
#function: generate random sequences that can be devided by 3
#return: list of unchecked sequences
def generate_sequences(distribution,sequence_length):
    sequences = []
    count = 0
    for single_dis in distribution:
        seq = []
        for i in range(0,int(round(single_dis["A"]*sequence_length[count]))):
            seq.append("A")
        for i in range(0,int(round(single_dis["T"]*sequence_length[count]))):
            seq.append("T")
        for i in range(0,int(round(single_dis["C"]*sequence_length[count]))):
            seq.append("C")
        for i in range(0,int(round(single_dis["G"]*sequence_length[count]))):
            seq.append("G")
        count += 1;
        random.shuffle(seq)
        s = ''.join(seq)
        while(len(s) %3 != 0):
            s += str(random.choice(["A","T","C","G"]))
        sequences.append(s)
    return sequences

#input: list of unchecked sequences
#function: check start and stop codon
#return: list of checked sequences
def check(sequences):
    new_sequences = []
    for seq in sequences:
        if(seq[0:3] != "ATG"):
            seq = "ATG" + seq[3:]
        for i in range(3,len(seq)-3,3):
            if (seq[i:i+3] == "TAA" 
            or seq[i:i+3] == "TAG" 
            or seq[i:i+3] == "TGA"):
                seq = seq[:i] + random.choice(["A","C","G"]) + seq[i+1:]
        if (seq[len(seq)-3:len(seq)] != "TAA"
        or seq[len(seq)-3:len(seq)] != "TAG"
        or seq[len(seq)-3:len(seq)] != "TGA"):
            seq = seq[:len(seq)-3] + random.choice(["TAA","TAG","TGA"])
        new_sequences.append(seq)
    return new_sequences

#input: list of sequences
#function: computedistribution
#return: list of distribution
def compute_distribution(sequences):
    distribution = []
    for seq in sequences:
        dictionary = {"A":float(0),"C":float(0), "T":float(0),"G":float(0)}
        for nucleotide in seq:
            dictionary[nucleotide] += 1
        for i in dictionary:
            dictionary[i] /= len(seq)
        distribution.append(dictionary)
    return distribution

#input: list of distributions, sequences, names, filename
#function: generate parameter file
#file format: line1: name1; line2: length; line3: distribution(0~1), etc..
def save_parameter(distributions,sequences,names,filename):
    f = open(filename,"w")
    for i in range(0,len(names)):
        f.write(names[i]+"\n")
        f.write(str(len(sequences[i]))+"\n")
        for dis in distributions[i]:
            f.write(dis+ " ")
            f.write(str(distributions[i][dis])+ " ")
        f.write("\n")
    f.close()
    
def print_sequences(sequences,names):
    for i in range(0,len(names)):
        print("%s %d bp" %(names[i], len(sequences[i])))
        for j in range(0, len(sequences[i]), length):
            print(sequences[i][j:j+length])
            
#input: filename(string)
#functions: load names and sequences
#return: list of names and sequences
def openfile(filename):
    file = open(filename)
    sequence = []
    name = []
    count = -1
    for line in file.readlines():
        if line[0] == '>':
            line = line.replace('\n','')
            name.append(line)
            sequence.append('')
            count += 1
        else:
            sequence[count] += line
    for i in range(0,count+1):
        sequence[i] = sequence[i].replace('\n','')
    return name,sequence

#input list of sequences lengths
#function: generate gaussian distribution and then random lengths
#return:list of random sequences lengths
def gaussian(seq_len):
    average_len = 0.0  
    for seq in seq_len:
        average_len += seq
    average_len /= len(seq_len)
    diviation = 0.0    
    for seq in seq_len:
        diviation += (average_len - seq)**2
    diviation = math.sqrt(diviation)
    new_len = [];
    for i in range(0,len(seq_len)):
        new_len.append(random.gauss(average_len,diviation))
    return new_len

#input: filename
#function: get lists of names,sequence_length,distributions for parameter file
#return:lists of names,sequence_length,distributions
def loadparams(filename):
    file = open(filename)
    names = []
    sequence_length = []
    distributions = []
    count = 0
    for line in file.readlines():
        line = line.replace('\n','')
        if count == 0:
            names.append(line)
            count += 1
        elif count == 1:
            sequence_length.append(int(line))
            count += 1
        elif count == 2:
            dictionary = {}
            s = line.split(" ")
            for i in range(0,8,2):
                dictionary[s[i]] = float(s[i+1])
            distributions.append(dictionary)
            count = 0
    return names,sequence_length,distributions


# entrance
if __name__ == '__main__':
    task = checkinput(sys.argv)
    if task == 1:
        distribution = [
        {"A":float(0.25),"C":float(0.25), "T":float(0.25),"G":float(0.25)}]
        sequence_length = [996]
        sequences = generate_sequences(distribution,sequence_length)
        sequences = check(sequences)
        distributions = compute_distribution(sequences)
        names = [">sequence1"]
        filename = "mypara.txt"
        save_parameter(distributions,sequences,names,filename)
        print_sequences(sequences,names)
    elif task == 2:
        names, sequences = openfile(sys.argv[2])
        distributions = compute_distribution(sequences)
        sequence_length = []
        for seq in sequences:
            sequence_length.append(len(seq))
        print(sequence_length)
        sequence_length = gaussian(sequence_length)
        sequences = generate_sequences(distributions,sequence_length)
        sequences = check(sequences)        
        print_sequences(sequences,names)
        if len(sys.argv) == 5:
            save_parameter(distributions,sequences,names,sys.argv[4])
    elif task == 3:
        names,sequence_length,distributions = loadparams(sys.argv[2])
        #use gaussian distribution to generate new lengths  
        #but the diviation will become rediculously large(several hundreds) 
        #thus the random sequences saperate in a very large range
        sequence_length = gaussian(sequence_length)
        sequences = generate_sequences(distributions,sequence_length)
        sequences = check(sequences)        
        print_sequences(sequences,names)
        if len(sys.argv) == 5:
            save_parameter(distributions,sequences,names,sys.argv[4])
    else:
        print("Please check your input!")