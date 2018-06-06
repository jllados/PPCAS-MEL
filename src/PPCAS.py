from pyspark import SparkConf, SparkContext
import sys
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import time
import math

def parse_input(file_name):
	name = None
	seq = ""
	count = 0
	seq_list = []
	
	with open(file_name) as f:
		for line in f:
			line = line.strip()
			if(len(line)>0):
				if(line[0]==">"):
					seq_list.append((count,name,seq))
					name=line[1:]
					seq= ""
					count = count + 1
				else: seq=seq+line
		seq_list.append((count,name,seq))
		seq_list.pop(0)
	return seq_list
	
def tc_header(out_name, seq_list):
	file_ = open(out_name + ".lib", 'w')
	file_.write('! TC_LIB_FORMAT_01\n')
	file_.write(str(len(seq_list)) + "\n")
	for i in xrange(len(seq_list)):
		file_.write(seq_list[i][1] + " " + str(len(seq_list[i][2])) + " " + seq_list[i][2] + "\n")
	file_.close()
	
def hdfs_cmd(out_name):
	file_2 = open("hdfs_comands.sh", 'w')
	file_2.write("hdfs dfs -cat " + hdfs_path + "/" + out_name + "/part-* >> " + out_name + ".lib\n")
	file_2.write("hdfs dfs -rm -r " + out_name + "\n")
	file_2.close()
	
def generate_tasks(n_seq, n_part):
	n_comb = float((float(n_seq) * (float(n_seq) - 1)) / 2)
	max_it = math.ceil(n_comb/n_part)
	it = 0

	tasks = []
	for i in range(0, n_seq):
		for j in range(i + 1, n_seq):
			if(it == 0):
				tasks.append((i, j, int(max_it)))
			
			it +=1
			if(max_it == it):
				n_comb -= max_it
				if(n_comb <= 0): break
				n_part -= 1
				max_it = math.ceil(n_comb/n_part)
				it = 0
	return tasks

def main(sc, seq, to_file, n_partitions, m_ip, hdfs_path, maxlib):
	
	def ProbaMatrix2CL(row):
		n_t1 = row[0]
		n_t2 = row[1]
		max_it = row[2]
		it = 0
		
		libstr = ""
		mapped_list = []
		try:
			for t1 in SEQ_LIST.value[n_t1:]:
				for t2 in SEQ_LIST.value[n_t2:]:
					lib = np.empty((len(t1[2])*len(t2[2]),3), dtype='int32')
					lib = lib.ravel()

					_libraryC = ctypes.CDLL("./PPCAS.so")
					array_1d = np.ctypeslib.ndpointer(ctypes.c_int,flags="C_CONTIGUOUS")
					_libraryC.proba_pair_wise.argtypes = [ctypes.c_char_p, ctypes.c_char_p, array_1d, ctypes.c_int]
					list_n = _libraryC.proba_pair_wise(t1[2], t2[2], lib, MAXCONS.value)
					lib = np.resize(lib,(list_n, 3))
					
					for line in lib:
						mapped_list.append((str(t1[0]) + "\t" + str(line[0].item()) + "\t" + str(t2[0]) + "\t" + str(line[1].item()) + "\t" + str(line[2].item())))
						
					if(to_file == 1):
						libstr += "Parwise " + str(t1[0]) + " " + str(t2[0]) + "\n"
						for line in lib:
							libstr += "\t" + str(line[0].item()) + "\t" + str(line[1].item()) + "\t" + str(line[2].item()) + "\n"
					elif(to_file == 2):
						libstr += "#" + str(t1[0]) + " " + str(t2[0]) + "\n"
						for line in lib:
							libstr += "\t" + str(line[0].item()) + "\t" + str(line[1].item()) + "\t" + str(line[2].item()) + "\t1" + "\t0\n"
					del lib
					it+=1
					if(it == max_it): raise StopIteration()
				n_t2 = t1[0] + 1
		except StopIteration:
			pass
		return libstr[:-1]

	out_name = seq + "_" + str(time.time())
	seq_list = parse_input(file_name)
	
	
	#SORT
	seq_list = sorted(seq_list, key=lambda constraint: constraint[1])
	for i in range(len(seq_list)):
		seq_list[i] = (i+1,seq_list[i][1],seq_list[i][2])
	
	SEQ_LIST = sc.broadcast(seq_list)
	TO_FILE = sc.broadcast(to_file)
	
	if(n_partitions):
		tasks = generate_tasks(len(seq_list), float(n_partitions))
	else:
		tasks = generate_tasks(len(seq_list), float(len(seq_list)))
	rdd_tasks = sc.parallelize(tasks, len(tasks))
	
	maxPairs = (len(seq_list) * (len(seq_list)-1)) / 2
	maxCons = maxlib/maxPairs
	MAXCONS = sc.broadcast(maxCons)
	
	if(to_file==1):
		rdd_tasks.map(ProbaMatrix2CL).saveAsTextFile("hdfs://" + m_ip + ":8020" + hdfs_path + "/" + out_name) #May use :9000
	elif(to_file==2):
		tc_header(out_name, seq_list)
		hdfs_cmd(out_name)
		rdd_tasks.map(ProbaMatrix2CL).saveAsTextFile("hdfs://" + m_ip + ":8020" + hdfs_path + "/" + out_name) #May use :9000
	else: 
		rdd_tasks.flatMap(ProbaMatrix2CL).collect()
		
	del seq_list[:]
	
if __name__ == "__main__":
	file_name = sys.argv[1]
	seq = file_name.rpartition('/')[2].split('.')[0]
	to_file = int(sys.argv[2])
	n_partitions = int(sys.argv[3])
	m_ip = sys.argv[4]
	hdfs_path = sys.argv[5]
	maxlib = int(sys.argv[6])
	
	APP_NAME = "PPCAS_" + seq 
	conf = SparkConf().setAppName(APP_NAME)
	sc   = SparkContext(conf=conf)
	main(sc, seq, to_file, n_partitions, m_ip, hdfs_path, maxlib)