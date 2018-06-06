#!/bin/bash

function usage {
	echo "Usage: $0 -i [sequence_file] -m [master_spark_ip]"
	echo "Optional arguments:"
	echo "	-b (library bound, default: none)"
	echo "	-h (HDFS path, default: /user/root)"
	echo "	-o (output, collection of files)"
	echo "	-p (number of partitions, default: number of sequences)"
	echo "	-t (T-Coffee output, single file)"
	exit 1
}

while getopts “i:m:otp:h:b:” OPTION
do
	case $OPTION in
		i)
			i=$OPTARG
			;;
		m)
			master_ip=$OPTARG
			;;
		p)
			n_partitions=$OPTARG
			;;
		h)
			hdfs_path=$OPTARG
			;;
		o)
			to_file=1
			;;
		t)
			to_file=2
			;;
		b)
			maxlib=$OPTARG
			;;
	esac
done

if [ -z "${i}" ] || [ -z "${master_ip}" ]; then
    usage
fi

if [ -z "${to_file}" ]; then
    to_file=0
fi

if [ -z "${n_partitions}" ]; then
    n_partitions=0
fi

if [ -z "${maxlib}" ]; then
    maxlib=0
fi

if [ -z "${hdfs_path}" ]; then
    hdfs_path=/user/root
fi

if [ ! -f /home/tcoffee/PPCAS_JOS/src/PPCAS.so ]; then
    echo "Shared Library don't exist, run 'make' first."
    exit 1
fi

#spark-submit --files src/PPCAS.so --executor-cores 1 --master spark://${master_ip}:7077 src/PPCAS.py $i $to_file $n_partitions $master_ip $hdfs_path $maxlib
spark-submit --files src/PPCAS.so --executor-cores 1 --master yarn src/PPCAS.py $i $to_file $n_partitions $master_ip $hdfs_path $maxlib
#--driver-memory 1g --executor-memory 1g --num-executors 80

if [ 2 -eq  ${to_file} ]; then
	bash hdfs_comands.sh
	rm -f hdfs_comands.sh
fi