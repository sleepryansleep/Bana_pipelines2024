input=$1
percentile=$2

cat $input | sort -n | awk -v p="$percentile" '{all[NR] = $0} END{print all[int(NR*p)]}'
