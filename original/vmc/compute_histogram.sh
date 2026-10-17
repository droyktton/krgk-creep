rawfile=$1
archivo=$rawfile.col$3
nbins=$2

#gawk -v n=$3 '(NR>10 && NF==4){print $n}' $rawfile > zzz
gawk -v n=$3 '($3<2){print $n}' $rawfile > zzz
cat zzz | sort -g > $archivo

min=$(head -n 1 $archivo)
max=$(tail -n 1 $archivo)
eve=$(wc -l $archivo | gawk '{print $1}')

echo $archivo, $min, $max, $eve

./compute_histogram $archivo $nbins 1.0 $max $eve

