for bins in 10 100 200 500 1000 2000 5000
do
./histogramador_local.sh lopt_sopt.dat 1 $bins
./histogramador_local.sh lopt_sopt.dat 2 $bins
cat histo.lopt_sopt.dat.sorted.sopt > histo_$bins.sopt
cat histo.lopt_sopt.dat.sorted.lopt > histo_$bins.lopt
done
