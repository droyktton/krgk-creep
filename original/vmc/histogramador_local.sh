# uso: Para hacer un histograma de la columna numero 'col' de un 'file' 
# ./histogramador.sh file col nbins

# uso: para hacer histogramas de lopt y sopt para todos los files del tipo L*f*.dat
# for file in L*f*.dat; do ./histogramador.sh $file 1 nbins; ./histogramador.sh $file 2 nbins; done

rawfile=$1
nbins=$3

if(($2==1)) 
then
ext="lopt"
fi
if(($2==2)) 
then
ext="sopt"
fi
if(($2==3)) 
then
ext="lrel"
fi
if(($2==4)) 
then
ext="srel"
fi
archivo=$rawfile.sorted.$ext
echo $(wc -l $rawfile)
echo $rawfile "->" $archivo 

#gawk -v n=$2 '(NR>10 && $3<2){print $n}' $rawfile > zzz
gawk -v n=$2 '(NR>10){print $n}' $rawfile > zzz
cat zzz | sort -g > $archivo
#echo "ordenamiento ok"

min=$(head -n 1 $archivo)
max=$(tail -n 1 $archivo)
eve=$(wc -l $archivo | gawk '{print $1}')
echo $archivo, $min, $max, $eve

./compute_histogram $archivo $nbins $min $max $eve

#gnuplot -persist << EOFGNU
#		reset
#		file='$archivo'
#		min=$min
#		max=$max
#		eve=$eve
#		nbin=$nbins		
#		load 'histograma.gnu'
#
#		reset
#		file='$archivo'
#		min=$min
#		max=$max
#		eve=$eve
#		nbin=$nbins		
#	        load 'histogramalogbin.gnu'
#EOFGNU

