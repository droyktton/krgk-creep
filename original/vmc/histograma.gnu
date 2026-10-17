# calcula maximo y minimo de los datos
#plot file u 0:1 w l
#min = GPVAL_DATA_Y_MIN
#max = GPVAL_DATA_Y_MAX

set table 's.dat'

#NUMBER OF BINS
#nbin=100000

min=1
#max=0.5

#print min, max

binwidth=(max-min)/nbin; 
bin(x,width)= width*floor(x/width) + binwidth/2.0; 
plot [][] file using (bin(($1),binwidth)):(1.0/eve) smooth freq with boxes

unset table

nom=sprintf("histo.%s",file)
print nom
set table nom
plot 's.dat' u ($1):($2) 
unset table


#set logs 
#plot [min:max][1:] nom u 1:2 w his 
#unset logs

