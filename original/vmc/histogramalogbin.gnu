# calcula maximo y minimo de los datos
#plot file u 0:1 w l
#min = GPVAL_DATA_Y_MIN
#max = GPVAL_DATA_Y_MAX

#logarithmic binning

set table 's.dat'

#NUMBER OF BINS
#nbin=100000

min=1
#max=0.5

#print min, max

binwidth=(log(max/min))/nbin; 
bin(x,width)= width*floor(x/width) + binwidth/2.0; 
plot [][] file using (bin(log($1/min),binwidth)):(1.0) smooth freq with boxes
unset table

# z=log(x), dH(z) = P(z) dz = P(x) dx -> P(z)/x = P(x) 
# dH(z)/dx = P(x) --> dH/dx = dH/(x*dz) = dH/(exp(z)*dz) 

nom=sprintf("histologbin.%s",file)
print nom
set table nom
plot 's.dat' u (exp($1)*min):($2/(exp($1)*min)) 
unset table


#set logs 
#plot [min:max][1:] nom u 1:2 w his 
#unset logs

