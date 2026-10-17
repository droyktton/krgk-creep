/*
Calcula el histograma de datos reales o enteros en un file de una columna

Para Compilar, es necesario tener instalado la GNU Scientific Library (GSL) -- development package (libgsl0-dev en ubuntu):

para valores enteros:
g++ compute_histogram.cpp -lgsl -lgslcblas -lm -o compute_histogram

para valores reales:
g++ -DDOUBLE compute_histogram.cpp -lgsl -lgslcblas -lm -o compute_histogram

Uso:

./compute_historgam file nbins min max events

donde 
-- "file" es de una columna y es el que contiene los datos
-- "nbins" es el numero de bins de los histogramas
-- "min" es el valor minimo a histogramear (no puede ser 0)
-- "max" es el valor maximo a histogramear
-- "events" es el numero de eventos a considerar del file

Resultado:
Dos ficheros de dos columnas, "histo_normal.nombrefile", "histo_logbin.nombrefile"
La primer columna es el centro del bin, la segunda es la frecuencia normalizada

*/


#include <gsl/gsl_histogram.h>
#include <iostream>
#include <fstream>
#include<cstdlib>
#include<cmath>

using namespace std;

void histogram_fprintf (FILE * histovelout, gsl_histogram *histovel, int events)
{
	double suma=gsl_histogram_sum(histovel); 
	size_t nrobins=gsl_histogram_bins(histovel);	

	int sum=0.0;
	for(int i=0;i<nrobins;i++) sum+=gsl_histogram_get(histovel,i);

	for(int i=0;i<nrobins;i++)
	{
		 double upper,lower,bin, bin2;
		 gsl_histogram_get_range(histovel, i, &lower, &upper);
		 bin=gsl_histogram_get(histovel,i)/(events*(upper-lower));
		 bin2=gsl_histogram_get(histovel,i)/(sum*(upper-lower));
		 fprintf(histovelout,"%g %g %g\n", 0.5*(lower+upper),bin,bin2); 
	}
	fprintf(histovelout,"\n\n");
}
int main(int argc, char **argv)
{
	if(argc!=6){
		cout << "uso:\n ./a.out file nbins min max events\n\n";
		exit(1);
	}
	
	char *filenom=argv[1];
	int NBINS=atoi(argv[2]);
	double MINIMO=atof(argv[3]);
	double MAXIMO=atof(argv[4]); 
	double EVENTS=atof(argv[5]); 

	cout << "file= " << argv[1] << endl;
	cout << "NBINS= " << NBINS << endl;
	cout << "MAX= " << MAXIMO << endl;
	cout << "MIN= " << MINIMO << endl;
	cout << "Events= " << EVENTS << "\n\n";

	if(MINIMO==0) cout << "warning: MIN no puede ser 0 para logarithmic binning\n";
	cout.flush();

	// histograma con logarithmic binning
	gsl_histogram * histologbin = gsl_histogram_alloc(NBINS);
	double *range=(double *)malloc((NBINS+1)*sizeof(double));
	for(int i=0;i<NBINS+1;i++) range[i]=MINIMO*powf(MAXIMO/MINIMO,i*1.0/NBINS);
	gsl_histogram_set_ranges (histologbin, range, NBINS+1);

	// histograma con uniform binning
	gsl_histogram * histo = gsl_histogram_alloc(NBINS);
	gsl_histogram_set_ranges_uniform (histo, MINIMO, MAXIMO);
	
	ifstream input(filenom);

	#ifndef DOUBLE
	int value;
	#else
	double value;
	#endif

	input >> value;
	cout << "reading... ";cout.flush();
	int counter=0;
	while(!input.eof())
	{
		//cout << value << "\n";
		gsl_histogram_increment (histo,value);
		gsl_histogram_increment (histologbin,value);
		input >> value;	
		counter++;
		if(counter% int(EVENTS/10) == 0) {
			cout <<  10*(int)(counter/int(EVENTS/10)) << "% -> ";
			cout.flush();
		}
	}
	cout << "ok\n";

	char nom[100];
	sprintf(nom,"histo_normal.%s",filenom);
	FILE *histo_out=fopen(nom,"w");
	cout << "outputs en " << nom << " y en ";

	sprintf(nom,"histo_logbin.%s",filenom);
	FILE *histologbin_out=fopen(nom,"w");
	cout << nom << "\n";

	//gsl_histogram_fprintf (histo_out, histo , "%g", "%g");
	//gsl_histogram_fprintf (histologbin_out, histologbin , "%g", "%g");

	histogram_fprintf(histo_out, histo, EVENTS);
	histogram_fprintf(histologbin_out, histologbin, EVENTS);

	fclose(histo_out);fclose(histologbin_out);
}

