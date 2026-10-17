set xlabel 'q'; set ylabel 'Sq'
set xra [0.05:3.14]
set title 'Lz=128, Lx=1024, M=Lx/Lz=8'
plot './ground_partial.fft.perio' w p t 'periodic M=8', \
(1./x)**(1+2*0.5) t 'zeta=0.5',\
'./ground_partial.fft.chain' w p t 'chain M=8'
