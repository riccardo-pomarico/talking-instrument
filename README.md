# Talking instrument

Corso: Sound Analysis, Synthesis and Processing

Autori: Riccardo Pomarico, Alice Marazzi.

Il progetto è stato sviluppato utilizzando Matlab e si focalizza sull'analisi del Linear Predictive Coding (LPC), nonché sull'implementazione di un Talking Instrument mediante l'utilizzo della sintesi incrociata basata su LPC. In fase iniziale, sono stati forniti un segnale vocale di test e una traccia monofonica di pianoforte.
È stata condotta un'analisi LPC a breve termine su entrambi i segnali, vocale e musicale. Successivamente, il segnale di errore residuo di ciascun segmento musicale è stato alimentato nel filtro all-pole calcolato da ciascun segmento vocale.
Per la risoluzione di LPC, è stato adottato un approccio di filtraggio di Wiener. Inizialmente, è stata calcolata la soluzione delle equazioni di Wiener-Hopf in forma chiusa, seguita dall'utilizzo del metodo iterativo dello steepest descent.
I filtri di sbiancamento e modellatura ottenuti nel dominio delle frequenze sono stati applicati, evitando l'utilizzo di convoluzione diretta o filtraggio nel dominio del tempo. Poiché il processo viene eseguito frame-by-frame, il segnale di output è stato calcolato mediante l'algoritmo di overlap-and-add (OLA).
Infine, il risultato della sintesi incrociata è stata salvata in un file audio WAV.
