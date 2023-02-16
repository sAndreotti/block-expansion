# Block-expansion
Algortimo per l'approssimazione di blocchi di aplotipi mediante la pBWT.

## Esecuzione
###### Dipendenze
- python3
- numpy
- guppy3 
- tqdm

###### Argomenti
- --panel <file> specifica il pannello di aplotipi
- --blocks <file> specifica il file dei blocchi 

###### Esempio esecuzione
> block-expansion.py --panel file_pannello --blocks file_blocchi
  
## Altro
il pannello viene preso in formato:
> 1 1 1 0 0 1 1 0 1 0 1 0 1 1 0  
> 0 0 1 0 0 0 1 1 0 1 0 1 0 0 0  
> 0 0 1 0 0 1 1 1 0 1 0 1 0 0 0  
> 0 1 0 0 0 1 0 0 1 1 1 0 1 0 0


mentre i blocchi nel formato:
> grandezza (start-fine,indice:#sequenza), [indice, indice, ...]
  
Grandezza: grandezza totale del blocco calcolata come 
>numero aplotpi x (fine-inizio)
  
Start: sito di inizio blocco
  
Fine: sito di fine blocco
  
Sequenza: caratteri di un aplotipo del blocco, in caso il blocco abbia un grado maggiore di 0 i siti nei quali sono presenti variazioni sono contrassegnati da una 'X'
>'00100' oppure '0X100'
[Indice, indice, ...]:lista degli indici del blocco
