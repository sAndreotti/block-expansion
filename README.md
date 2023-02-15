# block-expansion
Algortimo per l'approssimazione di blocchi di aplotipi mediante la pBWT.

##Esecuzione
###Dipendenze
- python3
- numpy
- guppy3 
- tqdm

###Argomenti
- --panel <file> specifica il pannello di aplotipi
- --blocks <file> specifica il file dei blocchi 

###Altro
il pannello viene preso in formato:
> 1 1 1 0 0 1 1 0 1 0 1 0 1 1 0  
> 0 0 1 0 0 0 1 1 0 1 0 1 0 0 0  
> 0 0 1 0 0 1 1 1 0 1 0 1 0 0 0  
> 0 1 0 0 0 1 0 0 1 1 1 0 1 0 0


mentre i blocchi nel formato:
grandezza (start-fine,linea:#sequenza), [indice, indice, ...]
