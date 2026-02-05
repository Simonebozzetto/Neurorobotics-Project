ciao ragazzi, vi condivido tutto quello che ho fatto e vi lascio questo file readme per capire cos'ho fatto!
Ecco a voi tutte le funzioni e script che ho fatto:

MAIN FILES:

preparation_all_data
preparation_grand_average
questi due script sono il main per calcolare rispettivamente PSD e ERD/ERS dei segnali di tutti i soggetti! in particolare all'interno vengono chiamate alcune funzioni create da me (spero siano anche comprensibili e intuitive dal nome)

grand_average_analysis
cal_feature_identification_extraction
evaluation_phase
questi tre sono gli altri main per rispettivamente la prima parte del grand average analysis (di cui non sono per niente sicuro si svolga in questo modo ed è incompleto), e la seconda parte di identificazione dei features e valutazione (più completa)

FUNCTIONS:

PSDcalculation (calcolo della PSD)
proc_spectrogram (funzione fornita dal professsore)
proc_pos2win (funzione fornita dal professsore)
preprocessing_1 (questa funzione riguarda il preprocessing necessario per calcolare ERD)
preprocessing (questa funzione riguarda il preprocessing necessario per calcolare PSD)
label_vectors_1
label_vectors (stesso significato delle due funzioni sopra)
Fisher (calcolo del Fisher score)
ERD_calculation (calcolo ERD)
cut_samples_concatenate
concatenation 
conc_to_trial

NOTA PER SVOLGERE QUANTO HO FATTO:
bisogna estrarre tutti i file GDF di tutti i soggetti nella stessa cartella, cartella nella quale ci saranno anche tutte le funzioni e main

Per pesantezza non vi carico anche direttamente i file che generano le funzioni ma verranno create direttamente nelle vostre cartelle se eseguite i codici

Ho anche già creato un file latex con già scritto il layout