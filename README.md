# SchemasDeductifsEnAMCD

Le présent dépôt contient l'intégralité des algorithmes implémentés dans le cadre de la rédaction du manuscrit de thèse intitulé : Explications en aide multicritère à la décision : schémas déductifs, algorithmes et expérimentations.

Il comporte 3 répertoires principaux :

A- Le répertoire CHAPITRES2et3

  Ce dossier contient principalement:
  
    - un fichier contenant l'ensemble des 15 algorithmes de calcul d'explication implémentés : *Reference_Functions.py*
    
    - un fichier exécutable (avec une fonction main) qui reproduit les expérimentations réalisées dans le Chapitre 3 : *Fichier_principal.py*
    
    - 3 répertoires qui fournissent, représentées sous la forme d'un séquence d'entiers, toutes les relations d'ordre linéaire additives définies sur 4, 5 et 6 items : *CBTO4,CBT05 et CBTO6*
    
    - un fichier contenant l'implémentation de l'algorithme, décrit en annexe, et permettant la génération de toutes les relations d'ordre linéaire additives définies sur 4, 5 et 6 items : *CBTO_Generator.py*


B- Le répertoire CHAPITRE4

  Ce dossier contient principalement:
  
    - un dossier *XP* qui contient les détails des explications calculées lors des expérimentations numériques du Chapitre 4.
    
    - un fichier contenant l'ensemble des algorithmes de calcul d'explication de déduction de la relation nécessaire : *Reference_Functions_Robuste.py*

    - un fichier exécutable (avec une fonction main) qui reproduit les expérimentations réalisées dans le Chapitre 4 : *Fichier_principal_Execution.py*

    - un fichier *FonctionElicitation.py* qui contient le code Python de Programme Linéaire en Nombres Entiers décrit dans la sous-section 4.4.4.


C- Le répertoire CHAPITRE5

  Ce dossier contient principalement:
  
    - un dossier *DONNEES* qui contient les instances générées dans le cadre des expérimentations numériques du Chapitre 5.
    
    - deux fichiers exécutables *Fichier_principal_Execution.py* et *Execution_preference_swaps.py* qui permettent de reproduire les expérimentations conduites dans ce chapitre.


**NB: Le code nécessite pour son exécution que soit installé le module gurobipy du solveur Gurobi dans sa version 9.1.0 ou une version postérieure**



    
