## Metabolic pathway visualisation with R shiny
***

### What is this app?

This is a Shiny application interface to the R package ipath (github.com/betsig/ipath).

The app allows a user to upload tables of pathway ids (e.g. KEGG KOs) corresponding to genes or proteins, and color pathways depending on correponding values for each gene/protein (e.g. log fold change).


<br> 

### How do I use it?

Navigate through the app using the tabs at the top.

#### Prepare input data
Only specific pathway ids (see table below) are able to be plotted on this metabolic pathway map. 
If annotations are not readily available for your data, we reccomend mapping gene or protein sequences to KEGG using a server such as KAAS. Alternatively, sequences can be blasted against uniprot/swissprot proteins using the pipeline outlined in *Getting pathway IDs from gene and protein sequences*.


| Data type	          | Prefix            | Example ID     |
| --------------------|-------------------|----------------|
| KEGG Pathways       | 	-	 | map00650 |
| KEGG Modules        | 	-	| M00071 |
| KEGG Reactions      | 	-	 | R01148 |
| KEGG Compounds      | 	-	 | C00003 |
| KEGG KOs | 	- | 	K01000 |
| STRING proteins	| - | 	224324.AQ_626 |
| KEGG proteins | 	-	 | aae:aq_626 |
| COGs/eggNOGG OGs | 	-	 | COG0007 |
| Enzyme EC numbers	 | E or EC	 | E2.4.1.82 |
| Uniprot IDs/ACCs | 	UNIPROT: | 	UNIPROT:Q93015 |
| NCBI GI IDs	 | NCBI-GI: | 	NCBI-GI:326314893 |
| CAS Registry Numbers	 | CAS: | 	CAS:7732-18-5 |
| PubChem CIDs	 | PubChem: | 	PubChem:3303 |
| ChEBI IDs	 | -	 | ChEBI:15377 |
| ChEMBL IDs | 	-	 | CHEMBL1098659 |
| PDB-CCD IDs	 | PDB-CCD: | 	PDB-CCD:HOH |
| 3DMET IDs	 | 3DMET:	 | 3DMET:B01124 |
| Nikkaji IDs	 | NIKKAJI: | 	NIKKAJI:J43.587B |
| KNApSAcK IDs | 	KNApSAcK: | 	KNApSAcK:C00001491 |
| LIPID MAPS IDs | 	LIPIDMAPS: | 	LIPIDMAPS:LMFA01060077 |
| Lipidbank IDs | 	LipidBank: | 	LipidBank:VCA0059 |


#### Upload data 
First, upload a dataset with pathway id using the __Upload CSV file with Pathway IDs__ file input. 
If you have a correponding dataset with a feature you want to colour you pathways by, upload it to the second __[Optional] Upload CSV file of corresponding gene/protein differential expression values__. If sequence ids are present in both datasets, these will be automatically joined together. Use the __Choose column name with sequence ids__ inputs to select which columns contain sequence ids if these are incorrectly inferred. 

#### View pathways
Note that mapping and plotting of pathways can take up to 20 seconds. Please be patient in this tab. 

__Choose variable containing pathway IDs__
First, select column name containing the pathway ids (e.g. "K01000"). 

__Choose variable to color pathways__
If you would like to colour the mapped pathways by a variable in your data, select the column here.

__Choose variable to change width of pathways__
If you would like to change the widths the mapped pathways by a variable in your data, select the column here. By default, all pathways will be the same width.

__Colour type__
Pick a colour palette type (discrete or continous).

__Colour palette__
Select from a range of color palettes.


#### Download PDF or SVG
Click the __Download SVG__ or __Download PDF__ button to download your pathway as and svg or pdf file. 


#### Getting pathway IDs from gene and protein sequences
See: COMING SOON.

<br>

