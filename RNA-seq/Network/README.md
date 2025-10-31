# Network

## network preparation
- obtain network from [https://skm.nib.si/](https://skm.nib.si/) downloads section
- use [PSS-tools: PSS-refinements.ipynb](https://github.com/NIB-SI/skm-tools/blob/af75f2cff9e0a23b92323806349707f097641cc6/examples/PSS-refinements.ipynb) to simplify network (e.g. remove complexed, metabolites, etc); optional (*)
- use [DiNARs' pre-processing app](https://github.com/NIB-SI/DiNAR/) shiny subApp

## prioritisation
- use PSS-translation-prioritisation script on output from DE analysis to translate potatoIDs to arabidopsisIDs to functional clusterIDs ([SKM: PSS](https://skm.nib.si/))

## visualisation
- [Cytoscape](https://cytoscape.org/) or
- [DiNAR](https://github.com/NIB-SI/DiNAR/) shiny app

__

(*) PSS-refinement
- remove complexes: 
```
['CML|HC-Pro' ,'EDS1|MPK3|PAD4' ,'CRT|ETR' ,'HSP90|RAR1|SGT1' ,'NDR1|RIN4' ,'RANGAP|Rx' ,'GPAphid2|RANGAP' ,'OBE1|WRKY17' ,'OBE1|WRKY11' ,'CO|OBE1' ,'MKS1|WRKY33' ,'WRKY30|WRKY53' ,'R-gene|potyvirus' ,'EIN3(like)|JAZ' ,'DELLA|GA|GID1' ,'SARD1|TCP8' ,'TCP8|WRKY28' ,'NAC019|TCP8' ,'NPR1|WRKY18' ,'BIK1|PEPR1|PEP1' ,'ET|ETR']
```
- keep types:
```
['Complex', 'Metabolite', 'PlantAbstract', 'PlantCoding', 'PlantNonCoding', 'Process']
```
- remove 'uninteresting' metabolites, i.e keep metabolites :
```
['ROS', 'Ca2+']
```
- collapse rections:
```
["catalysis", "translocation"]
```
