
workingDir = getwd()
workingDir = "C:/Dropbox/Academic/Manitoba/Biological networks/data/GeneFunction_DomainInfo/"
#Pathway file comes from BioPlanet
dir2<-"C:/Dropbox/Academic/Manitoba/Biological networks/data/PathwayCrosstalkFilesForCuneyt/"
pathwayFile = paste(dir2,"bioplanet_pathway.csv",sep="");
bioPlanetsRdata<-"C:/Dropbox/Academic/Manitoba/Biological networks/data/BioPlanetNetworks.RData"

#Clusters come from Mark's experiments
clusterFile = paste(workingDir,"sites_by_cluster.txt",sep="");

#Files to be used in gene similarity computations
#I needed to do the following modifications on the original file:
#1 - The line endings contained many repeating \t. remove them. 
#2- Change line endings with "\t\t\r\n" to "\r\n"
#3- Change pathway name and remove commas in the name
#4- Change pathway name and gene list seperator from "\t\t" to ","
subcellularSimFile =paste(workingDir,"GO_Cellular_Component_2018.txt",sep="");
molecularSimFile = paste(workingDir,"GO_Molecular_Function_2018.txt",sep="");
proDomainSimFile = paste(workingDir,"InterPro_Domains_2019.txt",sep="");
bioProcessSimFile = paste(workingDir,"GO_Biological_Process_2018.txt",sep="");

#Target Fasn pathways. This list comes from Guolin.
targetFasn<-c("p73 transcription factor network", "Delta Np63 pathway", "Ghrelin pathway",
"SREBF and miR-33 in cholesterol and lipid homeostasis", "Metabolism", "Insulin signaling pathway",
"Fatty acid biosynthesis", "AMPK signaling", "Triglyceride biosynthesis", "Fatty acyl-CoA biosynthesis",
"Lipid and lipoprotein metabolism", "Metabolism of vitamins and cofactors",
"Integration of energy metabolism", "Vitamin B5 (pantothenate) metabolism",
"Fatty acid, triacylglycerol, and ketone body metabolism", "ChREBP activates metabolic gene expression")
    
  
