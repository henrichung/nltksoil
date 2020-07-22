#this script consists of the commands necessary to read in a biom file 
#with biom's python module and export as tables to input in R


###ON COMMAND LINE
pip install --user numpy
pip install --user cython
pip install --user biom-format

biom convert -i ../../../../work/idoerg/hchung/table.from_biom_w_taxonomy.biom -o table.from_biom.txt --to-tsv
biom convert -i ../../../../work/idoerg/hchung/table.from_biom_w_taxonomy.biom -o otu_table.classic.txt --to-tsv --table-type=‘OTU table’ --header-key taxonomy

#IN PYTHON
import biom
import os

os.chdir("Documents/nltksoil")
table = biom.load_table("data/emp_deblur_90bp.qc_filtered.biom") 

#get column names
a = table.ids()  
b = a.astype("U") 
b.tofile("data/emp_deblur_90bp.qc_filtered_names.csv", sep="\n", format="%s")  

#save as csv (works)
table.metadata_to_dataframe('observation')
df = table.to_dataframe("a")
df.to_csv('myoutput.gz', compression='gzip')

#MOVING FILES AROUND
mv home/hchung/Projects/soil/table.from_biom.txt work/idoerg/hchung
mv data/97_otu_taxonomy.txt ../../../../work/idoerg/hchung
#READING IN R
temp <- fread("../../../../work/idoerg/hchung/table.from_biom_w_taxonomy.txt")
saveRDS(temp, "table.from_biom_w_taxonomy.RDS")
