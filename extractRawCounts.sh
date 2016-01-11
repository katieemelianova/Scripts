# get read counts from all Salmon quant files and put them in a table with column names as the directory name they came from and feature names in the first column

# For each mapped library, Salmon produces a YourLibraryName.quant directory containing a file called quant.sf with feature, TPM and read count. Put all these directories into one directory together, along with this script and run.

# N.B. assumption is that all the libraries were mapped to same reference. 

# use the dirname to get the library name and extract 4th col from quant.sf file. write to a library specific named file.
for i in *quant
do
	tissueName=$(echo $i | cut -d'.' -f 1)
	tail -n+11 $i'/quant.sf' | cut -f 4 > $tissueName'.rawCounts'
done

# Get one directory of all directories (doesnt matter which) and extract the Feature names (these will be the same for all quant files)
randomDir=$(ls *quant| head -1 | cut -d':' -f 1)
tail -n+11 $randomDir'/quant.sf' | cut -f 1 > Features

# list all the quant directories to get a row of header names to use in final table
dirNames=$(ls -d *quant | cut -d'.' -f 1 | tr '\n' '\t')
echo ' '$dirNames > Header
cat Header | tr ' ' '\t' > Header2

# Put everthng together
paste Features *rawCounts > allRawCountsFeatures
cat Header2 allRawCountsFeatures > allSalmonRawCounts.table

# Remove intermediate files
rm Features Header Header2 *rawCounts allRawCountsFeatures