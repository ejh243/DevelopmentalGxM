'''
Calculate correlations between probes (methylation) and genes (transcription) to generate a probe X gene matrix
'''
import pandas as pd
import calculations as calc

step = 10000

### these data are sorted so that indidivuals are in the same order
methylation_matrix = "/mnt/data1/Helen/GExpxDNAm/Data/DNAm_HYOX.csv" 
rpkm_matrix = "/mnt/data1/Helen/GExpxDNAm/Data/GeneExp_2.csv"

output = "/mnt/data1/Helen/GExpxDNAm/Output/Cor_GExp_Hydroxy_Matrix"

# Methylation data
meth = pd.read_csv(methylation_matrix, sep=",", header=0, index_col = 0,  na_values = "NA")
meth = meth.dropna()
n_meth = len(meth)

print n_meth, " DNA methylation sites read"

start = range(0,n_meth, step)
stop = range(step, n_meth, step)
stop = stop + [n_meth-1]

print len(start), " sections to be processed"


# Read the expression data - need to ensure timepoint match the methylation data
exp = pd.read_csv(rpkm_matrix, header=0, index_col=0)
exp = exp.dropna()

print len(exp), " genes read"
exp_arr = exp.as_matrix()

for a,b in zip(start,stop):

    meth_sub = meth[a:b]

    print "DNA methylation sites", a, "-", b, " selected"

    # Convert to numpy arrays
    meth_arr = meth_sub.as_matrix()

    # perform the matrix-matrix correlation
    corr = calc.corr2_coeff(exp_arr, meth_arr)

    # Go back to a dataframe!
    df = pd.DataFrame(data = corr, index=exp.index, columns=meth_sub.index)
    df.to_csv(output + "-" + str(a) + "-" + str(b) + ".txt", sep='\t')
