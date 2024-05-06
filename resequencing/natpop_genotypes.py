#Use wormtable to call high quality SNPs from the vcf file and make a hapmap input file
from __future__ import division
import wormtable as wt
import sys
import numpy
import pdb

#Define a function to check if record is an indel or not
def indel(record):
    for i in record.ALT:
        if i:
            if len(i) > 1:
                return True
        else:
            return False
    if len(record.REF) > 1:
        return True
    else:
        return False

#Define a function to check if 
def is_het(genotype, separator):
    return (genotype.split(separator)[0] != genotype.split(separator)[1])

#Define a function to check if samples share the same genotype
def is_shared(genotypes):
    return(all(genotypes[i] == '1/1' for i in range(len(genotypes))) or all(genotypes[i] == '0/0' for i in range(len(genotypes))))

#Define a function to check if genotypes are equal
def GT_diff(genotype1, separator1, genotype2, separator2):
    return (genotype1.split(separator1) != genotype2.split(separator2))

#Define a function to check if genotype is alt
#def GT_alt(genotype, separator):

#This function checks that an array of bases is polymorphic
def is_polymorphic(bases):
    if 'NA' in bases:
        un = numpy.unique(bases)
        myind = numpy.where(un == 'NA') #are there NA's in the array?
        result = numpy.delete(un, myind)
        if len(result) >= 2: #two or more alleles
            return (True)
        else:
            return (False)
    else:
        un = numpy.unique(bases)
        if len(un) >= 2:
            return (True)
        else:
            return (False)

def all_NAs(bases):
    numna = sum(bases == 'NA')
    total = len(bases)
    if numna/total >= 0.9: #If 90% or more samples are NA's
        return (True)
    else:
        return (False)


aineisto = wt.open_table('./all.natpop.wt') # open the wormtable
#aineisto = wt.open_table('./small.wt')
f = open('natpop_hapmap.txt', 'w') #Open file mutations.txt for writing
f2 = open('natpop_info.txt', 'w')

#samples = ['10948', '10886', '10932', '1165', '4498', '8816', '3223', '8845', '10908', '847', '10904', '851', '1131', '8850', '8819', '4708', '4712', '6203', '4824', '8783', '8790', '3975', '10928', '10912', '4494', '3210', '10923', '10950', '10951', '10946', '3211', '10906', 'P4452', 'P4463', 'P4468', 'P4471', 'P4476', 'P4479', '10882', '10883', '10884', '10892', '10907', '10914', '10915', '10918', '10925', '10926', '10927', '10935', '10937', '10943', '10983', '3943', 'P4489', '5910', '4730', '1133', '4716'] #samples for natural populations

samples = ['10881', '10882', '10883', '10884', '10885', '10886', '10887', '10888', '10889', '10890', '10891', '10892', '10893', '10894', '10895', '10896', '10897', '10898', '10899', '10900', '10901', '10902', '10903', '10904', '10905', '10906', '10907', '10908', '10909', '10910', '10911', '10912', '10914', '10915', '10916', '10917', '10918', '10919', '10920', '10921', '10922', '10923', '10925', '10926', '10927', '10928', '10929', '10930', '10931', '10932', '10934', '10935', '10936', '10937', '10938', '10939', '10941', '10942', '10943', '10946', '10948', '10950', '10951', '10954', '10982', '10983', '1131', '1132', '1133', '1165', '2229', '3200', '3210', '3211', '3212', '3223', '3943', '3968', '3975', '4494', '4498', '4708', '4712', '4713', '4715', '4716', '4730', '4824', '5910', '5914', '6203', '7833', '847', '851', '8783', '8784', '8787', '8789', '8790', '8816', '8819', '8829', '8845', '8848', '8850', '8851', 'P4451', 'P4452', 'P4457', 'P4459', 'P4463', 'P4468', 'P4471', 'P4472', 'P4476', 'P4479', 'P4486', 'P4489', 'P4496']

##Need to write a header for the hapmap file
f.write('rs' + '\t'+ 'alleles' + '\t' + 'chrom' + '\t' + 'pos' + '\t' + "\t".join(samples) + '\n')

#samples.sort() #The ancestor should always be the first sample [0]
if len(samples) >2:
    print samples
else:
    print "Problem with sample list"
    sys.exit()

#Columns to be retrieved from wormtable, note that genotype are 'samplename.GT' so all these are generated
values = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO.MQ'] + [i+'.GT' for i in samples] + [i+'.DP' for i in samples] + [i+'.GQ' for i in samples] + [i+'.RGQ' for i in samples]
empty = ['', "./.", None, ".", 'None', '*'] #missing_value_character
chrs = ['1', '2', '3', '4', '5', '6', '7'] #Chromosomes to genotype

called_monomorphic = 0
called_polymorphic = 0
snp_counter = 0

#For each row of the wormtable
for row in aineisto.cursor(values):
    chromosome, position, ref, alt, qual, filter, MQ = row[:7]
    
    #print chromosome, position #For debugging
    #Reference genotype qualities
    RGQs = row[7+3*len(samples):7+4*len(samples)]
    DPs = row[7+len(samples):7+2*len(samples)]
    DPs = [0L if i is None else i for i in DPs] #Converts None to 0
    
    #Initialize the array where to write bases
    bases = numpy.empty(len(samples), dtype = object)
    #pdb.set_trace()
    if '.' in chromosome:
        chromosome = chromosome.split(".")[1] #Take chromosome_number
    
    #1. Check that position is polymorphic in the sample
    if alt not in empty and ref not in empty and chromosome in chrs:
        alt = alt.split(",")
        #if type(alt) != list:
        #    print "Error!"
        #    alt = [alt]
        #2. Check that polymorphic site is a SNP, and not something else    
        if len(alt) == 1 and len(alt[0]) == 1 and len(ref) == 1:
            alt=alt[0]
            #position, qual, MQ = int(position), float(qual), float(MQ) #This line ???
            GTs = row[7:7+len(samples)]
            GQs = row[7+2*len(samples):7+3*len(samples)]
            GQs = [0L if i is None else i for i in GQs]
            #3. Check that overall quality of the site is OK
            if numpy.mean(DPs) >= 5 and numpy.mean(GQs) >= 30 and MQ >= 40:
                #4. Check each genotype and write bases
                for i in range(len(samples)):
                    if DPs[i] >= 5 and GQs[i] >= 30 and (is_het(GTs[i], GTs[i][1]) == False) and GTs[i] not in empty:
                        if GTs[i][0] == '1':
                            bases[i] = alt
                        if GTs[i][0] == '0':
                            bases[i] = ref
                    else:
                        bases[i] = 'NA' #If genotype fails quality of heterozygosity check, it is NA
                #5. Check again that there are not too many NA's and site is still polymorphic, if not count it is a monomorphic site
                if is_polymorphic(bases) == True and all_NAs(bases) == False:
                    called_polymorphic += 1 #Increment counter
                    snp_counter += 1
                    #Write genotypes to hapmap file
                    alleles = ''.join([ref, alt]) #Write the two alleles as ref alt
                    #rs = ''.join(['snp', snp_counter])
                    #write bases to hapmap file
                    f.write('snp' + str(snp_counter) + '\t' + alleles + '\t' + str(chromosome) + '\t' + str(position) + '\t' + "\t".join(bases) + '\n')
                if is_polymorphic(bases) == False:
                    called_monomorphic += 1 #Increment number of called monomorphic bases
    #1. If site is monomorphic
    if chromosome in chrs:
        #For monomorphic sites check RGQ of ancestor
        RGQs = [0L if i is None else i for i in RGQs] #Converts None to 0
        #Alternative is to drop 0's [x for x in RGQs if z != None]
        if numpy.mean(RGQs) >= 30 and numpy.mean(DPs) >= 5:
            called_monomorphic += 1 #Increment number of called monomorphic bases
#Store information about called sites
f2.write("Number of called monomorphic sites: " + str(called_monomorphic) + '\n' + "Number of called polymorphic sites: " + str(called_polymorphic) + '\n' + "Total number of called sites: " + str(called_monomorphic + called_polymorphic))

f.close()
f2.close()
### Done        
