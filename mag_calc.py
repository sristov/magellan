# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hršak and Strahil Ristov

import sys
import os.path

if len(sys.argv) > 1:
    if os.path.isfile(sys.argv[1]):
        input1 = open(sys.argv[1], 'r')
    else:
        print(sys.argv[1], " does not exist\n")
        exit(1)
else:
    if os.path.isfile("pdg_in.csv"):
        input1 = open('pdg_in.csv', 'r')
    else:
        print("pdg_in.csv does not exist\n")
        exit(1)

FirstRefYear = 1970  # default years of birth for reference population; used if reference_years.txt file is not present
LastRefYear = 2015

Female_gender = "2"
Male_gender = "1"
missing_entry = "0"

IDlist = []
HaplotypedList = []
HaplotypeNamesList = []
FatherMap = {}
MotherMap = {}
YobMap = {}
GenderMap = {}
HaplotypeMap = {}

HAP_column = -1

linebuffer = input1.readline()
linebuffer.rstrip()  # NOT ALWAYS WORKING!?
if linebuffer[len(linebuffer) - 1] == '\n':
    linebuffer = linebuffer[:-1]
lineparts = linebuffer.split(",")

for i in range(len(lineparts)):
    if lineparts[i] == 'ID':
        ID_column = i
    if lineparts[i] == 'father':
        Father_column = i
    if lineparts[i] == 'mother':
        Mother_column = i
    if lineparts[i] == 'YOB':
        YOB_column = i
    if lineparts[i] == 'gender':
        GENDER_column = i
    if lineparts[i] == 'haplotype':
        HAP_column = i

while True:
    linebuffer = input1.readline()
    if not linebuffer:
        break
    linebuffer = linebuffer.rstrip()
    lineparts = linebuffer.split(",")
    if lineparts[ID_column] in IDlist:
        print("duplicate record error ID: ", lineparts[ID_column])
    else:
        IDlist.append(lineparts[ID_column])
    if len(lineparts[Father_column]) == 0:
        lineparts[Father_column] = '0'
    FatherMap[lineparts[ID_column]] = lineparts[Father_column]
    if len(lineparts[Mother_column]) == 0:
        lineparts[Mother_column] = '0'
    MotherMap[lineparts[ID_column]] = lineparts[Mother_column]
    if len(lineparts[YOB_column]) == 0:
        lineparts[YOB_column] = 'MISSING_YEAR'
    YobMap[lineparts[ID_column]] = lineparts[YOB_column]
    GenderMap[lineparts[ID_column]] = lineparts[GENDER_column]
    if HAP_column != -1:
        if lineparts[HAP_column]:
            HaplotypeMap[lineparts[ID_column]] = lineparts[HAP_column]
            HaplotypedList.append(lineparts[ID_column])
            if not lineparts[HAP_column] in HaplotypeNamesList:
                HaplotypeNamesList.append(lineparts[HAP_column])
input1.close()
# checking for errors in pedigree: cycles, gender consistency, and non-existent ancestors
MissingFatherMap = {}
MissingMotherMap = {}
AddedMalesList = []
AddedFemalesList = []

for individual in IDlist:
    if FatherMap[individual] == individual:
        print("Error: cycle in record, please use mag_verif to correct errors.\n")
        exit(31)
    if MotherMap[individual] == individual:
        print("Error: cycle in record, please use mag_verif to correct errors.\n")
        exit(32)
    if FatherMap[individual] != "0":
        if FatherMap[individual] in GenderMap:
            if GenderMap[FatherMap[individual]] != Male_gender:
                print("Error: gender inconsistency, please use mag_verif to correct errors.\n")
                exit(33)
    if MotherMap[individual] != "0":
        if MotherMap[individual] in GenderMap:
            if GenderMap[MotherMap[individual]] != Female_gender:
                print("Error: gender inconsistency, please use mag_verif to correct errors.\n")
                exit(34)
    if FatherMap[individual] != "0":
        if (YobMap[individual] != 'MISSING_YEAR') and (YobMap[FatherMap[individual]] != 'MISSING_YEAR'):
            if int(YobMap[individual]) < int(YobMap[FatherMap[individual]]):
                print("Warning: year of birth inconsistency, please use mag_verif to correct errors.\n")
    if MotherMap[individual] != "0":
        if (YobMap[individual] != 'MISSING_YEAR') and (YobMap[MotherMap[individual]] != 'MISSING_YEAR'):
            if int(YobMap[individual]) < int(YobMap[MotherMap[individual]]):
                print("Warning: year of birth inconsistency, please use mag_verif to correct errors.\n")
    if FatherMap[individual] != "0":
        if FatherMap[individual] not in IDlist:
            if FatherMap[individual] not in MissingFatherMap:
                MissingFatherMap[FatherMap[individual]] = individual
            else:
                MissingFatherMap.pop(FatherMap[individual])
                IDlist.append(FatherMap[individual])
                FatherMap[FatherMap[individual]] = missing_entry
                MotherMap[FatherMap[individual]] = missing_entry
                YobMap[FatherMap[individual]] = missing_entry
                GenderMap[FatherMap[individual]] = Male_gender
                AddedMalesList.append(FatherMap[individual])
    if MotherMap[individual] != "0":
        if MotherMap[individual] not in IDlist:
            if MotherMap[individual] not in MissingMotherMap:
                MissingMotherMap[MotherMap[individual]] = individual
            else:
                MissingMotherMap.pop(MotherMap[individual])
                IDlist.append(MotherMap[individual])
                FatherMap[MotherMap[individual]] = missing_entry
                MotherMap[MotherMap[individual]] = missing_entry
                YobMap[MotherMap[individual]] = missing_entry
                GenderMap[MotherMap[individual]] = Female_gender
                AddedFemalesList.append(MotherMap[individual])
for key in MissingFatherMap:
    FatherMap[MissingFatherMap[key]] = '0'
for key in MissingMotherMap:
    MotherMap[MissingMotherMap[key]] = '0'

if len(AddedMalesList) or len(AddedFemalesList) or len(MissingFatherMap) or len(MissingMotherMap):
    corr_log = open('autocorrection_log.txt', 'w')
    if len(AddedMalesList):
        corr_log.write('new records formed for male ancestors:\n')
        for male in AddedMalesList:
            corr_log.write('     ' + male + '\n')
    if len(AddedFemalesList):
        corr_log.write('new records formed for female ancestors:\n')
        for female in AddedFemalesList:
            corr_log.write('     ' + female + '\n')
    if len(MissingFatherMap):
        corr_log.write('removed unique and non-defined male ancestors:\n')
        for i in MissingFatherMap:
            corr_log.write('     ' + i + '\n')
    if len(MissingMotherMap):
        corr_log.write('removed unique and non-defined female ancestors:\n')
        for i in MissingMotherMap:
            corr_log.write('     ' + i + '\n')
    corr_log.close()

# test for conflicting haplotypes
for i in range(len(HaplotypedList)):
    MotherLine1 = []
    ID_string = HaplotypedList[i]
    MotherLine1.append(ID_string)
    while MotherMap[ID_string] != "0":
        ID_string = MotherMap[ID_string]
        MotherLine1.append(ID_string)

    for j in range(i + 1, len(HaplotypedList)):
        MotherLine2 = []
        ID_string = HaplotypedList[j]
        MotherLine2.append(ID_string)
        while MotherMap[ID_string] != "0":
            ID_string = MotherMap[ID_string]
            MotherLine2.append(ID_string)

        for k in range(len(MotherLine1)):
            CommonAncestorFlag = False
            for l in range(len(MotherLine2)):
                if MotherLine2[l] == MotherLine1[k]:
                    CommonAncestorFlag = True
                    if HaplotypeMap[MotherLine1[0]] != HaplotypeMap[MotherLine2[0]]:
                        print('Error: conflicting haplotypes, please use mag_verif to remove conflicts\n')
                        exit(35)
            if CommonAncestorFlag:
                break

if os.path.isfile('reference_years.txt'):
    input2 = open('reference_years.txt', 'r')
    FirstRefYear = int(input2.readline())
    LastRefYear = int(input2.readline())
    input2.close()

ReferencePopulationList = []
FemalesInReferencePopulationList = []
for individual in IDlist:
    if YobMap[individual] == 'MISSING_YEAR':
        continue
    if FirstRefYear <= int(YobMap[individual]) <= LastRefYear:
        ReferencePopulationList.append(individual)
        if GenderMap[individual] == Female_gender:
            FemalesInReferencePopulationList.append(individual)
if len(ReferencePopulationList) == 0:
    print('ERROR: There are no individuals in the reference population!')
    exit(1)

FounderDamsList = []
for individual in IDlist:
    if MotherMap[individual] == '0':
        if GenderMap[individual] == Female_gender:
            FounderDamsList.append(individual)

FounderDamLineInRefPopList = []
FounderDamLineWithSampledFemaleInRefPopList = []
FounderDamLineWithOnlyMalesInRefPopList = []
IndividualToFounderDamMap = {}
DamLineAllInRefPopCountMap = {}
DamLineFemaleInRefPopCountMap = {}
DamLineHaplotypeMap = {}

for individual in IDlist:
    current = individual
    while MotherMap[current] != '0':
        current = MotherMap[current]
    if GenderMap[current] == Female_gender:  # accounting for the case of founder sire in ref pop
        IndividualToFounderDamMap[individual] = current
    else:
        IndividualToFounderDamMap[individual] = '0'

for i in FounderDamsList:
    DamLineAllInRefPopCountMap[i] = 0
    DamLineFemaleInRefPopCountMap[i] = 0
for individual in ReferencePopulationList:
    current = IndividualToFounderDamMap[individual]
    if current == '0':
        continue
    if GenderMap[current] == Female_gender:  # accounting for the case of founder sire in ref pop
        DamLineAllInRefPopCountMap[current] += 1
        if GenderMap[individual] == Female_gender:
            DamLineFemaleInRefPopCountMap[current] += 1
            if current not in FounderDamLineInRefPopList:
                FounderDamLineInRefPopList.append(current)
            if current in FounderDamLineWithOnlyMalesInRefPopList:
                FounderDamLineWithOnlyMalesInRefPopList.remove(current)
        else:
            if current not in FounderDamLineWithOnlyMalesInRefPopList:
                if current not in FounderDamLineInRefPopList:
                    FounderDamLineWithOnlyMalesInRefPopList.append(current)
        if individual in HaplotypeMap:
            if current not in DamLineHaplotypeMap:
                DamLineHaplotypeMap[current] = HaplotypeMap[individual]
            if GenderMap[individual] == Female_gender:
                if current not in FounderDamLineWithSampledFemaleInRefPopList:
                    FounderDamLineWithSampledFemaleInRefPopList.append(current)
    else:
        print('panic1\n')
        exit(1)

DamLinesWithSamplesInRefPopHaplotypeMap = {}
DamLinesWithSamplesInRefPopCountMap = {}
DamLinesWithFemaleSamplesInRefPopCountMap = {}
for i in FounderDamLineInRefPopList:
    DamLinesWithSamplesInRefPopCountMap[i] = 0
    DamLinesWithFemaleSamplesInRefPopCountMap[i] = 0
for i in FounderDamLineWithOnlyMalesInRefPopList:
    DamLinesWithSamplesInRefPopCountMap[i] = 0
for individual in HaplotypedList:
    if individual in ReferencePopulationList:
        if not IndividualToFounderDamMap[individual] in DamLinesWithSamplesInRefPopHaplotypeMap:
            DamLinesWithSamplesInRefPopHaplotypeMap[IndividualToFounderDamMap[individual]] = HaplotypeMap[individual]
            DamLinesWithSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
            if GenderMap[individual] == Female_gender:
                DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
        else:
            DamLinesWithSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
            if GenderMap[individual] == Female_gender:
                DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1

OneSampleCount = 0
for fdam in DamLinesWithSamplesInRefPopCountMap:
    if DamLinesWithSamplesInRefPopCountMap[fdam] == 1:
        OneSampleCount += 1

SamplesInRefPopCount = 0
SamplesInRefPopFemaleCount = 0
for individual in HaplotypedList:
    if individual in ReferencePopulationList:
        SamplesInRefPopCount += 1
        if GenderMap[individual] == Female_gender:
            SamplesInRefPopFemaleCount += 1

# calculations of the probability of identity and the effective population size

# dam line

NumberOfFemalesInRefPop = len(FemalesInReferencePopulationList)
PI_founder_dams = 1/int(len(FounderDamsList))
PI_dam_line_in_refpop = 0
for fdam in FounderDamLineInRefPopList:
    PI_dam_line_in_refpop += (int(DamLineFemaleInRefPopCountMap[fdam])/NumberOfFemalesInRefPop)**2
if PI_founder_dams != 1:
    increase_in_identity_founder_dams = (PI_dam_line_in_refpop-PI_founder_dams)/(1-PI_founder_dams)
else:
    # only one dam line -> division by zero
    increase_in_identity_founder_dams = 9999

# haplotype line

DamLinesWithoutSamplesList = []
for fdam in FounderDamsList:
    if fdam not in DamLinesWithSamplesInRefPopHaplotypeMap:
        DamLinesWithoutSamplesList.append(fdam)

SumFounderHapLine = len(DamLinesWithSamplesInRefPopHaplotypeMap) + len(DamLinesWithoutSamplesList)
PI_haplotype_in_founders = 0
for hapname in HaplotypeNamesList:
    A = 0
    for fdam in DamLineHaplotypeMap:
        if DamLineHaplotypeMap[fdam] == hapname:
            A += 1
    PI_haplotype_in_founders += (A/SumFounderHapLine)**2
for i in range(len(DamLinesWithoutSamplesList)):
    PI_haplotype_in_founders += (1/SumFounderHapLine)**2

HaplotypeToFounderDamLinesMap = {}
for hapname in HaplotypeNamesList:
    AuxDamLineList = []
    for fdam in DamLinesWithSamplesInRefPopHaplotypeMap:
        if DamLinesWithSamplesInRefPopHaplotypeMap[fdam] == hapname:
            AuxDamLineList.append(fdam)
    HaplotypeToFounderDamLinesMap[hapname] = AuxDamLineList

PI_haplotype_in_refpop = 0
increase_in_identity_haplotype = 9999
if NumberOfFemalesInRefPop > 0:
    for hapname in HaplotypeNamesList:
        A = 0
        for fdam in HaplotypeToFounderDamLinesMap[hapname]:
            A += DamLineFemaleInRefPopCountMap[fdam]
        PI_haplotype_in_refpop += (A/NumberOfFemalesInRefPop)**2
    for fdam in DamLinesWithoutSamplesList:
        PI_haplotype_in_refpop += (DamLineFemaleInRefPopCountMap[fdam]/NumberOfFemalesInRefPop)**2
    if PI_haplotype_in_founders != 1:
        increase_in_identity_haplotype = (PI_haplotype_in_refpop-PI_haplotype_in_founders)/(1-PI_haplotype_in_founders)
    else:
        # only one haplotype line -> division by zero
        increase_in_identity_haplotype = 9999

# sire line

FounderSiresList = []
FounderSiresInRefPopList = []
IndividualToFounderSireMap = {}
SireLineInRefPopCountMap = {}

for individual in IDlist:
    if FatherMap[individual] == '0':
        if GenderMap[individual] == Male_gender:
            FounderSiresList.append(individual)

for individual in IDlist:
    current = individual
    while FatherMap[current] != '0':
        current = FatherMap[current]
    if GenderMap[current] == Male_gender:  # accounting for the case of founder dam in ref pop
        IndividualToFounderSireMap[individual] = current
    else:
        IndividualToFounderSireMap[individual] = '0'

for i in FounderSiresList:
    SireLineInRefPopCountMap[i] = 0
SireLineInRefPopCountMap['0'] = 0

for individual in ReferencePopulationList:
    if GenderMap[individual] != Male_gender:
        continue
    current = IndividualToFounderSireMap[individual]
    if current == '0':
        continue
    if GenderMap[current] == Male_gender:
        SireLineInRefPopCountMap[current] += 1
        if current not in FounderSiresInRefPopList:
            FounderSiresInRefPopList.append(current)
    else:
        print('panic2\n')
        exit(1)

PI_founder_sires = 1/len(FounderSiresList)
PI_sire_line_in_refpop = 0
NumberOfMalesInRefPop = len(ReferencePopulationList) - NumberOfFemalesInRefPop
for fsire in FounderSiresInRefPopList:
    PI_sire_line_in_refpop += (int(SireLineInRefPopCountMap[fsire])/NumberOfMalesInRefPop)**2
if PI_founder_sires != 1:
    increase_in_identity_founder_sires = (PI_sire_line_in_refpop-PI_founder_sires)/(1-PI_founder_sires)
else:
    # only one sire line -> division by zero
    increase_in_identity_founder_sires = 9999

output1 = open('OutputCalc_InputAndResults.txt', 'w')
output1.write('Records in the studbook = ' + '%s'%len(IDlist) + '\n')
output1.write('No. of individuals in the reference population (' + '%s'%FirstRefYear + ' - ' + '%s'%LastRefYear + ') = ' + '%s'%len(ReferencePopulationList) + '\n')
output1.write('No. of female individuals in the reference population = ' + '%s'%len(FemalesInReferencePopulationList) + '\n')
output1.write('No. of founder dams = ' + '%s'%len(FounderDamsList) + '\n')
output1.write('No. of founder dam lines in reference population = ' + '%s'%len(FounderDamLineInRefPopList) + '\n')
output1.write('No. of founder dam lines including lines with only males in reference population = ' + '%s'%(len(FounderDamLineInRefPopList)+len(FounderDamLineWithOnlyMalesInRefPopList)) + '\n')
output1.write('No. of founder dam lines in reference population with samples = ' + '%s'%len(DamLinesWithSamplesInRefPopHaplotypeMap) + '\n')
output1.write('No. of founder dam lines in reference population with only one sample = ' + '%s'%OneSampleCount + '\n')
output1.write('Total number of samples (haplotyped individuals) = ' + '%s'%len(HaplotypedList) + '\n')
output1.write('No. of samples in reference population = ' + '%s'%SamplesInRefPopCount + ':   dams = ')
output1.write('%s'%SamplesInRefPopFemaleCount + ',  sires = ' + '%s'%(SamplesInRefPopCount-SamplesInRefPopFemaleCount) + '\n')
output1.write('----------------------------------------------------------------\n\nDam lines:\n\n')
output1.write('Probability of identity in founder dams = ' + '%s'%PI_founder_dams + '\n')
output1.write('Probability of identity of a dam line in reference population = ' + '%s'%PI_dam_line_in_refpop + '\n')
if increase_in_identity_founder_dams == 9999:
    output1.write('Increase in identity = N/A\n')
    output1.write('Effective dam line size = N/A \n\n')
else:
    output1.write('Increase in identity = ' + '%s'%increase_in_identity_founder_dams + '\n')
    output1.write('Effective dam line size = ' '%s'%(1/increase_in_identity_founder_dams) + '\n\n')
output1.write('Haplotype lines:\n\n')
output1.write('Probability of identity of a haplotype line in founder population = ' + '%s'%PI_haplotype_in_founders + '\n')
output1.write('Probability of identity of a haplotype line in reference population = ' + '%s'%PI_haplotype_in_refpop + '\n')
if increase_in_identity_haplotype == 9999:
    output1.write('Increase in identity = N/A\n')
    output1.write('Effective haplotype line size = N/A \n\n')
else:
    output1.write('Increase in identity = ' + '%s'%increase_in_identity_haplotype + '\n')
    output1.write('Effective haplotype line size = ' '%s'%(1/increase_in_identity_haplotype) + '\n\n')
output1.write('*Sire lines:\n\n')
output1.write('No. of founder sires = ' + '%s'%len(FounderSiresList) + '\n')
output1.write('No. of founder sire lines in reference population = ' + '%s'%len(FounderSiresInRefPopList) + '\n\n')
output1.write('Probability of identity in founder sires = ' + '%s'%PI_founder_sires + '\n')
output1.write('Probability of identity of a sire line in reference population = ' + '%s'%PI_sire_line_in_refpop + '\n')
if increase_in_identity_founder_sires == 9999:
    output1.write('Increase in identity = N/A\n')
    output1.write('Effective sire line size = N/A \n\n')
else:
    output1.write('Increase in identity = ' + '%s'%increase_in_identity_founder_sires + '\n')
    output1.write('Effective sire line size = ' + '%s'%(1/increase_in_identity_founder_sires) + '\n\n')
output1.close()
