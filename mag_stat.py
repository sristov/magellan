# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn 1.0.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik and Strahil Ristov

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

# test for conflicting haplotypes
for i in range(len(HaplotypedList)):
    MotherLine1 = []
    ID_string = HaplotypedList[i]
    MotherLine1.append(ID_string)
    while MotherMap[ID_string] != "0":
        ID_string = MotherMap[ID_string]
        MotherLine1.append(ID_string)

    for j in range(i+1, len(HaplotypedList)):
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

ReferencePopulationList = []
FemalesInReferencePopulationList = []
for individual in IDlist:
    if YobMap[individual] == 'MISSING_YEAR':
        continue
    if FirstRefYear <= int(YobMap[individual]) <= LastRefYear:
        ReferencePopulationList.append(individual)
        if GenderMap[individual] == Female_gender:
            FemalesInReferencePopulationList.append(individual)

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
IndividualsMissingFounderDamList = []
IndividualsMissingFounderDamInRefPopCount = 0

for individual in IDlist:
    current = individual
    while MotherMap[current] != '0':
        current = MotherMap[current]
    if GenderMap[current] == Female_gender:  # accounting for the case of founder sire
        IndividualToFounderDamMap[individual] = current
    else:
        IndividualToFounderDamMap[individual] = '0'
        IndividualsMissingFounderDamList.append(individual)
        if individual in ReferencePopulationList:
            IndividualsMissingFounderDamInRefPopCount += 1

for i in FounderDamsList:
    DamLineAllInRefPopCountMap[i] = 0
    DamLineFemaleInRefPopCountMap[i] = 0
for individual in ReferencePopulationList:
    current = IndividualToFounderDamMap[individual]
    if current == '0':
        continue
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
        if IndividualToFounderDamMap[individual] not in DamLinesWithSamplesInRefPopHaplotypeMap:
            DamLinesWithSamplesInRefPopHaplotypeMap[IndividualToFounderDamMap[individual]] = HaplotypeMap[individual]
            DamLinesWithSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
            if GenderMap[individual] == Female_gender:
                DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
        else:
            DamLinesWithSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
            if GenderMap[individual] == Female_gender:
                DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1

OneSampleCount = 0
for i in DamLinesWithSamplesInRefPopCountMap:
    if DamLinesWithSamplesInRefPopCountMap[i] == 1:
        OneSampleCount += 1

SamplesInRefPopCount = 0
SamplesInRefPopFemaleCount = 0
for i in range(len(HaplotypedList)):
    if HaplotypedList[i] in ReferencePopulationList:
        SamplesInRefPopCount += 1
        if GenderMap[HaplotypedList[i]] == Female_gender:
            SamplesInRefPopFemaleCount += 1

DamLineMembershipCountMap = {}
for i in FounderDamsList:
    DamLineMembershipCountMap[i] = 0
for individual in IDlist:
    if IndividualToFounderDamMap[individual] != '0':
        DamLineMembershipCountMap[IndividualToFounderDamMap[individual]] += 1

damlineout1 = open('OutputStat_DamLinesWithFemalesInRefPop.txt', 'w')
for fdam in FounderDamLineInRefPopList:
    damlineout1.write(fdam + ': no. of desc. in ref. pop. = ' + '%s' % DamLineAllInRefPopCountMap[fdam])
    if fdam in DamLinesWithSamplesInRefPopHaplotypeMap:
        damlineout1.write('; no. of samples = ' + '%s' % DamLinesWithSamplesInRefPopCountMap[fdam] + ': F = ' + '%s' % DamLinesWithFemaleSamplesInRefPopCountMap[fdam] + ', M = ' + '%s' % (DamLinesWithSamplesInRefPopCountMap[fdam]-DamLinesWithFemaleSamplesInRefPopCountMap[fdam]) + '; haplotype = ' + DamLinesWithSamplesInRefPopHaplotypeMap[fdam] + '\n')
    else:
        damlineout1.write('\n')

damlineout2 = open('OutputStat_DamLinesWithOnlyMalesInRefPop.txt', 'w')
if len(FounderDamLineWithOnlyMalesInRefPopList) == 0:
    damlineout2.write('There are no founder dams with only male descendants in reference population.\n')
for fdam in FounderDamLineWithOnlyMalesInRefPopList:
    damlineout2.write(fdam + ': no. of male desc. in ref. pop. = ' + '%s' % DamLineAllInRefPopCountMap[fdam])
    if fdam in DamLinesWithSamplesInRefPopHaplotypeMap:
        damlineout1.write('; no. of samples = ' + '%s' % DamLinesWithSamplesInRefPopCountMap[fdam] +'; haplotype = ' + DamLinesWithSamplesInRefPopHaplotypeMap[fdam] + '\n')
    else:
        damlineout1.write('\n')

damlineout3 = open('OutputStat_DamLineMembershipAllInRefPop.txt', 'w')
damlineout3.write('founder dam:individual in dam line and in reference population\n')
damlineout4 = open('OutputStat_DamLineMembershipFemaleOnlyInRefPop.txt', 'w')
damlineout4.write('founder dam:female individual in dam line and in reference population\n')
damlineout5 = open('OutputStat_DamLineMembership_1.txt', 'w')
damlineout6 = open('OutputStat_DamLineMembership_2.txt', 'w')
damlineout6.write('founder dam:individual in dam line\n')
for fdam in sorted(DamLineMembershipCountMap, key=DamLineMembershipCountMap.__getitem__, reverse=True):
    if fdam in DamLineHaplotypeMap:
        damlineout5.write("founder dam: " + fdam + "     number of individuals in dam line = " + "%s" % (DamLineMembershipCountMap[fdam]) + "       haplotype = " + DamLineHaplotypeMap[fdam] + "\n")
    else:
        damlineout5.write("founder dam: " + fdam + "     number of individuals in dam line = " + "%s" % (DamLineMembershipCountMap[fdam]) + "       haplotype = N/A\n")
    for individual in IDlist:
        if IndividualToFounderDamMap[individual] == fdam:
            damlineout5.write(individual + "\n")
            damlineout6.write(fdam + ":" + individual + "\n")
            if individual in ReferencePopulationList:
                damlineout3.write(fdam + ':' + individual + '\n')
                if GenderMap[individual] == Female_gender:
                    damlineout4.write(fdam + ':' + individual + '\n')
