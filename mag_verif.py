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
ConflictingFatherList = []
ConflictingMotherList = []
FemaleFathersList = []
MaleMothersList = []
YoungerFathersList = []
YoungerMothersList = []

for individual in IDlist:
    if FatherMap[individual] == individual:
        ConflictingFatherList.append(individual)

    if MotherMap[individual] == individual:
        ConflictingMotherList.append(individual)

    if FatherMap[individual] != "0":
        if FatherMap[individual] in GenderMap:
            if GenderMap[FatherMap[individual]] != Male_gender:
                FemaleFathersList.append(individual)

    if MotherMap[individual] != "0":
        if MotherMap[individual] in GenderMap:
            if GenderMap[MotherMap[individual]] != Female_gender:
                MaleMothersList.append(individual)

    if FatherMap[individual] != "0":
        if (YobMap[individual] != 'MISSING_YEAR') and (YobMap[FatherMap[individual]] != 'MISSING_YEAR'):
            if int(YobMap[individual]) < int(YobMap[FatherMap[individual]]):
                YoungerFathersList.append(individual)

    if MotherMap[individual] != "0":
        if (YobMap[individual] != 'MISSING_YEAR') and (YobMap[MotherMap[individual]] != 'MISSING_YEAR'):
            if int(YobMap[individual]) < int(YobMap[MotherMap[individual]]):
                YoungerMothersList.append(individual)

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

FatalError = False
if len(ConflictingFatherList):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('Process is terminated. Please correct the following errors:\n')
    for individual in ConflictingFatherList:
        error_log.write('Individual ' + individual + 'is repeated as its father.\n')
    print("error1, please see ERROR_ALERT.TXT")
    error_log.close()
    FatalError = True

if len(ConflictingMotherList):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('Process is terminated. Please correct the following errors:\n')
    for individual in ConflictingMotherList:
        error_log.write('Individual ' + individual + 'is repeated as its mother.\n')
    error_log.close()
    print("error2, please see ERROR_ALERT.TXT")
    FatalError = True

if len(FemaleFathersList):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('Process is terminated. Please correct the following errors:\n')
    for individual in ConflictingMotherList:
        error_log.write('Process is terminated. Please correct the following error:\n')
        error_log.write(
            "Individual " + individual + 'has a father \n' + FatherMap[individual] + 'that is defined as female.\n')
    error_log.close()
    print("error3, please see ERROR_ALERT.TXT")
    FatalError = True

if len(MaleMothersList):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('Process is terminated. Please correct the following errors:\n')
    for individual in ConflictingMotherList:
        error_log.write('Process is terminated. Please correct the following error:\n')
        error_log.write(
            "Individual " + individual + 'has a mother \n' + MotherMap[individual] + 'that is defined as male.\n')
    error_log.close()
    print("error4, please see ERROR_ALERT.TXT")
    FatalError = True

if len(YoungerFathersList):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('WARNING! Please correct the following errors:\n')
    for individual in YoungerFathersList:
        error_log.write(
            'Individual ' + individual + ' has younger father ' + FatherMap[individual] + '\n')
    error_log.close()
    print("warning1, please see ERROR_ALERT.TXT")

if len(YoungerMothersList):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('WARNING! Please correct the following errors:\n')
    for individual in YoungerMothersList:
        error_log.write(
            'Individual ' + individual + ' has younger mother ' + MotherMap[individual] + '\n')
    error_log.close()
    print("warning2, please see ERROR_ALERT.TXT")

def checkConflicts(ParentMap,GenderMap,LineType):
    OneGenderList = [key for key, value in GenderMap.items() if value == LineType]
    ConflictMap = {}
    for individual in OneGenderList:
        ID_string = individual
        LineageList = [individual]
        while ParentMap[ID_string] != '0':
            if ParentMap[ID_string] in LineageList:
                ConflictMap[individual] = ParentMap[ID_string]
                LineageList.append(ParentMap[ID_string])
                break
            ID_string = ParentMap[ID_string]
            LineageList.append(ID_string)

    return ConflictMap

FatherLineConflicts = checkConflicts(FatherMap,GenderMap,Male_gender)
MotherLineConflicts = checkConflicts(MotherMap,GenderMap,Female_gender)

if len(FatherLineConflicts):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('Please correct the following errors:\n')
    for individual in FatherLineConflicts:
        error_log.write("Male individual {0} has a repeated ancestor {1} in its lineage:\n".format(
            individual, FatherLineConflicts[individual]))
        error_log.write("\n\n")
    error_log.close()
    print("error5, please see ERROR_ALERT.TXT")
    FatalError = True

if len(MotherLineConflicts):
    error_log = open('ERROR_ALERT.TXT', 'w')
    error_log.write('Please correct the following errors:\n')
    for individual in MotherLineConflicts:
        error_log.write("Female individual {0} has a repeated ancestor {1} in its lineage:\n".format(
            individual, MotherLineConflicts[individual]))
        error_log.write("\n\n")
    print("error6, please see ERROR_ALERT.TXT")
    error_log.close()
    FatalError = True

if FatalError:
    exit(30)

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

out_conflicts = open('OutputVerif_ConflictingIndividuals.txt', 'w')
out_conflicts_detail = open('OutputVerif_DetailedConflictingIndividuals.txt', 'w')
out_summary = open('OutputVerif_Summary.txt', 'w')
out_misplaced = open('OutputVerif_MisplacedBranches.txt', 'w')

MatchCount = 0
MismatchCount = 0
InformativeList = []
UnitsInConflictsList = []
ConflictCountMap = {}
ConflictIndividualsMap = {}

for individual in HaplotypedList:
    ConflictCountMap[individual] = 0
    ConflictIndividualsMap[individual] = []
ConflictMatrix = [[0 for i in range(len(HaplotypedList))] for j in range(len(HaplotypedList))]

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

        # find common ancestor
        for k in range(len(MotherLine1)):
            CommonAncestorFlag = False
            for l in range(len(MotherLine2)):
                if MotherLine2[l] == MotherLine1[k]:
                    pathlength = k + l
                    for m in range(k):
                        if MotherLine1[m] not in InformativeList:
                            InformativeList.append(MotherLine1[m])
                    for m in range(l):
                        if MotherLine2[m] not in InformativeList:
                            InformativeList.append(MotherLine2[m])
                    if HaplotypeMap[MotherLine1[0]] == HaplotypeMap[MotherLine2[0]]:
                        MatchCount += 1
                    else:
                        ConflictCountMap[MotherLine1[0]] += 1
                        ConflictCountMap[MotherLine2[0]] += 1
                        if MotherLine1[k] not in ConflictIndividualsMap[MotherLine1[0]]:
                            ConflictIndividualsMap[MotherLine1[0]].append(MotherLine1[k])
                        if MotherLine2[l] not in ConflictIndividualsMap[MotherLine2[0]]:
                            ConflictIndividualsMap[MotherLine2[0]].append(MotherLine2[l])
                        ConflictMatrix[i][j] = 1
                        ConflictMatrix[j][i] = 1
                        MismatchCount += 1
                        if MotherLine1[0] not in UnitsInConflictsList:
                            UnitsInConflictsList.append(MotherLine1[0])
                        if MotherLine2[0] not in UnitsInConflictsList:
                            UnitsInConflictsList.append(MotherLine2[0])
                    CommonAncestorFlag = True
                    break
            if CommonAncestorFlag:
                break
if len(UnitsInConflictsList) > 0:
    out_conflicts.write('ID of the conflicting unit; Haplotype; No. of conflicts\n\n')
    out_conflicts_detail.write('ID of the conflicting unit; Haplotype; Conflicting members\n\n')
else:
    out_conflicts.write('There are no conflicts in the pedigree.\n')

ConflictCountDecMap = {}
for individual in HaplotypedList:
    ConflictCountDecMap[individual] = ConflictCountMap[individual]
ConflictingUnitsList = []
MisplacedUnitsList = []
maxConflictUnitIndex = 0
while True:
    maxConflictCount = 0
    for i in range(len(HaplotypedList)):
        if ConflictCountDecMap[HaplotypedList[i]] > maxConflictCount:
            maxConflictCount = ConflictCountDecMap[HaplotypedList[i]]
            maxConflictUnit = HaplotypedList[i]
            maxConflictUnitIndex = i
    if maxConflictCount == 0:
        break
    ConflictingUnitsList.append(maxConflictUnit)
    ConflictCountDecMap[maxConflictUnit] = 0
    out_conflicts.write(
        maxConflictUnit + '; ' + HaplotypeMap[maxConflictUnit] + '; ' + '%s' % ConflictCountMap[maxConflictUnit] + '\n')
    out_conflicts_detail.write(
        maxConflictUnit + '; ' + HaplotypeMap[maxConflictUnit])
    for unit in ConflictIndividualsMap[maxConflictUnit]:
        out_conflicts_detail.write('; {0}'.format(unit))
    out_conflicts_detail.write('\n')
    for j in range(len(HaplotypedList)):
        if ConflictMatrix[maxConflictUnitIndex][j] == 1:
            ConflictCountDecMap[HaplotypedList[j]] -= 1

def fillLowestList(InputList, MotherMap):
    ResultList = []
    for i in InputList:
        ResultList.append(i)
    for i in InputList:
        current = i
        while True:
            if MotherMap[current] != "0":
                current = MotherMap[current]
                if current in ResultList:
                    ResultList.remove(current)
            else:
                break
    return ResultList

out_misplaced.write('Misplaced branches in the pedigree:\n\n')
LowestConflictingUnitsList = fillLowestList(ConflictingUnitsList, MotherMap)
MisplacedUnitsCount = 0
MisplacedBranchesCount = 0
for individual in LowestConflictingUnitsList:
    MisplacedBranchHaplotypedList = []
    current = individual
    MisplacedBranchHaplotypedList.append(current)
    while MotherMap[current] != "0":
        current = MotherMap[current]
        if current in HaplotypedList:
            if HaplotypeMap[current] == HaplotypeMap[individual]:
                MisplacedBranchHaplotypedList.append(current)
                MisplacedUnitsCount += 1
            else:
                break
    if len(MisplacedBranchHaplotypedList) > 1:
        MisplacedBranchesCount += 1
        for i in MisplacedBranchHaplotypedList:
            out_misplaced.write(i + '   haplotype = ' + HaplotypeMap[i] + '\n')
        out_misplaced.write('\n' * 2)
MisplacedUnitsCount += MisplacedBranchesCount
if MisplacedUnitsCount == 0:
    out_misplaced.write('There are no misplaced branches in the pedigree.\n')

out_summary.write('No. of individuals in pedigree' + ' ' * 25 + '%s' % len(IDlist) + '\n')
out_summary.write('No. of distinct haplotypes' + ' ' * 29 + '%s' % len(HaplotypeNamesList) + '\n')
out_summary.write('No. of individuals with sequenced haplotype' + ' ' * 12 + '%s' % len(HaplotypedList) + '\n')
out_summary.write('No. of pairwise mismatches' + ' ' * 29 + '%s' % MismatchCount + '\n')
pruned = len(ConflictingUnitsList) - MisplacedUnitsCount + MisplacedBranchesCount
if len(ConflictingUnitsList) > 0:
    out_summary.write(
        'Maximal no. of conflicts / individual' + ' ' * 18 + '%s' % ConflictCountMap[ConflictingUnitsList[0]] + '/' +
        ConflictingUnitsList[0] + '\n')
else:
    out_summary.write('Maximal no. of conflicts / individual' + ' ' * 18 + '0 / NA\n')
out_summary.write('No. of informative individuals' + ' ' * 25 + '%s' % len(InformativeList) + '\n')
out_summary.write('No. of individuals participating in conflicts' + ' ' * 10 + '%s' % len(UnitsInConflictsList) + '\n')
out_summary.write('No. of conflicting or misplaced individuals' + ' ' * 12 + '%s' % len(ConflictingUnitsList) + '\n')
out_summary.write('Pruned no. of conflicting individuals' + ' ' * 18 + '%s' % pruned + '\n')
out_summary.write('-----------------------------------------------------------\n')
if len(ConflictingUnitsList) > 0:
    HC = len(ConflictingUnitsList) / len(HaplotypedList) * 100
    IC = len(ConflictingUnitsList) / len(InformativeList) * 100
    MISPLACED = (MisplacedUnitsCount - MisplacedBranchesCount) / len(ConflictingUnitsList) * 100
    out_summary.write('HC index = ' + '%s' % HC + '\n' + 'IC index = ' + '%s' % IC + '\n' + 'MISPLACED index = ' + '%s' % MISPLACED + '\n')
else:
    out_summary.write('There are no conflicts in the pedigree.\n')

out_conflicts.close()
out_conflicts_detail.close()
out_summary.close()
out_misplaced.close()
