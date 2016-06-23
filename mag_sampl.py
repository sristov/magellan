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
HowManyToSequence = 100  # default number of sequencings; used if planned_number_of_sequencings.txt file is not present

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
AvailableMap = {}

HAP_column = -1
AVAILABLE_column = -1

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
    if lineparts[i] == 'available':
        AVAILABLE_column = i

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
    if AVAILABLE_column != -1:
        if lineparts[AVAILABLE_column] == '1':
            AvailableMap[lineparts[ID_column]] = '1'
        else:
            AvailableMap[lineparts[ID_column]] = '0'

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

# if available data not provided: set available = ref pop
if AVAILABLE_column == -1:
    for individual in ReferencePopulationList:
        AvailableMap[individual] = '1'

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
    if GenderMap[current] == Female_gender:  # accounting for the case of founder sire in ref pop
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
        if IndividualToFounderDamMap[individual] not in DamLinesWithSamplesInRefPopHaplotypeMap:
            DamLinesWithSamplesInRefPopHaplotypeMap[IndividualToFounderDamMap[individual]] = HaplotypeMap[individual]
            DamLinesWithSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
            if GenderMap[individual] == Female_gender:
                DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
        else:
            DamLinesWithSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1
            if GenderMap[individual] == Female_gender:
                DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderDamMap[individual]] += 1

# calculating the target number of samples per dam line

if os.path.isfile('planned_number_of_sequencings.txt'):
    input3 = open('planned_number_of_sequencings.txt', 'r')
    HowManyToSequence = int(input3.readline())

TargetPerDamLineMap = {}
RemainingToDoPerDamLineMap = {}
RefPopSize = int(len(ReferencePopulationList))
AdjustedRefPopSize = RefPopSize - IndividualsMissingFounderDamInRefPopCount
PreviouslySequencedInRefPop = 0
Remainder1map = {}

for dam in DamLinesWithSamplesInRefPopCountMap:
    PreviouslySequencedInRefPop += DamLinesWithSamplesInRefPopCountMap[dam]
MixedList = FounderDamLineInRefPopList + FounderDamLineWithOnlyMalesInRefPopList
NoPreviousSamplesCount = 0
for dam in MixedList:
    if DamLinesWithSamplesInRefPopCountMap[dam] == 0:
        NoPreviousSamplesCount += 1
MixedCountMap = {}
for dam in MixedList:
    MixedCountMap[dam] = DamLineAllInRefPopCountMap[dam]

if NoPreviousSamplesCount >= HowManyToSequence:
    IntermediateSum = HowManyToSequence
    for dam in sorted(DamLineAllInRefPopCountMap, key=DamLineAllInRefPopCountMap.__getitem__, reverse=True):
        if DamLinesWithSamplesInRefPopCountMap[dam] == 0:
            RemainingToDoPerDamLineMap[dam] = 1
            IntermediateSum -= 1
        if IntermediateSum == 0:
            break
    TargetSum = 0
else:
    IntermediateSum = HowManyToSequence - NoPreviousSamplesCount
    TargetSum = PreviouslySequencedInRefPop + HowManyToSequence - len(MixedList)
    for dam in sorted(MixedCountMap, key=MixedCountMap.__getitem__, reverse=True):
        if IntermediateSum > 0:
            TargetPerDamLineMap[dam] = TargetSum * (DamLineAllInRefPopCountMap[dam]-1)/(AdjustedRefPopSize-len(MixedList))
            Remainder1map[dam] = TargetPerDamLineMap[dam] - int(TargetPerDamLineMap[dam])
            TargetPerDamLineMap[dam] = int(TargetPerDamLineMap[dam])
            TargetPerDamLineMap[dam] += 1
            if DamLinesWithSamplesInRefPopCountMap[dam] > TargetPerDamLineMap[dam]:
                RemainingToDoPerDamLineMap[dam] = 0
                Remainder1map[dam] = 0
            else:
                if TargetPerDamLineMap[dam] > DamLineAllInRefPopCountMap[dam]:
                    TargetPerDamLineMap[dam] = DamLineAllInRefPopCountMap[dam]
                RemainingToDoPerDamLineMap[dam] = TargetPerDamLineMap[dam] - DamLinesWithSamplesInRefPopCountMap[dam]
                if DamLinesWithSamplesInRefPopCountMap[dam] > 0:
                    IntermediateSum -= RemainingToDoPerDamLineMap[dam]
                else:
                    IntermediateSum -= (RemainingToDoPerDamLineMap[dam] - 1)
                if IntermediateSum < 0:
                    RemainingToDoPerDamLineMap[dam] += IntermediateSum
                    IntermediateSum = 0
        else:
            if DamLinesWithSamplesInRefPopCountMap[dam] == 0:
                if dam in RemainingToDoPerDamLineMap:
                    if RemainingToDoPerDamLineMap[dam] == 0:
                        RemainingToDoPerDamLineMap[dam] = 1
                else:
                    RemainingToDoPerDamLineMap[dam] = 1

for dam in sorted(Remainder1map, key=Remainder1map.__getitem__, reverse=True):
    if IntermediateSum == 0:
        break
    if DamLineAllInRefPopCountMap[dam] > (DamLinesWithSamplesInRefPopCountMap[dam] + RemainingToDoPerDamLineMap[dam]):
        TargetPerDamLineMap[dam] += 1
        if DamLinesWithSamplesInRefPopCountMap[dam] > TargetPerDamLineMap[dam]:
            continue
        RemainingToDoPerDamLineMap[dam] += 1
        IntermediateSum -= 1

nonassigned = 0
if IntermediateSum > 0:
    count = 0
    for i in IndividualsMissingFounderDamList:
        if i in ReferencePopulationList:
            count += 1
    if IntermediateSum > count:
        print('\nNot enough individuals for such an ambitious sampling plan!\n')
    if count > IntermediateSum:
        nonassigned = IntermediateSum
    else:
        nonassigned = count
    for dam in MixedList:
        if (RemainingToDoPerDamLineMap[dam] + DamLinesWithSamplesInRefPopCountMap[dam]) > DamLineAllInRefPopCountMap[dam]:
            RemainingToDoPerDamLineMap[dam] = DamLineAllInRefPopCountMap[dam] - DamLinesWithSamplesInRefPopCountMap[dam]

output1 = open('OutputSampl_DetailedInfo.txt', 'w')
output1.write('reference population size = %s\ndam lines in ref.pop. = %s\n' % (RefPopSize, len(MixedList)))
output1.write('IndividualsMissingFounderDamInRefPopCount = %s\n' % IndividualsMissingFounderDamInRefPopCount)
output1.write('PreviouslySequencedInRefPop = %s\nHowManyToSequence = %s\n' % (PreviouslySequencedInRefPop, HowManyToSequence))
output1.write('Number of lines without previous samples = %s\n' % NoPreviousSamplesCount)
output1.write('TargetSum (that remains to be proportionally divided among dam lines) = %s\n\n' % TargetSum)
output1.write('FOUNDER DAM:\n')
for dam in MixedList:
    if dam not in TargetPerDamLineMap:
        TargetPerDamLineMap[dam] = 0
    if dam not in RemainingToDoPerDamLineMap:
        RemainingToDoPerDamLineMap[dam] = 0
    output1.write(dam + ' ->  in ref.pop.: %s   targ: %s   prev: %s  todo:%s\n' % (DamLineAllInRefPopCountMap[dam], TargetPerDamLineMap[dam], DamLinesWithSamplesInRefPopCountMap[dam], RemainingToDoPerDamLineMap[dam]))
if nonassigned > 0:
    output1.write('\nThere are %s remaining non-assigned planned samplings\n' % IntermediateSum)
    output1.write('The list of %s individuals in reference population that can be sampled but that do not have the founder dam:\n' % count)
    for individual in IndividualsMissingFounderDamList:
        if individual in ReferencePopulationList:
            output1.write(individual + '\n')

# calculating the available number of samples per dam line (out of the target number of samples)

AvailablePerDamLineCountMap = {}
RemainingDiffAvailableMap = {}
AvailableAllCount = 0
AvailableRealCount = 0
AvailableDifferenceCount = 0
RemainingToDoAvailablePerDamLineMap = {}
if AVAILABLE_column != -1:
    for dam in MixedList:
        AvailablePerDamLineCountMap[dam] = 0
        RemainingToDoAvailablePerDamLineMap[dam] = RemainingToDoPerDamLineMap[dam]
    for individual in ReferencePopulationList:
        if AvailableMap[individual] == '1':
            AvailableAllCount += 1
            if IndividualToFounderDamMap[individual] != '0':
                if individual not in HaplotypedList:
                    AvailableRealCount += 1
                    AvailablePerDamLineCountMap[IndividualToFounderDamMap[individual]] += 1
    for dam in MixedList:
        RemainingDiffAvailableMap[dam] = None
        if AvailablePerDamLineCountMap[dam] < RemainingToDoPerDamLineMap[dam]:
            RemainingDiffAvailableMap[dam] = RemainingToDoPerDamLineMap[dam] - AvailablePerDamLineCountMap[dam]
            AvailableDifferenceCount += RemainingDiffAvailableMap[dam]
            RemainingToDoAvailablePerDamLineMap[dam] = AvailablePerDamLineCountMap[dam]
else:
    for dam in MixedList:
        RemainingToDoAvailablePerDamLineMap[dam] = RemainingToDoPerDamLineMap[dam]

if AvailableDifferenceCount > 0:
    output2 = open('OutputSampl_AvailabilityRestrictions.txt', 'w')
    output2.write('reference population size = %s    dam lines in ref. pop. = %s    PreviouslySequencedInRefPop = %s    HowManyToSequence = %s\n\n' % (RefPopSize, len(MixedList), PreviouslySequencedInRefPop, HowManyToSequence))
    output2.write('AvailableAllCount = %s    AvailableRealCount = %s    sum of DIFF (planned - available) = %s    total of both available and planned: %s\n\n' % (AvailableAllCount, AvailableRealCount, AvailableDifferenceCount, HowManyToSequence-AvailableDifferenceCount))
    for dam in sorted(MixedCountMap, key=MixedCountMap.__getitem__, reverse=True):
        output2.write('dam: %s   in ref.pop: %s  targ: %s  prev: %s  todo: %s  avail: %s  DIFF: %s  A_TODO: %s' % (dam, DamLineAllInRefPopCountMap[dam], TargetPerDamLineMap[dam], DamLinesWithSamplesInRefPopCountMap[dam], RemainingToDoPerDamLineMap[dam], AvailablePerDamLineCountMap[dam], RemainingDiffAvailableMap[dam], RemainingToDoAvailablePerDamLineMap[dam]))
        if AvailablePerDamLineCountMap[dam] > RemainingToDoAvailablePerDamLineMap[dam]:
            output2.write('  REMAINING AVAILABLE: %s\n' % (AvailablePerDamLineCountMap[dam] - RemainingToDoAvailablePerDamLineMap[dam]))
        else:
            output2.write('\n')


# determining individuals for sampling

def distanceInFounderDamLine(first_id, second_id):
    MotherLine1 = []
    MotherLine2 = []
    ID_string = first_id
    MotherLine1.append(ID_string)
    while MotherMap[ID_string] != '0':
        ID_string = MotherMap[ID_string]
        MotherLine1.append(ID_string)
    ID_string = second_id
    MotherLine2.append(ID_string)
    while MotherMap[ID_string] != '0':
        ID_string = MotherMap[ID_string]
        MotherLine2.append(ID_string)
    for k in range(len(MotherLine1)):
        for l in range(len(MotherLine2)):
            if MotherLine2[l] == MotherLine1[k]:
                return k + l

output3 = open('OutputSampl_IndividualsForSampling.txt', 'w')
for dam in MixedList:
    CountDown = RemainingToDoAvailablePerDamLineMap[dam]
    if CountDown == 0:
        continue
    CandidateList = []
    PreviousList = []
    for individual in ReferencePopulationList:
        if IndividualToFounderDamMap[individual] == dam:
            if AvailableMap[individual] == '1':
                if individual in HaplotypedList:
                    PreviousList.append(individual)
                else:
                    CandidateList.append(individual)
    if len(PreviousList) == 0:
        Distance = 0
        TopCandidate = None
        for individual in CandidateList:
            SumDistance = 0
            for ind_prev in CandidateList:
                SumDistance += distanceInFounderDamLine(individual, ind_prev)
            if Distance != 0:
                if SumDistance <= Distance:
                    Distance = SumDistance
                    TopCandidate = individual
            else:
                TopCandidate = individual
        PreviousList.append(TopCandidate)
        CandidateList.remove(TopCandidate)
        CountDown -= 1
    while CountDown > 0:
        Distance = 0
        TopCandidate = None
        for individual in CandidateList:
            SumDistance = 0
            for ind_prev in PreviousList:
                SumDistance += distanceInFounderDamLine(individual, ind_prev)
            if SumDistance >= Distance:
                Distance = SumDistance
                TopCandidate = individual
        PreviousList.append(TopCandidate)
        CandidateList.remove(TopCandidate)
        CountDown -= 1
    output3.write('founder dam %s : ' % dam)
    for individual in PreviousList:
        if individual not in HaplotypedList:
            output3.write(individual + ', ')
    output3.write('\n\n')


