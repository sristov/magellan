# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

import sys
import numpy as np
from tkinter import messagebox

def mag_stat(IDlist,
             ParentMap,
             YobMap,
             GenderMap,
             HaplotypeMap,
             HaplotypedList,
             FirstRefYear,
             LastRefYear,
             Lineage,
             mode):

    Female_gender = "2"
    Male_gender = "1"

    if Lineage == 'paternal':
        Gender = Male_gender
    elif Lineage == 'maternal':
        Gender = Female_gender

    all_yob = np.unique(list(YobMap.values()))
    all_yob = np.delete(all_yob, np.where(all_yob == 'MISSING_YEAR'))
    all_yob = all_yob.astype(int)

    if FirstRefYear == 0:
        FirstRefYear = int(min(all_yob))
    if LastRefYear == 0:
        LastRefYear = int(max(all_yob))

    ReferencePopulationList = []
    FemalesInReferencePopulationList = []
    for individual in IDlist:
        if YobMap[individual] == 'MISSING_YEAR':
            continue
        if FirstRefYear <= int(YobMap[individual]) <= LastRefYear:
            if Lineage == 'paternal':
                if GenderMap[individual] == Male_gender:
                    ReferencePopulationList.append(individual)
            elif Lineage == 'maternal':
                ReferencePopulationList.append(individual)
                if GenderMap[individual] == Female_gender:
                    FemalesInReferencePopulationList.append(individual)

    if len(ReferencePopulationList) == 0:
        if mode == "gui":
            messagebox.showerror("ERROR", "ERROR: There are no individuals in the reference population!")
            return
        elif mode == "cl":
            sys.stderr.write("ERROR: There are no individuals in the reference population!\n")
            sys.exit()

    FounderList = [i for i in IDlist if (ParentMap[i] == '0') and GenderMap[i] == Gender]

    FounderLineInRefPopList = set()
    FounderLineWithSampledInRefPopList = set()
    FounderDamLineWithOnlyMalesInRefPopList = set()
    IndividualToFounderMap = {}
    LineAllInRefPopCountMap = {}
    DamLineFemaleInRefPopCountMap = {}
    LineHaplotypeMap = {}
    IndividualsMissingFounderList = []
    IndividualsMissingFounderInRefPopCount = 0

    for individual in IDlist:
        current = individual
        try:
            while ParentMap[current] != '0':
                current = ParentMap[current]
        except KeyError:
            pass

        try:
            # accounting for the case of founder sire
            if GenderMap[current] == Gender:  
                IndividualToFounderMap[individual] = current
            else:
                IndividualToFounderMap[individual] = '0'
                IndividualsMissingFounderList.append(individual)
                if individual in ReferencePopulationList:
                    IndividualsMissingFounderInRefPopCount += 1
        except KeyError:
            IndividualToFounderMap[individual] = current
            IndividualsMissingFounderList.append(individual)
            if individual in ReferencePopulationList:
                IndividualsMissingFounderInRefPopCount += 1
            pass

    for i in FounderList:
        LineAllInRefPopCountMap[i] = 0
        if Lineage == 'maternal':
            DamLineFemaleInRefPopCountMap[i] = 0

    for individual in ReferencePopulationList:
        current = IndividualToFounderMap[individual]
        if current == '0':
            continue
        LineAllInRefPopCountMap[current] += 1
        if GenderMap[individual] == Gender:
            if Lineage == 'maternal':
                DamLineFemaleInRefPopCountMap[current] += 1
            FounderLineInRefPopList.add(current)
            if Lineage == 'maternal':
                FounderDamLineWithOnlyMalesInRefPopList.discard(current)
        else:
            if current not in FounderLineInRefPopList:
                FounderDamLineWithOnlyMalesInRefPopList.add(current)
        if individual in HaplotypeMap:
            if current not in LineHaplotypeMap:
                LineHaplotypeMap[current] = HaplotypeMap[individual]
            if GenderMap[individual] == Gender:
                FounderLineWithSampledInRefPopList.add(current)

    LinesWithSamplesInRefPopHaplotypeMap = {}
    LinesWithSamplesInRefPopCountMap = {}
    DamLinesWithFemaleSamplesInRefPopCountMap = {}
    for i in FounderLineInRefPopList:
        LinesWithSamplesInRefPopCountMap[i] = 0
        if Lineage == 'maternal':
            DamLinesWithFemaleSamplesInRefPopCountMap[i] = 0
    for i in FounderDamLineWithOnlyMalesInRefPopList:
        LinesWithSamplesInRefPopCountMap[i] = 0
    for individual in HaplotypedList:
        if individual in ReferencePopulationList:
            if IndividualToFounderMap[individual] not in LinesWithSamplesInRefPopHaplotypeMap:
                LinesWithSamplesInRefPopHaplotypeMap[IndividualToFounderMap[individual]] = HaplotypeMap[individual]
                LinesWithSamplesInRefPopCountMap[IndividualToFounderMap[individual]] += 1
                if GenderMap[individual] == Female_gender:
                    DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderMap[individual]] += 1
            else:
                LinesWithSamplesInRefPopCountMap[IndividualToFounderMap[individual]] += 1
                if GenderMap[individual] == Female_gender:
                    DamLinesWithFemaleSamplesInRefPopCountMap[IndividualToFounderMap[individual]] += 1

    SamplesInRefPopCount = 0
    SamplesInRefPopFemaleCount = 0
    for individual in HaplotypedList:
        if individual in ReferencePopulationList:
            SamplesInRefPopCount += 1
            if GenderMap[individual] == Female_gender:
                SamplesInRefPopFemaleCount += 1

    LineMembershipCountMap = {}
    for i in FounderList:
        LineMembershipCountMap[i] = 0
    for individual in IDlist:
        if IndividualToFounderMap[individual] != '0':
            LineMembershipCountMap[IndividualToFounderMap[individual]] += 1

    if Lineage == 'maternal':
        line = 'dam'
        with open('OutputStat_DamLinesWithFemalesInRefPop.txt', 'w') as damlineout1, \
             open('OutputStat_DamLinesWithOnlyMalesInRefPop.txt', 'w') as damlineout2:

            for fdam in FounderLineInRefPopList:
                damlineout1.write(f'{fdam}: no. of desc. in ref. pop. = {LineAllInRefPopCountMap[fdam]}\n')
                if fdam in LinesWithSamplesInRefPopHaplotypeMap:
                    damlineout1.write(f'; no. of samples = {LinesWithSamplesInRefPopCountMap[fdam]}:')
                    damlineout1.write(f' F = {DamLinesWithFemaleSamplesInRefPopCountMap[fdam]},; haplotype = {LinesWithSamplesInRefPopHaplotypeMap[fdam]}')
                    damlineout1.write(f' M = {(LinesWithSamplesInRefPopCountMap[fdam]-DamLinesWithFemaleSamplesInRefPopCountMap[fdam])};')
                    damlineout1.write(f' haplotype = {LinesWithSamplesInRefPopHaplotypeMap[fdam]}\n')
                else:
                    damlineout1.write('\n')
 
            if len(FounderDamLineWithOnlyMalesInRefPopList) == 0:
                damlineout2.write('There are no founder dams with only male descendants in reference population.\n')
            for fdam in FounderDamLineWithOnlyMalesInRefPopList:
                damlineout2.write(f'{fdam}: no. of male desc. in ref. pop. = {LineAllInRefPopCountMap[fdam]}')
                if fdam in LinesWithSamplesInRefPopHaplotypeMap:
                    damlineout1.write('; no. of samples = {DamLinesWithSamplesInRefPopCountMap[fdam]}; haplotype = {LinesWithSamplesInRefPopHaplotypeMap[fdam]}\n')
                else:
                    damlineout1.write('\n')
    elif Lineage == 'paternal':
        line = 'sire'
        with open('OutputStat_SireLinesWithMalesInRefPop.txt', 'w') as sirelineout1:

            for fsire in FounderLineInRefPopList:
                sirelineout1.write(f'{fsire}: no. of desc. in ref. pop. = {LineAllInRefPopCountMap[fsire]}\n')
            if fsire in LinesWithSamplesInRefPopHaplotypeMap:
                sirelineout1.write(f'; no. of samples = {LinesWithSamplesInRefPopCountMap[fsire]}')
                sirelineout1.write(f'; haplotype = {LinesWithSamplesInRefPopHaplotypeMap[fsire]}\n')
            else:
                sirelineout1.write('\n')

    with open(f'OutputStat_{line.capitalize()}LineMembershipAllInRefPop.txt', 'w') as damlineout3, \
         open(f'OutputStat_{line.capitalize()}LineMembership_1.txt', 'w') as damlineout5, \
         open(f'OutputStat_{line.capitalize()}LineMembership_2.txt', 'w') as damlineout6:

        damlineout3.write(f'founder {line}:individual in dam line and in reference population\n')
        damlineout6.write(f'founder {line}:individual in {line} line\n')

        for fdam in sorted(LineMembershipCountMap, key=LineMembershipCountMap.__getitem__, reverse=True):
            if fdam in LineHaplotypeMap:
                damlineout5.write(f"founder {line}: {fdam}     number of individuals in {line} line = {LineMembershipCountMap[fdam]}       haplotype = {LineHaplotypeMap[fdam]}\n")
            else:
                damlineout5.write(f"founder {line}: {fdam}     number of individuals in {line} line = {LineMembershipCountMap[fdam]}       haplotype = N/A\n")

            offspring = [key for key, value in IndividualToFounderMap.items() if value == fdam]
            for individual in offspring:
                damlineout5.write(f"{individual}\n")
                damlineout6.write(f"{fdam}:{individual}\n")
                if individual in ReferencePopulationList:
                    damlineout3.write(f'{fdam}:{individual}\n')

    if Lineage == 'maternal':
        with open('OutputStat_DamLineMembershipFemaleOnlyInRefPop.txt', 'w') as damlineout4:
            damlineout4.write('founder dam:female individual in dam line and in reference population\n')
            for fdam in sorted(LineMembershipCountMap, key=LineMembershipCountMap.__getitem__, reverse=True):
                offspring = [key for key, value in IndividualToFounderMap.items() if value == fdam]
                for individual in offspring:
                    if (individual in ReferencePopulationList) and (GenderMap[individual] == Female_gender):
                        damlineout4.write(f'{fdam}:{individual}\n')
    return