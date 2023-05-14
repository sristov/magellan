# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

from tkinter import messagebox
import numpy as np
import sys
from operator import countOf

Female_gender = "2"
Male_gender = "1"


def SetIndividualToFounderMap(IDlist, ParentMap, GenderMap, Gender):

    IndividualToFounderMap = {}
    for individual in IDlist:
        current = individual
        while ParentMap[current] != '0':
            current = ParentMap[current]
        # accounting for the case of founder sire in ref pop
        if GenderMap[current] == Gender:
            IndividualToFounderMap[individual] = current
        else:
            IndividualToFounderMap[individual] = '0'

    return IndividualToFounderMap


def MapSamplesInRefPop(ReferencePopulationList,FoundersList,IndividualToFounderMap,GenderMap,HaplotypeMap,Lineage,Gender,mode):

    FounderLineInRefPopList = set()
    FounderDamLineWithOnlyMalesInRefPopList = set()
    LineGenderInRefPopCountMap = {}
    LineHaplotypeMap = {}

    for i in FoundersList:
        if Gender == Gender:
            LineGenderInRefPopCountMap[i] = 0

    for individual in ReferencePopulationList:
        current = IndividualToFounderMap[individual]
        if current == '0':
            continue
        # accounting for the case of founder sire in ref pop
        if GenderMap[current] == Gender:
            if GenderMap[individual] == Gender:
                LineGenderInRefPopCountMap[current] += 1
                FounderLineInRefPopList.add(current)
            if Lineage == 'maternal' and Gender == Female_gender:
                if GenderMap[individual] == Female_gender:
                    FounderDamLineWithOnlyMalesInRefPopList.discard(current)
                else:
                    if current not in FounderLineInRefPopList:
                        FounderDamLineWithOnlyMalesInRefPopList.add(current)
            if not (Lineage == 'maternal' and Gender == Male_gender):
                if individual in HaplotypeMap:
                    if current not in LineHaplotypeMap:
                        LineHaplotypeMap[current] = HaplotypeMap[individual]
        else:
            message = f"ERROR1: GenderMap[{current}] != {Gender}"
            if mode == "gui":
                messagebox.showerror("ERROR", message)
                return
            elif mode == "cl":
                sys.stderr.write(message+'\n')
                sys.exit()


    return FounderLineInRefPopList, FounderDamLineWithOnlyMalesInRefPopList, LineGenderInRefPopCountMap, LineHaplotypeMap

def CountSamplesInRefPop(ReferencePopulationList,
                         HaplotypedList,
                         FounderLineInRefPopList,
                         FounderDamLineWithOnlyMalesInRefPopList,
                         IndividualToFounderMap,
                         HaplotypeMap,
                         Lineage):

    LinesWithSamplesInRefPopHaplotypeMap = {}
    LinesWithSamplesInRefPopCountMap = {}
    for i in FounderLineInRefPopList:
        LinesWithSamplesInRefPopCountMap[i] = 0
    if Lineage == 'maternal':
        for i in FounderDamLineWithOnlyMalesInRefPopList:
            LinesWithSamplesInRefPopCountMap[i] = 0
    for individual in HaplotypedList:
        if individual in ReferencePopulationList:
            if not IndividualToFounderMap[individual] in LinesWithSamplesInRefPopHaplotypeMap:
                LinesWithSamplesInRefPopHaplotypeMap[IndividualToFounderMap[individual]] = HaplotypeMap[individual]
                LinesWithSamplesInRefPopCountMap[IndividualToFounderMap[individual]] += 1
            else:
                LinesWithSamplesInRefPopCountMap[IndividualToFounderMap[individual]] += 1

    return LinesWithSamplesInRefPopHaplotypeMap, LinesWithSamplesInRefPopCountMap


def CalculateEffectivePopSize(NumberOfGenderInRefPop,FoundersLength,FounderLineInRefPopList,LineGenderInRefPopCountMap):
        
        try:
            PI_founders = 1/FoundersLength
        except ZeroDivisionError:
            PI_founders = 0

        PI_line_in_refpop = np.sum([(int(LineGenderInRefPopCountMap[fdam])/NumberOfGenderInRefPop)**2 for fdam in FounderLineInRefPopList])

        try:
            increase_in_identity_founders = (PI_line_in_refpop-PI_founders)/(1-PI_founders)
        except ZeroDivisionError:
            # only one dam line -> division by zero
            increase_in_identity_founders = 9999

        return PI_founders, PI_line_in_refpop, increase_in_identity_founders

def Calculate_PI(FounderList,
                 HaplotypeNamesList,
                 LinesWithSamplesInRefPopHaplotypeMap,
                 LineHaplotypeMap,
                 NumberOfGenderInRefPop,
                 LineGenderInRefPopCountMap):

    LinesWithoutSamplesList = [fdam for fdam in FounderList if fdam not in LinesWithSamplesInRefPopHaplotypeMap]

    SumFounderHapLine = len(LinesWithSamplesInRefPopHaplotypeMap) + len(LinesWithoutSamplesList)

    PI_haplotype_in_founders = 0
    for hapname in HaplotypeNamesList:
        A = len([key for key, value in LineHaplotypeMap.items() if value == hapname])
        PI_haplotype_in_founders += (A/SumFounderHapLine)**2
    PI_haplotype_in_founders += len(LinesWithoutSamplesList)*(1/SumFounderHapLine)**2

    HaplotypeToFounderLinesMap = {}
    for hapname in HaplotypeNamesList:
        AuxLineList = [fdam for fdam in LinesWithSamplesInRefPopHaplotypeMap if LinesWithSamplesInRefPopHaplotypeMap[fdam] == hapname]
        HaplotypeToFounderLinesMap[hapname] = AuxLineList

    PI_haplotype_in_refpop = 0
    increase_in_identity_haplotype = 9999
    if NumberOfGenderInRefPop > 0:
        for hapname in HaplotypeNamesList:
            A = np.sum([LineGenderInRefPopCountMap[fdam] for fdam in HaplotypeToFounderLinesMap[hapname]])
            PI_haplotype_in_refpop += (A/NumberOfGenderInRefPop)**2

        for fdam in LinesWithoutSamplesList:
            PI_haplotype_in_refpop += (LineGenderInRefPopCountMap[fdam]/NumberOfGenderInRefPop)**2

        try:
            increase_in_identity_haplotype = (PI_haplotype_in_refpop-PI_haplotype_in_founders)/(1-PI_haplotype_in_founders)
        except ZeroDivisionError:
            # only one haplotype line -> division by zero
            increase_in_identity_haplotype = 9999

    return PI_haplotype_in_founders, PI_haplotype_in_refpop, increase_in_identity_haplotype

def mag_calc(IDlist,
             FatherMap,
             MotherMap,
             YobMap,
             GenderMap,
             HaplotypeMap,
             HaplotypedList,
             HaplotypeNamesList,
             FirstRefYear,
             LastRefYear,
             Lineage,
             mode):

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
        if int(FirstRefYear) <= int(YobMap[individual]) <= int(LastRefYear):
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

    if Lineage == 'maternal':
        line = 'dam'
    elif Lineage == 'paternal':
        line = 'sire'

    FounderDamsList = [i for i in IDlist if MotherMap[i] == '0' and GenderMap[i] == Female_gender]
    FounderSiresList = [i for i in IDlist if FatherMap[i] == '0' and GenderMap[i] == Male_gender]

    IndividualToFounderDamMap = {}
    IndividualToFounderSireMap = {}


    for individual in IDlist:
        current = individual
        while MotherMap[current] != '0':
            current = MotherMap[current]
        # accounting for the case of founder sire in ref pop
        if GenderMap[current] == Female_gender:
            IndividualToFounderDamMap[individual] = current
        else:
            IndividualToFounderDamMap[individual] = '0'
        current = individual
        while FatherMap[current] != '0':
            current = FatherMap[current]
        # accounting for the case of founder dam in ref pop
        if GenderMap[current] == Male_gender:
            IndividualToFounderSireMap[individual] = current
        else:
            IndividualToFounderSireMap[individual] = '0'

    if Lineage == 'maternal':
        FounderDamLineInRefPopList, \
        FounderDamLineWithOnlyMalesInRefPopList, \
        DamLineFemaleInRefPopCountMap, \
        DamLineHaplotypeMap = MapSamplesInRefPop(ReferencePopulationList,
                                                 FounderDamsList,
                                                 IndividualToFounderDamMap,
                                                 GenderMap,
                                                 HaplotypeMap,
                                                 Lineage,
                                                 Female_gender,
                                                 mode)

        FounderSiresInRefPopList, \
        DummyList1, \
        SireLineInRefPopCountMap, \
        DummyMap1 = MapSamplesInRefPop(ReferencePopulationList,
                                       FounderSiresList,
                                       IndividualToFounderSireMap,
                                       GenderMap,HaplotypeMap,
                                       Lineage,Male_gender,mode)
        SireLineInRefPopCountMap['0'] = 0

        DamLinesWithSamplesInRefPopHaplotypeMap, \
        DamLinesWithSamplesInRefPopCountMap = CountSamplesInRefPop(ReferencePopulationList,
                                                                   HaplotypedList,
                                                                   FounderDamLineInRefPopList,
                                                                   FounderDamLineWithOnlyMalesInRefPopList,
                                                                   IndividualToFounderDamMap,
                                                                   HaplotypeMap,Lineage)

        DummyList1.clear()
        DummyMap1.clear()

        OneSampleCount = countOf(DamLinesWithSamplesInRefPopCountMap.values(), 1)

    elif Lineage == 'paternal':
        FounderSiresInRefPopList, \
        FounderSireLineWithSampledMaleInRefPopList, \
        SireLineAllInRefPopCountMap, \
        SireLineHaplotypeMap = MapSamplesInRefPop(ReferencePopulationList, 
                                                  FounderSiresList,
                                                  IndividualToFounderSireMap,
                                                  GenderMap,
                                                  HaplotypeMap,
                                                  Lineage,Male_gender,mode)

        SireLinesWithSamplesInRefPopHaplotypeMap, \
        SireLinesWithSamplesInRefPopCountMap = CountSamplesInRefPop(ReferencePopulationList,
                                                                    HaplotypedList,FounderSiresInRefPopList,
                                                                    FounderSireLineWithSampledMaleInRefPopList,
                                                                    IndividualToFounderSireMap,
                                                                    HaplotypeMap,Lineage)

        OneSampleCount = countOf(SireLinesWithSamplesInRefPopCountMap.values(), 1)


    SamplesInRefPop = np.intersect1d(HaplotypedList, ReferencePopulationList, assume_unique=True)
    SamplesInRefPopCount = len(SamplesInRefPop)
    if Lineage == 'maternal':
        SamplesInRefPopFemaleCount = len(np.intersect1d(SamplesInRefPop, [key for key, value in GenderMap.items() if value == Female_gender]))
    elif Lineage == 'paternal':
        SamplesInRefPopFemaleCount = 0

    # calculations of the probability of identity and the effective population size

    if Lineage == 'maternal':
        # dam line

        PI_founder_dams, \
        PI_dam_line_in_refpop,\
        increase_in_identity_founder_dams = CalculateEffectivePopSize(len(FemalesInReferencePopulationList),
                                                                      len(FounderDamsList),
                                                                      FounderDamLineInRefPopList, 
                                                                      DamLineFemaleInRefPopCountMap)

        PI_founder_sires, \
        PI_sire_line_in_refpop,\
        increase_in_identity_founder_sires = CalculateEffectivePopSize(len(ReferencePopulationList)-len(FemalesInReferencePopulationList),
                                                                       len(FounderSiresList),
                                                                       FounderSiresInRefPopList,
                                                                       SireLineInRefPopCountMap)


        PI_haplotype_in_founders, \
        PI_haplotype_in_refpop, \
        increase_in_identity_haplotype = Calculate_PI(FounderDamsList,
                                                      HaplotypeNamesList,
                                                      DamLinesWithSamplesInRefPopHaplotypeMap,
                                                      DamLineHaplotypeMap,
                                                      len(FemalesInReferencePopulationList),
                                                      DamLineFemaleInRefPopCountMap)

    elif Lineage == 'paternal':

        PI_founder_sires, \
        PI_sire_line_in_refpop,\
        increase_in_identity_founder_sires = CalculateEffectivePopSize(len(ReferencePopulationList),
                                                                       len(FounderSiresList),
                                                                       FounderSiresInRefPopList,
                                                                       SireLineAllInRefPopCountMap)

        PI_haplotype_in_founders, \
        PI_haplotype_in_refpop, \
        increase_in_identity_haplotype = Calculate_PI(FounderSiresList,
                                                      HaplotypeNamesList,
                                                      SireLinesWithSamplesInRefPopHaplotypeMap,
                                                      SireLineHaplotypeMap,
                                                      len(ReferencePopulationList),
                                                      SireLineAllInRefPopCountMap)


    with open('OutputCalc_InputAndResults.txt', 'w') as output1:

        output1.write(f'Records in the studbook = {len(IDlist)}\n')
        output1.write(f'No. of individuals in the reference population ({FirstRefYear} - {LastRefYear}) = {len(ReferencePopulationList)}\n')
        output1.write(f'No. of female individuals in the reference population = {len(FemalesInReferencePopulationList)}\n')
        output1.write(f'No. of founder dams = {len(FounderDamsList)}\n')
        if Lineage == 'maternal':
            output1.write(f'No. of founder dam lines in reference population = {len(FounderDamLineInRefPopList)}\n')
            output1.write(f'No. of founder dam lines including lines with only males in reference population')
            output1.write(f' = {(len(FounderDamLineInRefPopList)+len(FounderDamLineWithOnlyMalesInRefPopList))}\n')
            output1.write(f'No. of founder {line} lines in reference population with samples = {len(DamLinesWithSamplesInRefPopHaplotypeMap)}\n')
        elif Lineage == 'paternal':
            output1.write(f'No. of founder sire lines in reference population = {len(FounderSiresInRefPopList)}\n')
            output1.write(f'No. of founder sire lines in reference population with samples = {len(SireLinesWithSamplesInRefPopHaplotypeMap)}\n')
        output1.write(f'No. of founder {line} lines in reference population with only one sample = {OneSampleCount}\n')
        output1.write(f'Total number of samples (haplotyped individuals) = {len(HaplotypedList)}\n')
        output1.write(f'No. of samples in reference population = {SamplesInRefPopCount}:   ')
        output1.write(f'dams = {SamplesInRefPopFemaleCount},  sires = {SamplesInRefPopCount-SamplesInRefPopFemaleCount}\n')
        output1.write('----------------------------------------------------------------\n\n')
        if Lineage == 'maternal':
            output1.write('Dam lines:\n\n')
            output1.write(f'Probability of identity in founder dams = {PI_founder_dams}\n')
            output1.write(f'Probability of identity of a dam line in reference population = {PI_dam_line_in_refpop}\n')
            if increase_in_identity_founder_dams == 9999:
                output1.write('Increase in identity = N/A\n')
                output1.write('Effective dam line size = N/A \n\n')
            else:
                output1.write(f'Increase in identity = {increase_in_identity_founder_dams}\n')
                output1.write(f'Effective dam line size = {(1/increase_in_identity_founder_dams)}\n\n')
            output1.write('Haplotype lines:\n\n')
            output1.write(f'Probability of identity of a haplotype line in founder population = {PI_haplotype_in_founders}\n')
            output1.write(f'Probability of identity of a haplotype line in reference population = {PI_haplotype_in_refpop}\n')
            if increase_in_identity_haplotype == 9999:
                output1.write('Increase in identity = N/A\n')
                output1.write('Effective haplotype line size = N/A \n\n')
            else:
                output1.write(f'Increase in identity = {increase_in_identity_haplotype}\n')
                output1.write(f'Effective haplotype line size = {(1/increase_in_identity_haplotype)}\n\n')

        output1.write('*Sire lines:\n\n')
        output1.write(f'No. of founder sires = {len(FounderSiresList)}\n')
        output1.write(
            f'No. of founder sire lines in reference population = {len(FounderSiresInRefPopList)}\n\n')
        output1.write(
            f'Probability of identity in founder sires = {PI_founder_sires}\n')
        output1.write(
            f'Probability of identity of a sire line in reference population = {PI_sire_line_in_refpop}\n')
        if increase_in_identity_founder_sires == 9999:
            output1.write('Increase in identity = N/A\n')
            output1.write('Effective sire line size = N/A \n\n')
        else:
            output1.write(f'Increase in identity = {increase_in_identity_founder_sires}\n')
            output1.write(f'Effective sire line size = {(1/increase_in_identity_founder_sires)}\n\n')
        if Lineage == 'paternal':
            output1.write('Haplotype lines:\n\n')
            output1.write(f'Probability of identity of a haplotype line in founder population = {PI_haplotype_in_founders}\n')
            output1.write(f'Probability of identity of a haplotype line in reference population = {PI_haplotype_in_refpop}\n')
            if increase_in_identity_haplotype == 9999:
                output1.write('Increase in identity = N/A\n')
                output1.write('Effective haplotype line size = N/A \n\n')
            else:
                output1.write(f'Increase in identity = {increase_in_identity_haplotype}\n')
                output1.write(f'Effective haplotype line size = {(1/increase_in_identity_haplotype)}\n\n')

    return
