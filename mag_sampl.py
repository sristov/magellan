# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

import sys
import time
import numpy as np
from tkinter import messagebox
from fast_sampler import fastSampler

Female_gender = "2"
Male_gender = "1"

class MagSampl():

    def __init__(self,
                 IDlist,
                 ParentMap,
                 YobMap,
                 GenderMap,
                 HaplotypeMap,
                 HaplotypedList,
                 HaplotypeNamesList,
                 AvailableMap,
                 Lineage):

        self.IDlist = IDlist
        self.Lineage = Lineage
        self.ParentMap = ParentMap
        if self.Lineage == 'paternal':
           self.Gender = Male_gender
        elif self.Lineage == 'maternal':
           self.Gender = Female_gender
        
        self.YobMap = YobMap
        self.GenderMap = GenderMap
        self.HaplotypeMap = HaplotypeMap
        self.HaplotypedList = HaplotypedList
        self.HaplotypeNamesList = HaplotypeNamesList
        self.AvailableMap = AvailableMap
        

        pass

    def setRefPop(self,FirstRefYear,LastRefYear,mode):

        all_yob = np.unique(list(self.YobMap.values()))
        all_yob = np.delete(all_yob, np.where(all_yob == 'MISSING_YEAR'))
        all_yob = all_yob.astype(int)

        if FirstRefYear == 0:
            FirstRefYear = int(min(all_yob))
        if LastRefYear == 0:
            LastRefYear = int(max(all_yob))

        self.ReferencePopulationList = []
        self.FemalesInReferencePopulationList = []
        for individual in self.IDlist:
            if self.YobMap[individual] == 'MISSING_YEAR':
                continue
            if int(FirstRefYear) <= int(self.YobMap[individual]) <= int(LastRefYear):
                if self.Lineage == 'paternal':
                    if self.GenderMap[individual] == Male_gender:
                        self.ReferencePopulationList.append(individual)
                elif self.Lineage == 'maternal':
                    self.ReferencePopulationList.append(individual)
                    if self.GenderMap[individual] == Female_gender:
                        self.FemalesInReferencePopulationList.append(individual)

        if len(self.ReferencePopulationList) == 0:
            if mode == "gui":
                messagebox.showerror("ERROR", "ERROR: There are no individuals in the reference population!")
                return 1
            elif mode == "cl":
                sys.stderr.write("ERROR: There are no individuals in the reference population!\n")
                sys.exit()

        AVAILABLE_column = 0
        # if available data not provided: set available = ref pop
        if len(self.AvailableMap) == 0:
            AVAILABLE_column = -1
            for individual in self.ReferencePopulationList:
                self.AvailableMap[individual] = '1'

        return AVAILABLE_column

    def setFounderPop(self,mode):

        self.FounderList = [individual for individual in self.IDlist if self.ParentMap[individual] == '0' and self.GenderMap[individual] == self.Gender]

        self.FounderLineInRefPopList = set()
        self.FounderLineWithSampledInRefPopList = set()
        self.FounderDamLineWithOnlyMalesInRefPopList = set()
        self.IndividualToFounderMap = {}
        self.LineAllInPopCountMap = {}
        self.LineAllInRefPopCountMap = {}
        self.DamLineFemaleInRefPopCountMap = {}
        self.LineHaplotypeMap = {}
        self.IndividualsMissingFounderList = []
        self.IndividualsMissingFounderInRefPopCount = 0

        for individual in self.IDlist:
            current = individual
            while self.ParentMap[current] != '0':
                current = self.ParentMap[current]
            if self.GenderMap[current] == self.Gender:  # accounting for the case of founder sire in ref pop
                self.IndividualToFounderMap[individual] = current
            else:
                self.IndividualToFounderMap[individual] = '0'
                self.IndividualsMissingFounderList.append(individual)
                if individual in self.ReferencePopulationList:
                    self.IndividualsMissingFounderInRefPopCount += 1

        for i in self.FounderList:
            self.LineAllInPopCountMap[i] = 0
            self.LineAllInRefPopCountMap[i] = 0
            if self.Lineage == 'maternal':
                self.DamLineFemaleInRefPopCountMap[i] = 0

        # Count total population in each line
        for individual in self.IDlist:
            current = self.IndividualToFounderMap[individual]
            if current == '0':
                continue
            self.LineAllInPopCountMap[current] += 1

        for individual in self.ReferencePopulationList:
            current = self.IndividualToFounderMap[individual]
            if current == '0':
               continue
            # accounting for the case of founder sire in ref pop
            if self.GenderMap[current] == self.Gender: 
                self.LineAllInRefPopCountMap[current] += 1
                if self.GenderMap[individual] == self.Gender:
                    if self.Lineage == 'maternal':
                        self.DamLineFemaleInRefPopCountMap[current] += 1
                    self.FounderLineInRefPopList.add(current)
                    if self.Lineage == 'maternal':
                        self.FounderDamLineWithOnlyMalesInRefPopList.discard(current)
                else:
                    if current not in self.FounderLineInRefPopList:
                        self.FounderDamLineWithOnlyMalesInRefPopList.add(current)
                if individual in self.HaplotypeMap:
                    if current not in self.LineHaplotypeMap:
                        self.LineHaplotypeMap[current] = self.HaplotypeMap[individual]
                    if self.GenderMap[individual] == self.Gender:
                        self.FounderLineWithSampledInRefPopList.add(current)
            else:
                if mode == "gui":
                    messagebox.showerror(f"ERROR", "ERROR: GenderMap[{current}] != {self.Gender}!")
                    return
                elif mode == "cl":
                    sys.stderr.write("ERROR: GenderMap[{current}] != {self.Gender}!\n")
                    sys.exit()

        self.LinesWithSamplesInRefPopHaplotypeMap = {}
        self.LinesWithSamplesInRefPopCountMap = {}
        self.DamLinesWithFemaleSamplesInRefPopCountMap = {}

        for i in self.FounderLineInRefPopList:
            self.LinesWithSamplesInRefPopCountMap[i] = 0
            if self.Lineage == 'maternal':
                self.DamLinesWithFemaleSamplesInRefPopCountMap[i] = 0
        for i in self.FounderDamLineWithOnlyMalesInRefPopList:
            self.LinesWithSamplesInRefPopCountMap[i] = 0
        for individual in self.HaplotypedList:
            if individual in self.ReferencePopulationList:
                if self.IndividualToFounderMap[individual] not in self.LinesWithSamplesInRefPopHaplotypeMap:
                    self.LinesWithSamplesInRefPopHaplotypeMap[self.IndividualToFounderMap[individual]] = self.HaplotypeMap[individual]
                    self.LinesWithSamplesInRefPopCountMap[self.IndividualToFounderMap[individual]] += 1
                    if self.GenderMap[individual] == Female_gender:
                        self.DamLinesWithFemaleSamplesInRefPopCountMap[self.IndividualToFounderMap[individual]] += 1
                else:
                    self.LinesWithSamplesInRefPopCountMap[self.IndividualToFounderMap[individual]] += 1
                    if self.GenderMap[individual] == Female_gender:
                        self.DamLinesWithFemaleSamplesInRefPopCountMap[self.IndividualToFounderMap[individual]] += 1

        pass

    def calcTargetNumSamplesPerLine(self,HowManyToSequence,mode):

        self.TargetPerLineMap = {}
        self.RemainingToDoPerLineMap = {}
        self.RefPopSize = int(len(self.ReferencePopulationList))
        self.AdjustedRefPopSize = self.RefPopSize - self.IndividualsMissingFounderInRefPopCount
        self.PreviouslySequencedInRefPop = 0
        self.Remainder1map = {}

        for individual in self.LinesWithSamplesInRefPopCountMap:
            self.PreviouslySequencedInRefPop += self.LinesWithSamplesInRefPopCountMap[individual]
        self.MixedList = self.FounderLineInRefPopList.union(self.FounderDamLineWithOnlyMalesInRefPopList)
        self.NoPreviousSamplesCount = 0
        for individual in self.MixedList:
            if self.LinesWithSamplesInRefPopCountMap[individual] == 0:
                self.NoPreviousSamplesCount += 1
        self.MixedCountMap = {}
        for individual in self.MixedList:
            self.MixedCountMap[individual] = self.LineAllInRefPopCountMap[individual]

        if self.NoPreviousSamplesCount >= int(HowManyToSequence):
            self.IntermediateSum = int(HowManyToSequence)
            for individual in sorted(self.LineAllInRefPopCountMap, key=self.LineAllInRefPopCountMap.__getitem__, reverse=True):
                if self.LinesWithSamplesInRefPopCountMap[individual] == 0:
                    self.RemainingToDoPerLineMap[individual] = 1
                    self.IntermediateSum -= 1
                if self.IntermediateSum == 0:
                    break
            self.TargetSum = 0
        else:
            self.IntermediateSum = int(HowManyToSequence) - self.NoPreviousSamplesCount
            self.TargetSum = self.PreviouslySequencedInRefPop + int(HowManyToSequence) - len(self.MixedList)
            for individual in sorted(self.MixedCountMap, key=self.MixedCountMap.__getitem__, reverse=True):
                if self.IntermediateSum > 0:
                    self.TargetPerLineMap[individual] = self.TargetSum * (self.LineAllInRefPopCountMap[individual]-1)/(self.AdjustedRefPopSize-len(self.MixedList))
                    self.Remainder1map[individual] = self.TargetPerLineMap[individual] - int(self.TargetPerLineMap[individual])
                    self.TargetPerLineMap[individual] = int(self.TargetPerLineMap[individual])
                    self.TargetPerLineMap[individual] += 1
                    if self.LinesWithSamplesInRefPopCountMap[individual] > self.TargetPerLineMap[individual]:
                        self.RemainingToDoPerLineMap[individual] = 0
                        self.Remainder1map[individual] = 0
                    else:
                        if self.TargetPerLineMap[individual] > self.LineAllInRefPopCountMap[individual]:
                            self.TargetPerLineMap[individual] = self.LineAllInRefPopCountMap[individual]
                        self.RemainingToDoPerLineMap[individual] = self.TargetPerLineMap[individual] - self.LinesWithSamplesInRefPopCountMap[individual]
                        if self.LinesWithSamplesInRefPopCountMap[individual] > 0:
                            self.IntermediateSum -= self.RemainingToDoPerLineMap[individual]
                        else:
                            self.IntermediateSum -= (self.RemainingToDoPerLineMap[individual] - 1)
                        if self.IntermediateSum < 0:
                            self.RemainingToDoPerLineMap[individual] += self.IntermediateSum
                            self.IntermediateSum = 0
                else:
                    if self.LinesWithSamplesInRefPopCountMap[individual] == 0:
                        if individual in self.RemainingToDoPerLineMap:
                            if self.RemainingToDoPerLineMap[individual] == 0:
                                self.RemainingToDoPerLineMap[individual] = 1
                        else:
                            self.RemainingToDoPerLineMap[individual] = 1

        for dam in sorted(self.Remainder1map, key=self.Remainder1map.__getitem__, reverse=True):
            if self.IntermediateSum == 0:
                break
            if self.LineAllInRefPopCountMap[dam] > (self.LinesWithSamplesInRefPopCountMap[dam] + self.RemainingToDoPerLineMap[dam]):
                self.TargetPerLineMap[dam] += 1
                if self.LinesWithSamplesInRefPopCountMap[dam] > self.TargetPerLineMap[dam]:
                    continue
                self.RemainingToDoPerLineMap[dam] += 1
                self.IntermediateSum -= 1

        nonassigned = 0
        if self.IntermediateSum > 0:
            count = 0
            for i in self.IndividualsMissingFounderList:
                if i in self.ReferencePopulationList:
                    count += 1
            if self.IntermediateSum > count:
                if mode == "gui":
                    messagebox.showwarning("Warning", "Warning, not enough individuals for such an ambitious sampling plan!")
                elif mode == "cl":
                    sys.stderr.write("Warning, not enough individuals for such an ambitious sampling plan!\n")
            if count > self.IntermediateSum:
                nonassigned = self.IntermediateSum
            else:
                nonassigned = count
            for individual in self.MixedList:
                if (self.RemainingToDoPerLineMap[individual] + self.LinesWithSamplesInRefPopCountMap[individual]) > self.LineAllInRefPopCountMap[individual]:
                    self.RemainingToDoPerLineMap[individual] = self.LineAllInRefPopCountMap[individual] - self.LinesWithSamplesInRefPopCountMap[individual]

        with open('OutputSampl_DetailedInfo.txt', 'w') as output1:
            output1.write(f"reference population size = {self.RefPopSize}\ndam/sire lines in ref.pop. = {len(self.MixedList)}\n")
            output1.write(f"IndividualsMissingFounderInRefPopCount = {self.IndividualsMissingFounderInRefPopCount}\n")
            output1.write(f"PreviouslySequencedInRefPop = {self.PreviouslySequencedInRefPop}\n{HowManyToSequence = }\n")
            output1.write(f"Number of lines without previous samples = {self.NoPreviousSamplesCount}\n")
            output1.write(f"TargetSum (that remains to be proportionally divided among dam/sire lines) = {self.TargetSum}\n\n")
            if self.Lineage == 'paternal':
                output1.write('FOUNDER SIRE:\n')
            elif self.Lineage == 'maternal':
                output1.write('FOUNDER DAM:\n')
            for individual in self.MixedList:
                if individual not in self.TargetPerLineMap:
                    self.TargetPerLineMap[individual] = 0
                if individual not in self.RemainingToDoPerLineMap:
                    self.RemainingToDoPerLineMap[individual] = 0
                output1.write(f"{individual} ->  in ref.pop.: {self.LineAllInRefPopCountMap[individual]}   targ: {self.TargetPerLineMap[individual]}   prev: {self.LinesWithSamplesInRefPopCountMap[individual]}  todo: {self.RemainingToDoPerLineMap[individual]}\n")

            if nonassigned > 0:
                output1.write(f"\nThere are {self.IntermediateSum} remaining non-assigned planned samplings\n")
                if self.Lineage == "maternal":
                    output1.write(f"The list of {count} individuals in reference population that can be sampled but that do not have the founder dam:\n")
                if self.Lineage == "paternal":
                    output1.write(f"The list of {count} individuals in reference population that can be sampled but that do not have the founder sire:\n")
                for individual in self.IndividualsMissingFounderList:
                    if individual in self.ReferencePopulationList:
                        output1.write(f"{individual}\n")

        pass

    def calcAvailNumOfSamplPerLine(self,HowManyToSequence,AVAILABLE_column):

        AvailablePerLineCountMap = {}
        RemainingDiffAvailableMap = {}
        AvailableAllCount = 0
        AvailableRealCount = 0
        AvailableDifferenceCount = 0
        self.RemainingToDoAvailablePerLineMap = {}
        if AVAILABLE_column != -1:
            for individual in self.MixedList:
                AvailablePerLineCountMap[individual] = 0
                self.RemainingToDoAvailablePerLineMap[individual] = self.RemainingToDoPerLineMap[individual]
            for individual in self.ReferencePopulationList:
                if self.AvailableMap[individual] == '1':
                    AvailableAllCount += 1
                    if self.IndividualToFounderMap[individual] != '0':
                        if individual not in self.HaplotypedList:
                            AvailableRealCount += 1
                            AvailablePerLineCountMap[self.IndividualToFounderMap[individual]] += 1
            for individual in self.MixedList:
                RemainingDiffAvailableMap[individual] = None
                if AvailablePerLineCountMap[individual] < self.RemainingToDoPerLineMap[individual]:
                    RemainingDiffAvailableMap[individual] = self.RemainingToDoPerLineMap[individual] - AvailablePerLineCountMap[individual]
                    AvailableDifferenceCount += RemainingDiffAvailableMap[individual]
                    self.RemainingToDoAvailablePerLineMap[individual] = AvailablePerLineCountMap[individual]
        else:
            for individual in self.MixedList:
                self.RemainingToDoAvailablePerLineMap[individual] = self.RemainingToDoPerLineMap[individual]

        if AvailableDifferenceCount > 0:
            with open('OutputSampl_AvailabilityRestrictions.txt', 'w') as output2:
                if self.Lineage == "maternal":
                    line = "dam"
                elif self.Lineage == "paternal":
                    line = "sire"
                output2.write(f"reference population size = {self.RefPopSize}    {line} lines in ref. pop. = {len(self.MixedList)}    {self.PreviouslySequencedInRefPop = }    {HowManyToSequence = }%s\n\n")
                output2.write(f"{AvailableAllCount = }    {AvailableRealCount = }    sum of DIFF (planned - available) = {AvailableDifferenceCount}    total of both available and planned: {HowManyToSequence-AvailableDifferenceCount}\n\n")
                for individual in sorted(self.MixedCountMap, key=self.MixedCountMap.__getitem__, reverse=True):
                    output2.write(f"{line}: {individual}   in ref.pop: {self.LineAllInRefPopCountMap[individual]}  targ: {self.TargetPerLineMap[individual]}  ")
                    output2.write(f"prev: {self.LinesWithSamplesInRefPopCountMap[individual]}  todo: {self.RemainingToDoPerLineMap[individual]}  ")
                    output2.write(f"avail: {AvailablePerLineCountMap[individual]}  DIFF: {RemainingDiffAvailableMap[individual]}  A_TODO: {self.RemainingToDoAvailablePerLineMap[individual]}")
                    if AvailablePerLineCountMap[individual] > self.RemainingToDoAvailablePerLineMap[individual]:
                        output2.write(f"  REMAINING AVAILABLE: {AvailablePerLineCountMap[individual] - self.RemainingToDoAvailablePerLineMap[individual]}\n")
                    else:
                        output2.write('\n')

        pass

    def selectIndividualsForSampling(self,SamplingMethod,mode):

        PercentageOfTotalPopMap = {}
        PercentageOfRefPopMap = {}
        if self.Lineage == "maternal":
            line = "dam"
        elif self.Lineage == "paternal":
            line = "sire"
        with open('OutputSampl_IndividualsForSampling.txt', 'w') as output3:
            output3.write(f"Selected individuals per {line}line:\n\n")
            chosen = fastSampler(self.ParentMap,
                                 self.ReferencePopulationList,
                                 self.IndividualToFounderMap,
                                 self.HaplotypedList,
                                 self.RemainingToDoAvailablePerLineMap,
                                 SamplingMethod,mode)
            for individual in chosen:
                PercentageOfTotalPopMap[individual] = (self.LineAllInPopCountMap[individual] / len(self.IDlist)) * 100
                PercentageOfRefPopMap[individual] = (self.LineAllInRefPopCountMap[individual] / len(self.ReferencePopulationList)) * 100
                output3.write(f"founder {line} {individual} : ")
                for choice in chosen[individual]:
                    output3.write(f"{choice}, ")
                output3.write('\n\n')

            output3.write(f"Percentages of {line}lines for sampling in total population:\n")
            for dam in sorted(PercentageOfTotalPopMap, key=PercentageOfTotalPopMap.__getitem__, reverse=True):
                output3.write(f"{line} {dam} : {PercentageOfTotalPopMap[dam]:.2f} %\n")
            output3.write('\n\n')
            output3.write(f"Percentages of {line}lines for sampling in reference population:\n")
            for dam in sorted(PercentageOfRefPopMap, key=PercentageOfRefPopMap.__getitem__, reverse=True):
                output3.write(f"{line} {dam} : {PercentageOfRefPopMap[dam]:.2f} %\n")
            output3.write('\n\n')
        pass

    def __call__(self,FirstRefYear,LastRefYear,FastMethod,HowManyToSequence,mode):
        AVAILABLE_column = self.setRefPop(FirstRefYear,LastRefYear,mode)
        if AVAILABLE_column == 1:
            return
        self.setFounderPop(mode)
        self.calcTargetNumSamplesPerLine(HowManyToSequence,mode)
        self.calcAvailNumOfSamplPerLine(HowManyToSequence,AVAILABLE_column)
        self.selectIndividualsForSampling(FastMethod,mode)
        return
