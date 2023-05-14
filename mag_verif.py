# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

from tkinter import messagebox
import sys
import numpy as np

def listLineToFounder(individual,ParentMap):
    LineToFounderList = list()
    # add all new members to the founder line
    ID_string = individual
    LineToFounderList.append(ID_string)
    while ParentMap[ID_string] != "0":
        ID_string = ParentMap[ID_string]
        LineToFounderList.append(ID_string)

    return LineToFounderList

# list offspring too
# Generations indicate how many generations after individual to add
# for all offspring, put Generations to -1
def listOffspringLine(individual,ParentMap,GenderMap,Generations):
    ## add all new members to the last child
    OffspringLineList = list()
    ID_strings = [individual]
    OffspringLineList.extend(ID_strings)
    i = 0
    while len(ID_strings):
        i += 1
        ID_strings = [key for key, value in ParentMap.items() if value in ID_strings]
        OffspringLineList.extend(ID_strings)
        if i == Generations:
            break

    #exclude all females for paternal lineage
    if GenderMap[individual] == '1':
        OffspringLineList[:] = [ID1 for ID1 in OffspringLineList if GenderMap[ID1] == '1']
    return OffspringLineList

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

def check_errors(IDlist,FatherMap,MotherMap,YobMap,GenderMap,mode):
# checking for errors in pedigree: cycles, gender consistency, and non-existent ancestors

    Female_gender = "2"
    Male_gender = "1"
    missing_entry = "0"
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
            try:
                if int(YobMap[individual]) < int(YobMap[FatherMap[individual]]):
                    YoungerFathersList.append(individual)
            except ValueError:
                # YobMap for individual or father is not an int
                pass
            except KeyError:
                # YobMap[father] does not exist
                pass

        if MotherMap[individual] != "0":
            try:
                if int(YobMap[individual]) < int(YobMap[MotherMap[individual]]):
                    YoungerMothersList.append(individual)
            except ValueError:
                # YobMap for individual or mother is not an int
                pass
            except KeyError:
                # YobMap[mother] does not exist
                pass

        if FatherMap[individual] != "0":
            if FatherMap[individual] not in IDlist:
                if FatherMap[individual] not in MissingFatherMap:
                    MissingFatherMap[FatherMap[individual]] = individual
                else:
                    MissingFatherMap.pop(FatherMap[individual])
                    IDlist.append(FatherMap[individual])
                    FatherMap[FatherMap[individual]] = missing_entry
                    MotherMap[FatherMap[individual]] = missing_entry
                    YobMap[FatherMap[individual]] = 'MISSING_YEAR'
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
                    YobMap[MotherMap[individual]] = 'MISSING_YEAR'
                    GenderMap[MotherMap[individual]] = Female_gender
                    AddedFemalesList.append(MotherMap[individual])

    for key in MissingFatherMap:
        FatherMap[MissingFatherMap[key]] = '0'
    for key in MissingMotherMap:
        MotherMap[MissingMotherMap[key]] = '0'

    if len(AddedMalesList) or len(AddedFemalesList) or len(MissingFatherMap) or len(MissingMotherMap):
        with open('autocorrection_log.txt', 'w') as corr_log:
            if len(AddedMalesList):
                corr_log.write('new records formed for male ancestors:\n')
                for male in AddedMalesList:
                    corr_log.write(f"     {male}\n")
            if len(AddedFemalesList):
                corr_log.write('new records formed for female ancestors:\n')
                for female in AddedFemalesList:
                    corr_log.write(f"     {female}\n")
            if len(MissingFatherMap):
                corr_log.write('removed unique and non-defined male ancestors:\n')
                for i in MissingFatherMap:
                    corr_log.write(f"     {i}\n")
            if len(MissingMotherMap):
                corr_log.write('removed unique and non-defined female ancestors:\n')
                for i in MissingMotherMap:
                    corr_log.write(f"     {i}\n")


    FatalError = False
    WarningError = False
    if len(ConflictingFatherList):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('Process is terminated. Please correct the following errors:\n')
            for individual in ConflictingFatherList:
                error_log.write(f"Individual {individual} is repeated as its father.\n")
        FatalError = True

    if len(ConflictingMotherList):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('Process is terminated. Please correct the following errors:\n')
            for individual in ConflictingMotherList:
                error_log.write(f"Individual {individual} is repeated as its mother.\n")
        FatalError = True

    if len(FemaleFathersList):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('Process is terminated. Please correct the following errors:\n')
            for individual in ConflictingMotherList:
                error_log.write('Process is terminated. Please correct the following error:\n')
                error_log.write(f"Individual {individual} has a father \n{FatherMap[individual]} that is defined as female.\n")
        FatalError = True

    if len(MaleMothersList):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('Process is terminated. Please correct the following errors:\n')
            for individual in ConflictingMotherList:
                error_log.write('Process is terminated. Please correct the following error:\n')
                error_log.write(f"Individual {individual} has a mother \n {MotherMap[individual]} that is defined as male.\n")
        FatalError = True

    if len(YoungerFathersList):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('WARNING! Following inconsistencies detected:\n')
            for individual in YoungerFathersList:
                error_log.write(f"Individual {individual} has younger father {FatherMap[individual]}\n")
        WarningError = True

    if len(YoungerMothersList):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('WARNING! Following inconsistencies detected:\n')
            for individual in YoungerMothersList:
                error_log.write(f"Individual {individual} has younger mother {MotherMap[individual]}\n")
        WarningError = True

    FatherLineConflicts = checkConflicts(FatherMap,GenderMap,Male_gender)
    MotherLineConflicts = checkConflicts(MotherMap,GenderMap,Female_gender)

    if len(FatherLineConflicts):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('Please correct the following errors:\n')
            for individual in FatherLineConflicts:
                error_log.write(f"Male individual {individual} has a repeated ancestor {FatherLineConflicts[individual]} in its lineage:\n")
                error_log.write("\n\n")
        FatalError = True

    if len(MotherLineConflicts):
        with open('ERROR_ALERT.TXT', 'w') as error_log:
            error_log.write('Please correct the following errors:\n')
            for individual in MotherLineConflicts:
                error_log.write(f"Female individual {individual} has a repeated ancestor {MotherLineConflicts[individual]} in its lineage:\n")
                error_log.write("\n\n")
        FatalError = True

    if FatalError:
        if mode == "gui":
            messagebox.showerror("Fatal error", "Fatal errors in pedigree, please see ERROR_ALERT.TXT.")
        elif mode == "cl":
            sys.stderr.write("Fatal errors in pedigree, please see ERROR_ALERT.TXT.\n")
            sys.exit()
    if WarningError:
        if mode == "gui":
            messagebox.showwarning("Warning", "Warning: YOB incosistency, please see ERROR_ALERT.TXT.")
        elif mode == "cl":
            sys.stderr.write("Warning: YOB incosistency, please see ERROR_ALERT.TXT.\n")

    return IDlist, FatherMap, MotherMap, YobMap, GenderMap, FatalError

def fillLowestList(InputList,ParentMap):
    ResultList = InputList
    for i in InputList:
        current = i
        while True:
            if ParentMap[current] != "0":
                current = ParentMap[current]
                if current in ResultList:
                    ResultList.remove(current)
            else:
                break
    return ResultList

def write_seq_diff(ConflictUnit,ConflictingUnits,HaplotypeMap,mode):
    SequenceLength = len(HaplotypeMap[ConflictUnit])
    IDLength = len(ConflictUnit)
    for unit in ConflictingUnits:
        if len(HaplotypeMap[ConflictUnit]) != len(HaplotypeMap[unit]):
            message = f"Sequence length of {ConflictUnit} is different from the sequence length of {unit}."
            if mode == "gui":
                messagebox.showerror("ERROR", message)
                return
            elif mode == "cl":
                sys.stderr.write(f"ERROR: {message}\n")
                return
    with open('OutputVerif_Mutations.txt', 'a') as out_mutations:
        out_mutations.write(f"position\t{ConflictUnit}")
        for unit in ConflictingUnits:
            out_mutations.write(f";\t{unit}")
        out_mutations.write("\n")
        for i in range(SequenceLength):
            out_mutations.write(f"{i+1:>8}:\t{HaplotypeMap[ConflictUnit][i]:>{IDLength}}")
            for unit in ConflictingUnits:
                if HaplotypeMap[ConflictUnit][i] != HaplotypeMap[unit][i]:
                    out_mutations.write(f";\t{HaplotypeMap[unit][i]:>{IDLength}}")
                else:
                    empty = " "
                    out_mutations.write(f";\t{empty:>{IDLength}}")
            out_mutations.write("\n")
        out_mutations.write("\n")
    return

def check_haplotype_conflicts(IDlist,ParentMap,HaplotypeMap,HaplotypedList,HaplotypeNamesList,VerifObject,mode):

    MatchCount = 0
    MismatchCount = 0
    InformativeList = set()
    UnitsInConflictsList = set()
    ConflictCountMap = {}

    for individual in HaplotypedList:
        ConflictCountMap[individual] = 0
    ConflictMatrix = np.zeros((len(HaplotypedList),len(HaplotypedList)), dtype=int)

    for i in range(len(HaplotypedList)):
        ParentLine1 = []
        ID_string = HaplotypedList[i]
        ParentLine1.append(ID_string)
        while ParentMap[ID_string] != "0":
            ID_string = ParentMap[ID_string]
            ParentLine1.append(ID_string)

        for j in range(i + 1, len(HaplotypedList)):
            ParentLine2 = []
            ID_string = HaplotypedList[j]
            ParentLine2.append(ID_string)
            while ParentMap[ID_string] != "0":
                ID_string = ParentMap[ID_string]
                ParentLine2.append(ID_string)

            # find common ancestor
            for k in range(len(ParentLine1)):
                CommonAncestorFlag = False
                for l in range(len(ParentLine2)):
                    if ParentLine2[l] == ParentLine1[k]:
                        for m in range(k):
                            InformativeList.add(ParentLine1[m])
                        for m in range(l):
                            InformativeList.add(ParentLine2[m])
                        if HaplotypeMap[ParentLine1[0]] == HaplotypeMap[ParentLine2[0]]:
                            MatchCount += 1
                        else:
                            ConflictCountMap[ParentLine1[0]] += 1
                            ConflictCountMap[ParentLine2[0]] += 1
                            ConflictMatrix[i,j] = 1
                            ConflictMatrix[j,i] = 1
                            MismatchCount += 1
                            UnitsInConflictsList.add(ParentLine1[0])
                            UnitsInConflictsList.add(ParentLine2[0])
                        CommonAncestorFlag = True
                        break
                if CommonAncestorFlag:
                    break
    
    with open('OutputVerif_ConflictingIndividuals.txt', 'w') as out_conflicts, \
         open('OutputVerif_DetailedConflictingIndividuals.txt', 'w') as out_conflicts_detail:
        if len(UnitsInConflictsList) > 0:
            out_conflicts.write('ID of the conflicting unit; Haplotype; No. of conflicts\n\n')
            out_conflicts_detail.write('ID of the conflicting unit; Haplotype; Conflicting members\n\n')
            if VerifObject == "snpseq":
                with open('OutputVerif_Mutations.txt', 'w') as out_mutations:
                    out_mutations.write('Point mutations found in SNP sequences\n\n')
        else:
            out_conflicts.write('There are no conflicts in the pedigree.\n')
            out_conflicts_detail.write('There are no conflicts in the pedigree.\n')
            if VerifObject == "snpseq":
                with open('OutputVerif_Mutations.txt', 'w') as out_mutations:
                    out_mutations.write('No point mutations found in SNP sequences\n\n')

        ConflictCountDecMap = {}
        for individual in HaplotypedList:
            ConflictCountDecMap[individual] = ConflictCountMap[individual]
        ConflictingUnitsList = []
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
                f"{maxConflictUnit}; {HaplotypeMap[maxConflictUnit]}; {ConflictCountMap[maxConflictUnit]}\n")
            out_conflicts_detail.write(
                f"{maxConflictUnit}; {HaplotypeMap[maxConflictUnit]}; ")
            ConflictingUnits = []
            for idx in range(len(HaplotypedList)):
                if ConflictMatrix[maxConflictUnitIndex,idx] == 1:
                    out_conflicts_detail.write(f"{HaplotypedList[idx]} ")
                    ConflictingUnits.append(HaplotypedList[idx])
                    ConflictCountDecMap[HaplotypedList[idx]] -= 1
            out_conflicts_detail.write('\n')
            if VerifObject == "snpseq":
                write_seq_diff(maxConflictUnit,ConflictingUnits,HaplotypeMap,mode)

    with open('OutputVerif_MisplacedBranches.txt', 'w') as out_misplaced:
        out_misplaced.write('Misplaced branches in the pedigree:\n\n')
        LowestConflictingUnitsList = fillLowestList(ConflictingUnitsList, ParentMap)
        MisplacedUnitsCount = 0
        MisplacedBranchesCount = 0
        for individual in LowestConflictingUnitsList:
            MisplacedBranchHaplotypedList = []
            current = individual
            MisplacedBranchHaplotypedList.append(current)
            while ParentMap[current] != "0":
                current = ParentMap[current]
                if current in HaplotypedList:
                    if HaplotypeMap[current] == HaplotypeMap[individual]:
                        MisplacedBranchHaplotypedList.append(current)
                        MisplacedUnitsCount += 1
                    else:
                        break
            if len(MisplacedBranchHaplotypedList) > 1:
                MisplacedBranchesCount += 1
                for i in MisplacedBranchHaplotypedList:
                    out_misplaced.write(f"{i}   haplotype = {HaplotypeMap[i]}\n")
                out_misplaced.write('\n\n')
        MisplacedUnitsCount += MisplacedBranchesCount
        if MisplacedUnitsCount == 0:
            out_misplaced.write('There are no misplaced branches in the pedigree.\n')

    with open('OutputVerif_Summary.txt', 'w') as out_summary:
        out_summary.write(f'No. of individuals in pedigree{len(IDlist): >25}\n')
        out_summary.write(f'No. of distinct haplotypes{len(HaplotypeNamesList): >29}\n')
        out_summary.write(f'No. of individuals with sequenced haplotype{len(HaplotypedList): >12}\n')
        out_summary.write(f'No. of pairwise mismatches{MismatchCount: >29}\n')
        pruned = len(ConflictingUnitsList) - MisplacedUnitsCount + MisplacedBranchesCount

        if len(ConflictingUnitsList) > 0:
            out_summary.write(
                f'Maximal no. of conflicts / individual{ConflictCountMap[ConflictingUnitsList[0]]: >12}/{ConflictingUnitsList[0]}\n')
        else:
            out_summary.write('Maximal no. of conflicts / individual' + ' ' * 12 + '0 / NA\n')

        out_summary.write(f'No. of informative individuals{len(InformativeList): >25}\n')
        out_summary.write(f'No. of individuals participating in conflicts{len(UnitsInConflictsList): >10}\n')
        out_summary.write(f'No. of conflicting or misplaced individuals{len(ConflictingUnitsList): >12}\n')
        out_summary.write(f'Pruned no. of conflicting individuals{pruned: >18}\n')
        out_summary.write('-----------------------------------------------------------\n')
        if len(ConflictingUnitsList) > 0:
            HC = len(ConflictingUnitsList) / len(HaplotypedList) * 100
            IC = len(ConflictingUnitsList) / len(InformativeList) * 100
            MISPLACED = (MisplacedUnitsCount - MisplacedBranchesCount) / len(ConflictingUnitsList) * 100
            out_summary.write(f'HC index = {HC}\nIC index = {IC}\nMISPLACED index = {MISPLACED}\n')
        else:
            out_summary.write('There are no conflicts in the pedigree.\n')

    return

def impute_haplotype(IDlist,Lineage,GenderMap,ParentMap,HaplotypeMap,HaplotypedList,HaplotypeNamesList,Reliability,mode):

    if not len(HaplotypedList):
        if mode == "gui":
            messagebox.showwarning("Message", "No haplotyped individuals in the pedigree.")
        elif mode == "cl":
            sys.stderr.write("No haplotyped individuals in the pedigree.\n")
        return

    if Lineage == 'maternal':
        Gender = '2'
    elif Lineage == 'paternal':
        Gender = '1'

    #Initialize lists and maps
    ImputedHaplotypeMap = {}
    ImputedHaplotypedList = set()
    FounderToMembersMap = {}
    UnVerifiedList = set()
    VerifiedList = set()
    
    # find all founders
    FounderList = [key for key, value in ParentMap.items() if (value == '0') and (GenderMap[key] == Gender)]
    # list all children (non-parents)
    #ChildrenList = [key for key, value in ParentMap.items() if (key not in ParentMap.values()) and (ParentMap[key] != '0')]
    

    # Initialize a map to save all members of the founder line
    for ID in FounderList:
        Offspring = listOffspringLine(ID,ParentMap,GenderMap,-1)
        FounderToMembersMap[ID] = Offspring

    # loop over all founders
    for founder in FounderList:

        # find all haplotyped members of the founder line
        HaplotypedMembersList = list(set(HaplotypedList).intersection(set(FounderToMembersMap[founder])))

        # if there are no haplotyped members in the lineage, go to the next founder
        if not len(HaplotypedMembersList):
            continue

        # else if there is only 1 haplotyped member in the line, impute up to its founder only
        elif len(HaplotypedMembersList) == 1:
            if Reliability == 'low':
                ChildrenLine = listOffspringLine(HaplotypedMembersList[0],ParentMap,GenderMap,-1)
                for ID in ChildrenLine:
                    ImputedHaplotypeMap[ID] = HaplotypeMap[HaplotypedMembersList[0]]
                    ImputedHaplotypedList.add(ID)
                UnVerifiedList.update(ChildrenLine[1:])
            elif Reliability == 'high':
                ImputedHaplotypeMap[HaplotypedMembersList[0]] = HaplotypeMap[HaplotypedMembersList[0]]
                ImputedHaplotypedList.add(HaplotypedMembersList[0])
            VerifiedList.add(HaplotypedMembersList[0])
            continue


        # else subset the HaplotypeMap only for founder line members
        LineHaplotypeMap = {k: HaplotypeMap[k] for k in HaplotypeMap.keys() & HaplotypedMembersList}
        
        # if there is only 1 unique haplotype in the lineage, there is no conflict
        if len(set(LineHaplotypeMap.values())) == 1:
            Haplotype = list(LineHaplotypeMap.values())[0]
            # Make a single loop over haplotyped members and 
            # keep the first haplotyped member as a reference
            ID1 = HaplotypedMembersList[0]
            AncestorsID1 = listLineToFounder(ID1,ParentMap)
            OffspringID1 = listOffspringLine(ID1,ParentMap,GenderMap,-1)
            # loop over offspring and impute to unverified if there
            # are no haplotyped members among them
            if not (set(OffspringID1[1:]) & set(HaplotypedMembersList)):
                for ID in OffspringID1[1:]:
                    if Reliability == 'low':
                        UnVerifiedList.add(ID)
                        ImputedHaplotypeMap[ID] = Haplotype
                        ImputedHaplotypedList.add(ID)
            for ID2 in HaplotypedMembersList[1:]:
                AncestorsID2 = listLineToFounder(ID2,ParentMap)
                OffspringID2 = listOffspringLine(ID2,ParentMap,GenderMap,-1)
                node = [member for member in AncestorsID1 if member in AncestorsID2][0]
                for ID in AncestorsID1:
                    VerifiedList.add(ID)
                    ImputedHaplotypeMap[ID] = Haplotype
                    ImputedHaplotypedList.add(ID)
                    if Reliability == 'low':
                        # add also other offspring of connecting nodes to unverified
                        # if there is no haplotyped individuals among them
                        Offspring1 = listOffspringLine(ID,ParentMap,GenderMap,1)
                        for IDx in Offspring1[1:]:
                            # we impute offspring of ID1 separately
                            if (IDx not in VerifiedList) and (IDx not in UnVerifiedList):
                                if IDx != ID1:
                                    Offspring2 = listOffspringLine(IDx,ParentMap,GenderMap,-1)
                                    if not (set(Offspring2[1:]) & set(HaplotypedMembersList)):
                                        for IDy in Offspring2:
                                            UnVerifiedList.add(IDy)
                                            ImputedHaplotypeMap[IDy] = Haplotype
                                            ImputedHaplotypedList.add(IDy)

                    if ID == node:
                        break
                for ID in AncestorsID2:
                    VerifiedList.add(ID)
                    ImputedHaplotypeMap[ID] = Haplotype
                    ImputedHaplotypedList.add(ID)
                    if Reliability == 'low':
                        # add also other offspring of connecting nodes to unverified
                        # if there is no haplotyped individuals among them
                        Offspring1 = listOffspringLine(ID,ParentMap,GenderMap,1)
                        for IDx in Offspring1[1:]:
                            # we impute offspring of ID2 separately
                            if (IDx not in VerifiedList) and (IDx not in UnVerifiedList):
                                if IDx != ID2:
                                    Offspring2 = listOffspringLine(IDx,ParentMap,GenderMap,-1)
                                    if not (set(Offspring2[1:]) & set(HaplotypedMembersList)):
                                        for IDy in Offspring2:
                                            UnVerifiedList.add(IDy)
                                            ImputedHaplotypeMap[IDy] = Haplotype

                    if ID == node:
                        break
                if Reliability == 'low':
                    if not (set(OffspringID2[1:]) & set(HaplotypedMembersList)):
                        for ID in OffspringID2[1:]:
                            if ID in VerifiedList:
                                break
                            UnVerifiedList.add(ID)
                            ImputedHaplotypeMap[ID] = Haplotype
                            ImputedHaplotypedList.add(ID)
            continue

        # if there are more unique haplotypes in the lineage, there is a conflict somewhere
        if len(set(LineHaplotypeMap.values())) > 1:
            # Loop over haplotyped member pairs and 
            # impute if there is no conflict in the pair
            for idx, ID1 in enumerate(HaplotypedMembersList):
                Haplotype1 = LineHaplotypeMap[ID1]
                AncestorsID1 = listLineToFounder(ID1,ParentMap)
                OffspringID1 = listOffspringLine(ID1,ParentMap,GenderMap,-1)
                for ID2 in HaplotypedMembersList[idx+1:]:
                    Haplotype2 = LineHaplotypeMap[ID2]
                    # if the two haplotypes accord, impute
                    if Haplotype1 == Haplotype2:
                        AncestorsID2 = listLineToFounder(ID2,ParentMap)
                        OffspringID2 = listOffspringLine(ID2,ParentMap,GenderMap,-1)
                        node = [member for member in AncestorsID1 if member in AncestorsID2][0]
                        for ID in AncestorsID1:
                            VerifiedList.add(ID)
                            ImputedHaplotypeMap[ID] = Haplotype1
                            ImputedHaplotypedList.add(ID)
                            if Reliability == 'low':
                                # add also other offspring of connecting nodes to unverified
                                # if there is no haplotyped individuals among them
                                Offspring1 = listOffspringLine(ID,ParentMap,GenderMap,1)
                                for IDx in Offspring1[1:]:
                                    # we impute offspring of ID1 separately
                                    if (IDx not in VerifiedList) and (IDx not in UnVerifiedList):
                                        if IDx != ID1:
                                            Offspring2 = listOffspringLine(IDx,ParentMap,GenderMap,-1)
                                            if not (set(Offspring2[1:]) & set(HaplotypedMembersList)):
                                                for IDy in Offspring2:
                                                    UnVerifiedList.add(IDy)
                                                    ImputedHaplotypeMap[IDy] = Haplotype1
                            if ID == node:
                                break
                        for ID in AncestorsID2:
                            VerifiedList.add(ID)
                            ImputedHaplotypeMap[ID] = Haplotype2
                            ImputedHaplotypedList.add(ID)
                            if Reliability == 'low':
                                # add also other offspring of connecting nodes to unverified
                                # if there is no haplotyped individuals among them
                                Offspring1 = listOffspringLine(ID,ParentMap,GenderMap,1)
                                for IDx in Offspring1[1:]:
                                    # we impute offspring of ID2 separately
                                    if (IDx not in VerifiedList) and (IDx not in UnVerifiedList):
                                        if IDx != ID2:
                                            Offspring2 = listOffspringLine(IDx,ParentMap,GenderMap,-1)
                                            if not (set(Offspring2[1:]) & set(HaplotypedMembersList)):
                                                for IDy in Offspring2:
                                                    UnVerifiedList.add(IDy)
                                                    ImputedHaplotypeMap[IDy] = Haplotype2
                            if ID == node:
                                break
                        if Reliability == 'low':
                            if not (set(OffspringID1[1:]) & set(HaplotypedMembersList)):
                                for ID in OffspringID1[1:]:
                                    UnVerifiedList.add(ID)
                                    ImputedHaplotypeMap[ID] = Haplotype1
                                    ImputedHaplotypedList.add(ID)
                            if not (set(OffspringID2[1:]) & set(HaplotypedMembersList)):
                                for ID in OffspringID2[1:]:
                                    if ID in VerifiedList:
                                        break
                                    UnVerifiedList.add(ID)
                                    ImputedHaplotypeMap[ID] = Haplotype2
                                    ImputedHaplotypedList.add(ID)
                    # if there is a conflict in haplotypes, impute only those two
                    else:
                        VerifiedList.add(ID1)
                        VerifiedList.add(ID2)
                        ImputedHaplotypeMap[ID1] = Haplotype1
                        ImputedHaplotypeMap[ID2] = Haplotype2
                        ImputedHaplotypedList.add(ID1)
                        ImputedHaplotypedList.add(ID2)

            continue

    with open('OutputVerif_Imputed.txt', 'w') as imput_out:
        imput_out.write("Imputed haplotypes:\n")
        for ID in ImputedHaplotypeMap:
            imput_out.write(f"{ID}:\t{ImputedHaplotypeMap[ID]}\n")
    if Reliability == 'low':
        with open('OutputVerif_ImputedVerified.txt', 'w') as verif_out:
            verif_out.write("Imputed verified individuals:\n")
            for ID in VerifiedList:
                verif_out.write(f"{ID}\n")
        with open('OutputVerif_ImputedUnVerified.txt', 'w') as unverif_out:
            unverif_out.write("Imputed unverified individuals:\n")
            for ID in UnVerifiedList:
                unverif_out.write(f"{ID}\n")

    return ImputedHaplotypeMap, VerifiedList
