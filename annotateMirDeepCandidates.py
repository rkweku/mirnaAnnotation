import csv
import datetime
import os
import sys

from Bio.Blast.Applications import NcbiblastnCommandline

sys.path.append(os.path.join(os.path.abspath(os.path.dirname(
    sys.argv[0])),'..','miRador/miRador'))

import setupMiRBase

libFolder = "miRDP21.1.2/AT_pub2_sRNA_noBypass_forSparta"
organism = "ath"
version = "CURRENT"
subjectSequencesFilename = "/home/rkweku/miRador/miRBase/miRBasePlantMirnas.fa"
dbFilename = "/home/rkweku/miRador/miRBase/miRBasePlantMirnas.db"
outputFilename = "mirDeep_AT_pub2_sRNA_noBypass_forSparta.csv"
fastaFilename = "mirDeep_AT_pub2_sRNA_noBypass_forSparta.fa"
annotatedFilename = "mirDeep_AT_pub2_sRNA_noBypass_forSparta-annotated.csv"
queryMirnasFilename = fastaFilename
preAnnotationFilename = outputFilename

def readMirdeepOutput(filename, precursorsDict):
    """Read the miRDeeP-P2 output

    Args:
        filename: Name of the file to be read
        precursorsDict: Continually passed dictionary to track already known
            miRNAs and their precursors

    Returns:
        Variable wholeFile containing the entire file

    """

    # Create empty array to hold full csv file
    wholeFile = []

    # loop through file using csv.reader and store in wholeFile
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter = "\t")
        for row in reader:
            wholeFile.append(row)
            chr = row[0]
            strand = row[1]
            abundance = int(row[2].split("_x")[1])
            position = int(row[4].split(".")[0])
            mirSeq = row[6].replace("T", "U")

            coordinates = (chr, strand, position)

            if(coordinates not in precursorsDict):
                precursorsDict[coordinates] = {}

            if(mirSeq not in precursorsDict[coordinates]):
                precursorsDict[coordinates][mirSeq] = {}

            precursorsDict[coordinates][mirSeq][filename] = abundance

    return(precursorsDict)

def mergeMirdeepFiles(precursorsDict, libFilenamesList, fastaFilename):
    """Merge the mirDeeP-P2 output files 

    Args:
        precursorsDictByLib: Dictionary of the precursors and their miRNA
            for each library
        libFilenamesList: A list of the library file

    Returns:
        Table of the unique miRNAs with abundances of each library where
        that candidate exists

    """

    counter = 1

    f = open(outputFilename, "w")
    g = open(fastaFilename, "w")

    # Write the header
    f.write("miR Name,Chr,Strand,miR Position,miR Sequence")

    for libFilename in libFilenamesList:
        libFilenameStripped = os.path.basename(libFilename)
        libFilenameStripped = libFilenameStripped.split("-15-0-10_filter_P_prediction.tsv")[0]
        f.write(",%s" % libFilenameStripped)

    f.write("\n")

    for coordinates in precursorsDict:
        chr = coordinates[0]
        strand = coordinates[1]
        if(strand == "+"):
            strand = "w"
        elif(strand == "-"):
            strand = "c"
        position = coordinates[2]

        for mirSeq, libList in precursorsDict[coordinates].items():
            mirSeq = mirSeq.replace("T", "U")
            f.write("candidate-%s,%s,%s,%s,%s" % (counter, chr, strand,
                position, mirSeq))
            g.write(">candidate-%s\n%s\n" % (counter, mirSeq))
            for library in libFilenamesList:
                if(library in libList):
                    abun = precursorsDict[coordinates][mirSeq][library]
                    f.write(",%s" % abun)
                else:
                    f.write(",0")

                if(library == libFilenamesList[-1]):
                    f.write("\n")

            counter += 1

    f.close()
    g.close()
 
def getFastaDate(filename):
    """Get the date of the input filename creation date

    Args:
        filename: The name of the subject sequence filaneme

    Returns:
        The date of the blast DB update

    """

    # Get the time in the struct_time format
    updateTime = datetime.datetime.fromtimestamp(os.path.getmtime
        (filename))

    return(updateTime)

def getLocalDBDate(dbName):
    """Get the date of the blast DB creation date

    Args:
        dbName: The name of the blast DB

    Returns:
        The date of the blast DB update

    """

    # Only check the creation date if the file even exists
    if(os.path.isfile('%s.nhr' % dbName)):
        # Get the time in the struct_time format
        localDatabaseTime = datetime.datetime.fromtimestamp(os.path.getmtime
            ('%s.nhr' % dbName))

        return(localDatabaseTime)

    # File doesn't exist so return year 1 for comparison sake.
    # Note: datetime.MINYEAR threw error hence the hardcoded date
    else:
        print("There is no file with name %s. Creating BLAST database."\
            % dbName)
        return(datetime.datetime(1,1,1))

def checkNeedUpdateByDate(subjectSequencesFilename, dbName):
    """Check to see if the BLAST database needs to be updated

    Args:
        subjectSequencesFilename: The name of the subject sequences file

    Returns:
        Bool value indicating if the BLAST DB needs to be updated

    """

    # Get the date that the FASTA file was last updated
    fastaFileDate = getFastaDate(subjectSequencesFilename)

    # Get the date that the database was updated last
    localDatabaseDate = getLocalDBDate(dbName)

    # If the expression database is newer than the BLAST database,
    # an update is needed
    if(fastaFileDate > localDatabaseDate):
        print("An update is needed to the BLAST database and will be created "\
        "now. This will take time so please be patient.")
        return(True)

    # If the expression database is older than the BLAST database,
    # no update is needed
    else:
        print("No update is needed to the BLAST database. Proceeding with "\
        "the BLAST.")
        return(False)


def createBlastDB(subjectSequencesFilename, dbFilename):
    """Create a local DB with the subject sequences for BLAST later

    Args:
        subjectSequencesFilename: The name of the file FASTA file holding 
            all the sequenes to be BLASTed against
        dbFilename: The name of the database file that will be written to

    """

    # Create a blast DB from the small RNAs
    os.system("makeblastdb -in %s -dbtype nucl -parse_seqids -out %s" %
        (subjectSequencesFilename, dbFilename))

def blastMirnas(subjectSequencesFilename, dbFilename,
        candidateSequencesFilename):
    """Create a subject database for the known miRNAs and BLAST
       all candidate miRNAs to those miRNAs to find the evolutionarily
       conserved miRNAS

    Args:
        flastFilename: The name of the FASTA file holding the subject
            sequences
        dbFilename: The name of the database file that will be written to
        candidateSequencesFilename: The name of the FASTA file holding the
            candidate sequences
        outputFilename: The name of the output file
    Returns:
        The BLAST results in plain text format

    """

    blastFilename = "blastResults.txt"

    # Run blastn-short, but set word size to 11 as 5 is too short IMO
    NcbiblastnCommandline(query=candidateSequencesFilename,
        db=dbFilename, task='blastn', outfmt=6,
        num_threads=8, strand='plus', out=blastFilename)()[0]

    return(blastFilename)

def addSequencesToOutput(querySequencesFilename, subjectSequencesFilename):
    """Add the full sequences to the XML file for each alignment

    Args:
        querySequencesFilename: Name of the file holding the query
            mirnas
        subjectSequencesFilename: Name of the file holding all subject 
            sequences and IDs

    """

    blastFilename = "blastResults.txt"

    # Parse query and subject sequences 
    querySequences = getSequencesFromFasta(querySequencesFilename)
    subjectSequences = getSequencesFromFasta(subjectSequencesFilename)

    # Read entire output file into to be traversed
    blastResults = readFile(blastFilename, '\t')

    # Reopen outputFile to add the sequences to it
    f_out = open(blastFilename, 'w')

    # Loop through the output file to find which row to add the sequence to
    for row in blastResults:
        # Get the query and subject IDs from the row
        queryID = row[0]
        subjectID = row[1]

        # Using the query and subject IDs, get the query and subject
        # sequences from the respective dictionaries
        querySeq = querySequences[queryID].replace("T", "U")
        subjectSeq = subjectSequences[subjectID].replace("T", "U")

        # Loop through the row to write the new row to the output file
        for entry in row:
            # As long as the entry is not the last entry, print the entry
            # with a tab after
            if entry != row[-1]:
                f_out.write('%s\t' % entry)

                # If the entry is either the queryID or subjectID, print the
                # corresponding sequence after it
                if(entry == queryID):
                    f_out.write('%s\t' % querySeq)

                if(entry == subjectID):
                    f_out.write('%s\t' % subjectSeq)

            # If the entry is the last element in the row, just write it to the
            # file without a tab
            else:
                f_out.write('%s\n' % entry)

    f_out.close()

def getSequencesFromFasta(filename):
    """Read a FASTA file into memory as a dictionary

    Args:
        filename: Name of the file to be read
    Returns:
        Variable wholeFile containing the entire file

    """

    # Create empty dictionary to hold full FASTA file
    queryMirnas = {}

    # loop through file using csv.reader and store in wholeFile
    f = open(filename, 'r')
    wholeFile = f.readlines()

    # Loop through the miRNAs and save the miRNA tags only
    for i in range(0,len(wholeFile),2):
        # Strip > and newline character from the ID
        seqID = wholeFile[i][1:-1].split(' ')[0].rstrip()
        # Strip the newline character from the sequence
        sequence = wholeFile[i+1]
        # Store the seqID in the dictionary with sequence as definition
        # Strip the newline character from the sequence
        queryMirnas[seqID] = sequence[:-1]

    return(queryMirnas)

def readFile(filename, separator):
    """Read a file into memory

    Args:
        filename: Name of the file to be read
        separator: Delimiter of the file

    Returns:
        Variable wholeFile containing the entire file

    """

    # Create empty array to hold full csv file
    wholeFile = []

    # loop through file using csv.reader and store in wholeFile
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter = separator)
        for row in reader:
            wholeFile.append(row)

    return(wholeFile)

def readBlastResults(filename):
    """Read in the results of the BLAST file after it has had
    its sequences added back to it

    Args:
        filename: The name of the blast file to be read in

    Returns:
        A dictionary of the BLAST results with the candidate miRNA
        name as the key and a dictionary of results as the value

    """

    blastDict = {}

    # Loop through the blast results fil and add the data to a dictionary
    with open(filename) as f:
        for line in f:
            # Remove the new lin character and split the line on tabs
            blastResult = line.split('\n')[0].split('\t')

            # Store each element of the blast result as a separate
            # variable
            queryName = blastResult[0]
            subjectName = blastResult[2]
            # Separate the subject organism identifiers from the miRNA
            # family
            subjectOrganismIdentifier = subjectName.split('-')[0]

            # Some miRNAs in miRBase have miR-### for some reason. Check
            # if this is the case, and if so, process it differently
            if("miR-" in subjectName):
                mirnaFamilyName = "miR-%s" % (subjectName.split('-')[2])
            else:
                mirnaFamilyName = subjectName.split('-')[1]

            # We only want the miRNA family idenifier, so remove anything
            # else from the name before adding to blastDict 
            if(not mirnaFamilyName[-1].isdigit()):
                counter = -1
                while(not mirnaFamilyName[counter].isdigit()):
                    counter -= 1 
                    trimmedMirnaFamilyName = mirnaFamilyName[0:counter+1]

                mirnaFamilyName = trimmedMirnaFamilyName

            # If the query name is not the blast dictionary, initialize
            # it as an empty dictionary
            if(queryName not in blastDict):
                blastDict[queryName] = {}

            # If the mirna family doesn't exist yet in the blastDict
            # subdictionary, initialize the list for organisms before
            # adding to it
            if(mirnaFamilyName not in blastDict[queryName]):
                blastDict[queryName][mirnaFamilyName] = []

            blastDict[queryName][mirnaFamilyName].append(blastResult)

    return(blastDict)

def createSimilarityDict(blastDict, organism):
    """Create a similarity dictionary to assist in the annotation of the
    candidate miRNAs that either sho w the equivalence of the candidate
    sequence with a known sequence, or that a candidate miRNA is
    highly similar to a known family in a number of species

    Args:
        blastDict: A dictionary of the blast results with the candidate
            miRNA name as the key and a dictionary of similar known
            miRNAs as the value. This dictionary has the miRNA family
            as a key and then a list of blast results as the value
        organism: First letter of genus and first 2 
            letters of species.

    Returns:
        A dictionary that showing the equivalence or near equivalence of
        a candidate miRNA to known miRNAs

    """

    similarityDict = {}

    for candidateMirna, candidateBlastDict in blastDict.items():
        identicalFlag = 0
        for mirnaFamily, blastResultsList in candidateBlastDict.items():
            for blastResult in blastResultsList:
                candidateSeq = blastResult[1]
                subjectName = blastResult[2]
                subjectSeq = blastResult[3]
                subjectOrganismIdentifier = subjectName.split('-')[0]
                percIdentity = float(blastResult[4])
                alignLength = int(blastResult[5])
                numMismatch = int(blastResult[6])
                numGaps = int(blastResult[7])

                numMatches = alignLength - numMismatch - numGaps

                # Right now, we are going to classify anything that has
                # fewer than 5 total differences between the candidate
                # sequence and a miRNA known miRNA as a member of that family
                if(len(candidateSeq) - numMatches < 5):
                    # If candidateMirna or mirnaFamily haven't yet been
                    # initialized in similarityDict, initiailize them
                    if(not identicalFlag and candidateMirna not in similarityDict):
                        similarityDict[candidateMirna] = {}
                    if(not identicalFlag and mirnaFamily not in
                            similarityDict[candidateMirna]):
                        similarityDict[candidateMirna][mirnaFamily] = []

                    # If the candidate sequence and subject sequences are the
                    # same (and of course the same organism), then set the
                    # the value of the similarityDict just to the miRNA name
                    # and skip the rest of the blast results for this sequence
                    if(candidateSeq == subjectSeq and organism ==
                            subjectOrganismIdentifier):
                        # If the similarity dictionary for this candidate
                        # miRNA has not been reinitialized as a list, reset
                        # it entirely to be a list
                        if(isinstance(similarityDict[candidateMirna], dict)):
                            similarityDict[candidateMirna] = []

                        similarityDict[candidateMirna].append(subjectName)
                        # No need to continue with this candidate miRNA, so
                        # break from the loop to inspect the next candidate
                        identicalFlag = 1

                    # If the exact sequence hasn't been found, add the 
                    # organism identifier to the similarity dictionary for
                    # this miRNA family. Also, skip this case if a tag has
                    # been identified as being identical to this candidate
                    if(not identicalFlag and subjectOrganismIdentifier not in
                            similarityDict[candidateMirna][mirnaFamily]):
                        similarityDict[candidateMirna][mirnaFamily].append(
                            subjectOrganismIdentifier)

    return(similarityDict)

def annotateIdenticalCandidates(similarityDict, mirBaseDict, identicalList,
        line):
    """Helper function to annotate the candidate miRNAs that have identical
    sequences to ones that have already been identified.

    Args:
        mirBaseDict: Dictionary of miRBase miRNAs and their coordinates for
            this organism, if available (will be an empty dictionary if not)
        identicalList: The list of miRBase miRNAs with the same sequence
            as the candidate miRNA
        line: The full line from the pre-annotated file that will be
            modified to provide the new annotation

    Returns:
        Flag to indicate if a positional match was found for this candidate
        and the update line with the proper annotation

    """

    annotatedFlag = False
    mirName = line[0]
    chrName = line[1]
    strand = line[2]
    position = int(line[3])
    mirSeq = line[4]

    # Remove "chr" if it exists in the chromosome name
    if("chr" in chrName.lower()):
        chrName = chrName.lower().replace("chr", "")

    # If mirBaseDict is populated, that means that we have positional
    # information for this organism and thus can generate the most accurate
    # annotations for this organism
    if(mirBaseDict):
        # Loop through all identical miRNA sequences
        for identicalMirna in identicalList:
            # It turns out that there can be annotated miRNAs in miRBase
            # that do not exist in the gff file, so do a check to ensure
            # that the identical miRNA exists in mirBaseDict prior to
            # entering this loop
            if(identicalMirna not in mirBaseDict):
                continue

            # Loop through all coordinates that this specific miRNA
            # can be found
            for coordinates in mirBaseDict[identicalMirna]:
                mirBaseChr = coordinates[0]
                # Remove "chr" if it exists in the chromosome name
                if("chr" in mirBaseChr.lower()):
                    mirBaseChr = mirBaseChr.lower().replace("chr", "")

                mirBaseStrand = coordinates[1]

                # Need to convert strand from +/- to w/c
                if(mirBaseStrand == "+"):
                    mirBaseStrand = "w"
                    mirBasePosition = int(coordinates[2])
                elif(mirBaseStrand == "-"):
                    mirBaseStrand = "c"
                    mirBasePosition = int(coordinates[2])
                # If there strand is not + or -, something is wrong with this
                # miRBase entry. We will not exit the run, but we will report
                # the issue to the user and continue to the next tag
                else:
                    print("Unrecognized strand of miRBase entry. We will "\
                        "skip this entry, but please check with the miRBase "\
                        " file %s.gff3 and miRNA name %s" % (organism,
                        identicalMirna))
                    continue

                # If the coordinates of the candidate miRNA meet the
                # coordinates of the miRBase miRNA, then this candidate miRNA
                # is this miRBase miRNA and we will change the annotation to
                # represent that.
                # I'm having trouble with mirDeeP-P2's coordinate system. I'm
                # going to use a +- 2pb window to determine position
                # equivalence
                if(chrName == mirBaseChr and strand == mirBaseStrand and
                        position - mirBasePosition >= -2 and
                        position - mirBasePosition <= 2):
                    line[0] = identicalMirna
                    line.append("Known")
                    annotatedFlag = True

    # If we did not find an annotated miRNA at this same position, we will
    # annotate it in the final file as being identical to the following known
    # miRNAs, but not at the same position
    if(not annotatedFlag):
        toAdd = ""
        # Loop through all identical tags and copy a line for it so that all
        # known miRBase miRNAs matching this read can be found
        line.append("Identical to the following known miRNAs at different "\
            "positions:")
        for identicalMirna in similarityDict[mirName]:
            toAdd = "%s%s " % (toAdd, identicalMirna)
        line.append(toAdd)

    return(annotatedFlag, line)

def annotateMirDeepCandidates(similarityDict, organism, mirBaseDict,
        numLibs):
    """Using the similarityDict, annotate the data in the pre-annotated file
    with the proper miRNA name, with the already existing family(ies) that
    the candidate miRNA likely belongs to, or the conservation of an existing
    existing family that has not yet been identified in this organism but
    is present in another. If a candidate miRNA is novel, it MUST have been
    identified in more than one library to stay

    Args:
        similarityDict: Dictionary of the candidate miRNAs and their
            potential annotations if there is no equivalence, or just
            simply its known name if it already exists
        organism: First letter of genus and first two
            of the species
        mirBaseDict: Dictionary of miRBase miRNAs and their coordinates for
            this organism, if available (will be an empty dictionary if not)
        numLibs: The total number of libraries in the analysis

    """

    foundKnownDict = {}

    preAnnotatedFile = readFile(preAnnotationFilename, ',')

    f = open(annotatedFilename, 'w')
    fastaOut = open("%s.fa" % os.path.splitext(annotatedFilename)[0], "w")

    header = preAnnotatedFile[0]
    header.append("Classification Flag")

    # Get the index of the first library. Useful for
    # when we need to validate novel tags in more than 1 library,
    # but we won't use it for conserved families already found
    # in this species
    startIterIndex = header.index("miR Sequence") + 1
    endIndex = len(header) - 1

    for entry in header:
        if(entry == header[-1]):
            f.write("%s\n" % entry)
        else:
            f.write("%s," % entry)

    for line in preAnnotatedFile[1:]:
        mirName = line[0]
        libCount = 0
        annotatedFlag = False
        identicalFlag = False
        similarFlag = False

        # Check if the candidate miRNA name is in the similarityDict
        # to determine if it is completely novel or not
        if(mirName in similarityDict):
            # If the value in similarityDict of the candidate miRNA
            # is simply a list, that means that the candidate miRNA
            # is already present within this organism in mirBase, thus
            # replace the candidate name with the proper annotated name
            if(isinstance(similarityDict[mirName], list)):
                identicalFlag = True
                identicalList = similarityDict[mirName]
                copyList = []

                # Annotate the identical candidate miRNA. The tabs became
                # too deep for this relatively simple function
                annotatedFlag, line = annotateIdenticalCandidates(
                    similarityDict, mirBaseDict, identicalList,
                    line)

                # Write the line to a file once the candidate miRNA has
                # been renamed
                for i in range(len(line)):
                    if i == len(line) - 1:
                        f.write("%s\n" % line[i])
                    else:
                        f.write("%s," % line[i])

                # Write the sequence to the FASTA file
                fastaOut.write(">%s\n%s\n" % (line[0], line[4]))

            # If the candidate miRNA is similar to sequences in mirBase, but
            # not identical to anything in this organism, we will need to
            # write this information to the output file.
            else:
                # Loop through each miRNA family that the candidate sequence
                # had fewer than 5 differences to
                for mirFamily, organismList in similarityDict[mirName].items():
                    # Initialize libCount to 0 as we are looping through
                    # miRNA families and we don't want to conserve the
                    # count on the next family
                    libCount = 0

                    # If the organism being studied has a very similar
                    # sequence to one that is already known for this organism
                    # in miRBase, we will tag it as a member of this family
                    # Add the annotation to just after the last library to
                    # ensure this flag will override any conserved family flag
                    # on this same line
                    if(organism in organismList):
                        similarFlag = True
                        line.insert(endIndex, "New member of existing family")
                        line.insert(endIndex + 1, "%s-%s-like" % (organism,
                            mirFamily))

                    # If the sequence is only similar to miRNAs found in
                    # other organismList, don't tag as a member of any family.
                    # Rather, just write the miRNA family name
                    else:
                        for i in range(startIterIndex, endIndex):
                            # If the miRNA has been identified in this
                            # library, increment libCount
                            if(float(line[i])):
                                libCount += 1
                        if(libCount > 1 and libCount > float(numLibs * .1)):
                            similarFlag = True
                            line.append("Conserved family of the following:")
                            line.append(mirFamily)

                            toWrite = ""
                            # Add the list of organismList with this same
                            # miRNA family that met the similarity
                            # requirement
                            for similarOrganism in sorted(organismList):
                                if(similarOrganism == sorted(
                                        organismList[-1])):
                                    toWrite += "%s" % similarOrganism
                                else:
                                    toWrite += "%s " % similarOrganism

                            line.append(toWrite)

                # If the number of libraries this miRNA was predicted in is
                # greater than 1 as well as greater than 10% of the given
                # libraries, then we will confirm the replication requirement
                if(similarFlag):
                    for i in range(len(line)):
                        if i == len(line) - 1:
                            f.write("%s\n" % line[i])
                        else:
                            f.write("%s," % line[i]) 

                    # Write the sequence to the FASTA file
                    fastaOut.write(">%s\n%s\n" % (line[0], line[4]))

        # If the candidate miRNA had no similar sequence, it is completely
        # novel by our tests and thus requires validation in more than one
        # library. Thus, here we will literate through the results and
        # remove candidate miRNAs that were identified in just one library
        else:
            # Identify how many libraries this miRNA has been predicted
            # in
            for i in range(startIterIndex, endIndex):
                if(float(line[i])):
                    libCount += 1

            # If the number of libraries this miRNA was predicted in is
            # greater than 1 as well as greater than 10% of the given
            # libraries, then we will confirm the replication requirement
            if(libCount > 1 and libCount > float(numLibs * .1)):
                line.append("Novel")

                for i in range(len(line)):
                    if i == len(line) - 1:
                        f.write("%s\n" % line[i])
                    else:
                        f.write("%s," % line[i])

                # Write the sequence to the FASTA file
                fastaOut.write(">%s\n%s\n" % (line[0], line[4]))

    f.close()
    fastaOut.close()

def main():
    libFilenamesList = []
    wholeFileList = []
    precursorsDict = {}

    # Call setup miRBase to initialize the miRBase folder and required
    # files as well as populate a dictionary with miRBase miRNAs for the
    # organism being studied
    mirBaseDict = setupMiRBase.setupMiRBase(organism, version)

    for file in os.listdir(libFolder):
        if file.endswith(".tsv"):
            libFilenamesList.append("%s/%s" % (libFolder,
                os.path.join(file)))

    for libFilename in libFilenamesList:
        precursorsDict = readMirdeepOutput(libFilename, precursorsDict)

    mergeMirdeepFiles(precursorsDict, libFilenamesList, fastaFilename)

    # Check to see if the BLAST database needs to be updated
    updateFlag = checkNeedUpdateByDate(
        subjectSequencesFilename, dbFilename)

    # Update the blast database if updateFlag is true
    if(updateFlag):
        # Create the database from the file holding the subject sequences
        dbName = createBlastDB(subjectSequencesFilename,
            dbFilename)

    # BLAST query miRNAs to known miRNAs
    print("Performing BLAST")
    blastFilename = blastMirnas(subjectSequencesFilename,
        dbFilename, queryMirnasFilename)

    # Add field for the subject and query sequences in the BLAST output
    # because these sequences are not within by default
    print("Adding sequences to output file")
    addSequencesToOutput(queryMirnasFilename,
        subjectSequencesFilename)

    # Read the BLAST results into a dictionary for quick querying
    blastDict = readBlastResults(blastFilename)

    # Create a dictionary of either the known miRNAs that the candidate
    # miRNAs are equal to, or the miRNA families and species that the
    # candidate miRNA is highly similar to
    similarityDict = createSimilarityDict(blastDict,
        organism)

    numLibs = len(libFilenamesList)

    # Properly annotate the candidate miRNAs with the data in similarityDict
    annotateMirDeepCandidates(similarityDict, organism, mirBaseDict, numLibs)

if __name__ == "__main__":
    main()
