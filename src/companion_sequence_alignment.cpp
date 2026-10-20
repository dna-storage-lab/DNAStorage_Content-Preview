#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cctype>
#include <cstdlib>

#include "edlib.h"

using namespace std;

#define DUMP_ALIGN_RESULT 0

const double threshold = 0.1;

// Convert a base to a lower-layer companion-sequence symbol:
// A/G -> A, T/C -> T; if the input is already 0/1, 0 -> A and 1 -> T.
// Treat N and other ambiguous IUPAC characters as N.
char toLowerLayerSymbol(char c)
{
    c = (char)toupper((unsigned char)c);
    if (c == 'A' || c == 'G' || c == '0')
        return 'A';
    if (c == 'T' || c == 'C' || c == '1')
        return 'T';
    if (c == 'N')
        return 'N';
    return '\0';
}

bool isSequenceLikeArgument(const string &s)
{
    for (char c : s)
    {
        if (isspace((unsigned char)c))
            continue;
        char u = (char)toupper((unsigned char)c);
        if (!(u == 'A' || u == 'T' || u == 'G' || u == 'C' || u == 'N' || u == '0' || u == '1'))
            return false;
    }
    return true;
}

// Read the lower-layer companion sequence:
// 1) First try argv[1] as a file path.
// 2) If the file cannot be opened but argv[1] looks like an A/T/0/1 sequence,
//    read it as a raw sequence string.
// Do not split FASTA input into upper and lower layers; produce one companion sequence.
string readCompanionSequence(const char *input)
{
    string seq;
    ifstream fin(input);

    if (fin.good())
    {
        string line;
        while (getline(fin, line))
        {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();

            // If FASTA input is supplied, skip headers and merge all sequence lines.
            if (line.empty() || line[0] == '>')
                continue;

            for (char c : line)
            {
                if (isspace((unsigned char)c))
                    continue;
                char x = toLowerLayerSymbol(c);
                if (x != '\0')
                    seq.push_back(x);
            }
        }
        fin.close();
    }
    else
    {
        string raw(input);
        if (!isSequenceLikeArgument(raw))
        {
            fprintf(stderr, "Failed to open %s! If argv[1] is a raw sequence, use only A/T/G/C/N or 0/1.\n", input);
            exit(1);
        }

        for (char c : raw)
        {
            if (isspace((unsigned char)c))
                continue;
            char x = toLowerLayerSymbol(c);
            if (x != '\0')
                seq.push_back(x);
        }
    }

    if (seq.empty())
    {
        fprintf(stderr, "No valid companion sequence was found in %s!\n", input);
        exit(1);
    }

    return seq;
}

// Reverse-complement a lower-layer companion sequence: A <-> T, N -> N.
// When A/T represent 0/1, the reverse strand must be reversed and complemented.
string reverseComplementLowerLayer(const string &in)
{
    string out;
    out.resize(in.size());

    for (size_t i = 0; i < in.size(); ++i)
    {
        char c = in[in.size() - 1 - i];
        if (c == 'A')
            out[i] = 'T';
        else if (c == 'T')
            out[i] = 'A';
        else
            out[i] = 'N';
    }

    return out;
}

// Reverse-complement the raw DNA read to output it in companion-sequence orientation.
string reverseComplementDNA(const string &in)
{
    string out;
    out.resize(in.size());

    for (size_t i = 0; i < in.size(); ++i)
    {
        char c = (char)toupper((unsigned char)in[in.size() - 1 - i]);
        if (c == 'A')
            out[i] = 'T';
        else if (c == 'T')
            out[i] = 'A';
        else if (c == 'G')
            out[i] = 'C';
        else if (c == 'C')
            out[i] = 'G';
        else
            out[i] = 'N';
    }

    return out;
}

string convertReadToLowerLayer(const string &read)
{
    string lower;
    lower.reserve(read.size());

    for (char c : read)
    {
        char x = toLowerLayerSymbol(c);
        if (x == '\0')
            x = 'N';
        lower.push_back(x);
    }

    return lower;
}

// Preserve the original INDEL-correction output behavior:
// match: copy the raw read base;
// op 1: output N;
// op 2: skip one base in the read;
// mismatch: output N and skip one base in the read.
string buildModifiedSequence(const string &read, int startLoc, const EdlibAlignResult &alignRst)
{
    string modified;
    modified.reserve((size_t)alignRst.alignmentLength + 1);

    int m = 0;
    for (int ti = 0; ti < alignRst.alignmentLength; ++ti)
    {
        unsigned char op = alignRst.alignment[ti];
        int pos = startLoc + m;

        if (op == 0) // match
        {
            if (pos >= 0 && pos < (int)read.size())
                modified.push_back(read[(size_t)pos]);
            else
                modified.push_back('N');
            m++;
        }
        else if (op == 1)
        {
            modified.push_back('N');
        }
        else if (op == 2)
        {
            m++;
        }
        else if (op == 3) // mismatch
        {
            modified.push_back('N');
            m++;
        }
    }

    return modified;
}

// Use this when the read is the query and the companion sequence is the target.
// This occurs when the read is shorter than the companion sequence.
// Edlib start/end values are then companion-sequence coordinates and cannot crop the read.
string buildModifiedSequenceReadAsQuery(const string &queryRead, const EdlibAlignResult &alignRst)
{
    string modified;
    modified.reserve((size_t)alignRst.alignmentLength + 1);

    int q = 0;
    for (int ti = 0; ti < alignRst.alignmentLength; ++ti)
    {
        unsigned char op = alignRst.alignment[ti];

        if (op == 0) // match
        {
            if (q >= 0 && q < (int)queryRead.size())
                modified.push_back(queryRead[(size_t)q]);
            else
                modified.push_back('N');
            q++;
        }
        else if (op == 1)
        {
            // Use N when the companion sequence has a position with no read base.
            modified.push_back('N');
        }
        else if (op == 2)
        {
            // Do not write extra read bases into the companion-sequence coordinates.
            q++;
        }
        else if (op == 3) // mismatch
        {
            modified.push_back('N');
            q++;
        }
    }

    return modified;
}

int main(int argc, char *argv[])
{
    time_t startTime, endTime;
    char timebuff[100];
    time(&startTime);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&startTime));
    // fprintf(stderr, "Running program: %s\n", timebuff);

    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s [Companion sequence/file] [reads FASTQ] [Match reads file]\n", argv[0]);
        fprintf(stderr, "Example 1: %s companion_sequence.txt reads.fastq matched_reads.txt\n", argv[0]);
        fprintf(stderr, "Example 2: %s ATTAATTTAA reads.fastq matched_reads.txt\n", argv[0]);
        exit(1);
    }

    const char *companion_sequence_input = argv[1];
    const char *file_fq = argv[2];
    const char *file_match = argv[3];

    //
    string companion_sequence = readCompanionSequence(companion_sequence_input);
    string companion_sequence_rc = reverseComplementLowerLayer(companion_sequence);
    int companion_sequence_length = (int)companion_sequence.size();

    // 
    printf("Number of companion sequences: 1\n");
    printf("Companion-sequence length: %d\n", companion_sequence_length);

    ifstream fpr_fq(file_fq);
    if (!fpr_fq.good())
    {
        fprintf(stderr, "Failed to open %s!\n", file_fq);
        exit(1);
    }

    ofstream fpw_match(file_match);
    if (!fpw_match.good())
    {
        fprintf(stderr, "Failed to open %s!\n", file_match);
        exit(1);
    }

#if DUMP_ALIGN_RESULT
    ofstream fpw_rslt("alignment_rslt_edlib.txt");
    if (!fpw_rslt.good())
    {
        fprintf(stderr, "Failed to open alignment_rslt_edlib.txt!\n");
        exit(1);
    }
#endif

    string header, curReads, plusLine, qualLine;
    int cnt = 0;
    int alignedNum = 0;

    while (getline(fpr_fq, header))
    {
        if (!getline(fpr_fq, curReads))
            break;
        if (!getline(fpr_fq, plusLine))
            break;
        if (!getline(fpr_fq, qualLine))
            break;

        if (!curReads.empty() && curReads.back() == '\r')
            curReads.pop_back();

        cnt++;


        int curLen = (int)curReads.size();
        // if (curLen < companion_sequence_length)
        //     continue;

        // Extract the lower-layer companion sequence from the FASTQ read: A/G -> A, T/C -> T.
        string lowerRead = convertReadToLowerLayer(curReads);

        // Match the behavior of the first Python program:
        // 1) if read_length < companion_sequence_length, align the read or its reverse complement
        //    against the companion sequence;
        // 2) otherwise, align the companion sequence against the read or its reverse complement.
        bool targetIsRead = false;
        int targetLen = 0;
        EdlibAlignResult alignRstF;
        EdlibAlignResult alignRstR;

        string lowerReadRC = reverseComplementLowerLayer(lowerRead);// Reverse-complement the read.

        if (curLen < companion_sequence_length)
        {
            

            alignRstF = edlibAlign(
                lowerRead.c_str(), curLen,
                companion_sequence.c_str(), companion_sequence_length,
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            alignRstR = edlibAlign(
                lowerReadRC.c_str(), curLen,
                companion_sequence.c_str(), companion_sequence_length,
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            targetIsRead = false; // start/end are companion-sequence coordinates.
            targetLen = companion_sequence_length;
        }
        else
        {
            alignRstF = edlibAlign(
                companion_sequence.c_str(), companion_sequence_length,
                lowerRead.c_str(), curLen,
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            alignRstR = edlibAlign(
                companion_sequence.c_str(), companion_sequence_length,
                lowerReadRC.c_str(), curLen,
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            targetIsRead = true; // start/end are read coordinates.
            targetLen = curLen;
        }

        if (alignRstF.status != EDLIB_STATUS_OK || alignRstR.status != EDLIB_STATUS_OK)
        {
            edlibFreeAlignResult(alignRstF);
            edlibFreeAlignResult(alignRstR);
            continue;
        }

        int fminED = alignRstF.editDistance;
        int rminED = alignRstR.editDistance;

        int optimal_companion_sequence_ID = 0; // There is one companion sequence, so its ID is 0.
        int optimal_ED = 0;
        int optimal_startLoc = 0;
        int optimal_endLoc = 0;
        int optimal_isForward = 1;
        int startLoc_in_raw_companion_sequence=0;
        int endLoc_in_raw_companion_sequence=0;
        const EdlibAlignResult *optimalAlign = NULL;

        if (fminED <= rminED)
        {
            optimal_ED = fminED;
            optimal_isForward = 1;
            optimal_startLoc = alignRstF.startLocations ? alignRstF.startLocations[0] : -1;
            optimal_endLoc = alignRstF.endLocations ? alignRstF.endLocations[0] : -1;
            optimalAlign = &alignRstF;
        }
        else
        {
            optimal_ED = rminED;
            optimal_isForward = 0;
            optimal_startLoc = alignRstR.startLocations ? alignRstR.startLocations[0] : -1;
            optimal_endLoc = alignRstR.endLocations ? alignRstR.endLocations[0] : -1;
            optimalAlign = &alignRstR;
        }

        if (optimal_ED >= 0 &&
            optimal_startLoc >= 0 && optimal_endLoc >= optimal_startLoc && optimal_endLoc < targetLen)
        {
            int actual_aligned_length = optimal_endLoc - optimal_startLoc + 1;
            double edit_rate = (double)optimal_ED / (double)actual_aligned_length;

            if (edit_rate <= threshold)
            {
                alignedNum++;

                string regionSeq;
                string matchSeq_modified;
                string companion_sequence_used;

                if (targetIsRead)
                {
                    startLoc_in_raw_companion_sequence= 0;
                    endLoc_in_raw_companion_sequence= companion_sequence_length - 1;
                    // When the companion sequence is aligned to the read, start/end are read coordinates.
                    if (optimal_isForward)// Keep the cropped read as the matched region on the forward strand.
                    {
                        regionSeq = curReads.substr((size_t)optimal_startLoc,
                                                    (size_t)actual_aligned_length);
                        matchSeq_modified = buildModifiedSequence(curReads, optimal_startLoc, *optimalAlign);
                        companion_sequence_used = companion_sequence;
                       
                    }
                    else
                    {
                        regionSeq = reverseComplementDNA(curReads).substr((size_t)optimal_startLoc,
                                                    (size_t)actual_aligned_length);
                        matchSeq_modified = buildModifiedSequence(reverseComplementDNA(curReads), optimal_startLoc, *optimalAlign);
                        companion_sequence_used = companion_sequence;
                        
                    }
                }
                else
                {
                    startLoc_in_raw_companion_sequence= optimal_startLoc;
                    endLoc_in_raw_companion_sequence= optimal_endLoc;
                    // When the read is aligned to the companion sequence, start/end are companion-sequence coordinates.
                    // Keep the full read as the matched region and correct reverse-strand reads in reverse-complement orientation.
                    companion_sequence_used = companion_sequence;
                    if (optimal_isForward)
                    {
                        matchSeq_modified = buildModifiedSequenceReadAsQuery(curReads, *optimalAlign);
                        regionSeq = curReads;
                        companion_sequence_used = companion_sequence.substr((size_t)optimal_startLoc, 
                                                (size_t)actual_aligned_length);// Output the cropped companion sequence.
                    }
                        
                    else
                    {
                        regionSeq = reverseComplementDNA(curReads);// Use the reverse-complemented read as the matched region.
                        matchSeq_modified = buildModifiedSequenceReadAsQuery(reverseComplementDNA(curReads), *optimalAlign);
                        companion_sequence_used = companion_sequence.substr((size_t)optimal_startLoc, 
                                                (size_t)actual_aligned_length);// Output the cropped companion sequence.
                    }
                        
                }

                fpw_match << optimal_companion_sequence_ID << '\t'
                          << optimal_ED << '\t'
                          <<startLoc_in_raw_companion_sequence<< '\t'
                          <<endLoc_in_raw_companion_sequence<< '\t'
                          << optimal_isForward << '\t'
                          << matchSeq_modified << '\t'
                          << regionSeq << '\t'
                          << companion_sequence_used << '\n';
                
                //  fpw_match << optimal_companion_sequence_ID << '\t'
                //           << optimal_ED << '\t'
                //           << optimal_isForward << '\t'
                //           << matchSeq_modified << '\t'
                //           << regionSeq << '\n';

#if DUMP_ALIGN_RESULT
                fpw_rslt << (cnt - 1) << '\t'
                         << optimal_companion_sequence_ID << '\t'
                         << optimal_ED << '\t'
                         << optimal_startLoc << '\t'
                         << optimal_endLoc << '\t'
                         << optimal_isForward << '\t'
                         << (targetIsRead ? "target_is_read" : "target_is_companion_sequence") << '\n';
#endif
            }
        }

        edlibFreeAlignResult(alignRstF);
        edlibFreeAlignResult(alignRstR);
    }

    fpr_fq.close();
    fpw_match.close();

#if DUMP_ALIGN_RESULT
    fpw_rslt.close();
#endif

    printf("\nTotal number of reads processed: %d\n", cnt);
    printf("Total number of fragment aligned: %d\n", alignedNum);
    fflush(stdout);

    time(&endTime);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&endTime));
    // fprintf(stderr, "\nFinished running: %s\n", timebuff);

    return 0;
}
