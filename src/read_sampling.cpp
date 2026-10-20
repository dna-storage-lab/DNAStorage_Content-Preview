#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

namespace {

struct FastqRecord {
    string header;
    string sequence;
    string separator;
    string quality;
};

void RemoveTrailingCarriageReturn(string& line) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
}

bool ReadFastqRecord(ifstream& input, FastqRecord& record, size_t record_number) {
    if (!getline(input, record.header)) {
        return false;
    }

    if (!getline(input, record.sequence) ||
        !getline(input, record.separator) ||
        !getline(input, record.quality)) {
        cerr << "Incomplete FASTQ record at read " << record_number << ".\n";
        exit(EXIT_FAILURE);
    }

    RemoveTrailingCarriageReturn(record.header);
    RemoveTrailingCarriageReturn(record.sequence);
    RemoveTrailingCarriageReturn(record.separator);
    RemoveTrailingCarriageReturn(record.quality);

    if (record.header.empty() || record.header[0] != '@') {
        cerr << "Invalid FASTQ header at read " << record_number << ".\n";
        exit(EXIT_FAILURE);
    }

    if (record.separator.empty() || record.separator[0] != '+') {
        cerr << "Invalid FASTQ separator at read " << record_number << ".\n";
        exit(EXIT_FAILURE);
    }

    if (record.sequence.size() != record.quality.size()) {
        cerr << "Sequence and quality lengths differ at read "
             << record_number << ".\n";
        exit(EXIT_FAILURE);
    }

    return true;
}

void WriteFastqRecord(ofstream& output, const FastqRecord& record) {
    output << record.header << '\n'
           << record.sequence << '\n'
           << record.separator << '\n'
           << record.quality << '\n';
}

string BuildOutputPath(const string& output_dir, size_t replicate) {
    ostringstream file_name;
    file_name << "replicate_" << setw(3) << setfill('0') << replicate << ".fastq";

    if (!output_dir.empty() &&
        (output_dir.back() == '/' || output_dir.back() == '\\')) {
        return output_dir + file_name.str();
    }

    return output_dir + "/" + file_name.str();
}

size_t ParsePositiveInteger(const char* value, const char* name) {
    char* end = NULL;
    const unsigned long parsed = strtoul(value, &end, 10);

    if (value[0] == '\0' || *end != '\0' || parsed == 0) {
        cerr << name << " must be a positive integer.\n";
        exit(EXIT_FAILURE);
    }

    return static_cast<size_t>(parsed);
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 6) {
        cerr << "Usage: " << argv[0]
             << " <input.fastq> <output_dir> <reads_per_replicate>"
             << " [replicate_count=100] [start_read=1]\n";
        return EXIT_FAILURE;
    }

    const string input_path = argv[1];
    const string output_dir = argv[2];
    const size_t reads_per_replicate =
        ParsePositiveInteger(argv[3], "reads_per_replicate");
    const size_t replicate_count = argc >= 5
        ? ParsePositiveInteger(argv[4], "replicate_count")
        : 100;
    const size_t start_read = argc >= 6
        ? ParsePositiveInteger(argv[5], "start_read")
        : 1;

    ifstream input(input_path.c_str(), ios::in | ios::binary);
    if (!input.is_open()) {
        cerr << "Failed to open FASTQ input: " << input_path << '\n';
        return EXIT_FAILURE;
    }

    FastqRecord record;
    size_t current_read = 1;

    for (; current_read < start_read; ++current_read) {
        if (!ReadFastqRecord(input, record, current_read)) {
            cerr << "The FASTQ file contains fewer than " << start_read
                 << " reads; sampling cannot start at the requested read.\n";
            return EXIT_FAILURE;
        }
    }

    size_t completed_replicates = 0;

    for (size_t replicate = 1; replicate <= replicate_count; ++replicate) {
        const string output_path = BuildOutputPath(output_dir, replicate);
        ofstream output(output_path.c_str(), ios::out | ios::binary);

        if (!output.is_open()) {
            cerr << "Failed to create sampled FASTQ: " << output_path << '\n';
            return EXIT_FAILURE;
        }

        bool complete = true;
        for (size_t i = 0; i < reads_per_replicate; ++i, ++current_read) {
            if (!ReadFastqRecord(input, record, current_read)) {
                complete = false;
                break;
            }
            WriteFastqRecord(output, record);
        }

        output.close();

        if (!complete) {
            remove(output_path.c_str());
            break;
        }

        ++completed_replicates;
    }

    cout << "Generated " << completed_replicates << " complete replicate(s), "
         << reads_per_replicate << " read(s) each, starting at FASTQ read "
         << start_read << ".\n";

    if (completed_replicates == 0) {
        cerr << "Not enough reads to generate one complete replicate.\n";
        return EXIT_FAILURE;
    }

    if (completed_replicates < replicate_count) {
        cerr << "Warning: requested " << replicate_count << " replicates but only "
             << completed_replicates << " complete replicates were available.\n";
    }

    return EXIT_SUCCESS;
}
