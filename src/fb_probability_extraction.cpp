#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "HalfMaker.h"

using namespace std;

struct ProcessResult {
    string new_read;
    vector<double> prob_1;
};

struct LineResult {
    string output_line;
    vector<double> prob_1;
    bool has_prob;
    string startLoc_in_raw_companion_sequence;
    string endLoc_in_raw_companion_sequence;
};

vector<string> SplitFields(const string& line) {
    istringstream iss(line);
    vector<string> fields;
    string field;

    while (iss >> field) {
        fields.push_back(field);
    }

    return fields;
}

ProcessResult ProcessTargetField(const string& companion_sequence_used, const string& regionSeq) {
    ProcessResult result;

    result.new_read = Companion_Sequence_Half_Maker_FB(
        static_cast<int>(companion_sequence_used.size()),
        companion_sequence_used,
        regionSeq,
        result.prob_1
    );

    return result;
}

void WriteProbabilityLine(
    ofstream& output,
    const string& startLoc_in_raw_companion_sequence,
    const string& endLoc_in_raw_companion_sequence,
    const vector<double>& prob_1,
    bool write_zero_prob
) {
    output << startLoc_in_raw_companion_sequence << '\t' << endLoc_in_raw_companion_sequence;

    for (size_t i = 0; i < prob_1.size(); ++i) {
        const double probability = write_zero_prob ? (1.0 - prob_1[i]) : prob_1[i];
        output << '\t' << probability;
    }

    output << '\n';
}

string JoinFields(const vector<string>& fields) {
    ostringstream joined;

    for (size_t i = 0; i < fields.size(); ++i) {
        if (i > 0) {
            joined << '\t';
        }
        joined << fields[i];
    }

    return joined.str();
}

LineResult ProcessOneLine(
    const string& line,
    size_t line_number,
    const string& input_path
) {
    LineResult line_result;
    line_result.has_prob = false;

    if (line.empty()) {
        return line_result;
    }

    vector<string> fields = SplitFields(line);
    const size_t required_field_count = 8;

    if (fields.size() < required_field_count) {
        cerr << "Line " << line_number << " in " << input_path
             << " has fewer than " << required_field_count
             << " fields; it will be copied without probabilities.\n";
        line_result.output_line = line;
        return line_result;
    }

    line_result.startLoc_in_raw_companion_sequence = fields[2];
    line_result.endLoc_in_raw_companion_sequence = fields[3];

    const string region_sequence = fields[6];
    string companion_sequence = fields[7];
    for (size_t i = 0; i < companion_sequence.size(); ++i) {
        companion_sequence[i] = companion_sequence[i] == 'A' ? '0' : '1';
    }

    const ProcessResult result = ProcessTargetField(companion_sequence, region_sequence);

    fields[6] = result.new_read;
    line_result.output_line = JoinFields(fields);
    line_result.prob_1 = result.prob_1;
    line_result.has_prob = true;
    return line_result;
}

void ProcessLineRange(
    const vector<string>& lines,
    vector<LineResult>& results,
    size_t start,
    size_t end,
    const string& input_path
) {
    for (size_t index = start; index < end; ++index) {
        results[index] = ProcessOneLine(lines[index], index + 1, input_path);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        cerr << "Usage: " << argv[0]
             << " <matched_reads.txt> <output_prefix> [thread_count]\n";
        return 1;
    }

    const string input_path = argv[1];
    const string output_prefix = argv[2];
    const string output_path = output_prefix + "_fb_result.txt";
    const string prob_1_path = output_prefix + "_prob1.txt";
    const string prob_0_path = output_prefix + "_prob0.txt";

    unsigned int thread_count = thread::hardware_concurrency();
    if (thread_count == 0) {
        thread_count = 1;
    }

    if (argc == 4) {
        istringstream parser(argv[3]);
        unsigned int requested_threads = 0;
        parser >> requested_threads;
        if (!parser || !parser.eof() || requested_threads == 0) {
            cerr << "thread_count must be a positive integer.\n";
            return 1;
        }
        thread_count = requested_threads;
    }

    ifstream input(input_path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open matched-read input: " << input_path << '\n';
        return 1;
    }

    vector<string> lines;
    string line;
    while (getline(input, line)) {
        lines.push_back(line);
    }

    cout << "Processing " << lines.size() << " aligned read(s) with up to "
         << thread_count << " thread(s).\n";

    vector<LineResult> results(lines.size());
    const size_t real_thread_count = lines.empty()
        ? 0
        : min(static_cast<size_t>(thread_count), lines.size());
    vector<thread> workers;

    if (real_thread_count > 0) {
        const size_t block_size =
            (lines.size() + real_thread_count - 1) / real_thread_count;

        for (size_t t = 0; t < real_thread_count; ++t) {
            const size_t start = t * block_size;
            const size_t end = min(start + block_size, lines.size());

            if (start >= end) {
                break;
            }

            workers.emplace_back(
                ProcessLineRange,
                cref(lines),
                ref(results),
                start,
                end,
                cref(input_path)
            );
        }
    }

    for (size_t t = 0; t < workers.size(); ++t) {
        workers[t].join();
    }

    ofstream corrected_output(output_path.c_str());
    ofstream prob_1_output(prob_1_path.c_str());
    ofstream prob_0_output(prob_0_path.c_str());

    if (!corrected_output.is_open() ||
        !prob_1_output.is_open() ||
        !prob_0_output.is_open()) {
        cerr << "Failed to create one or more Step 3 outputs for prefix: "
             << output_prefix << '\n';
        return 1;
    }

    prob_1_output << fixed << setprecision(10);
    prob_0_output << fixed << setprecision(10);

    for (size_t index = 0; index < results.size(); ++index) {
        corrected_output << results[index].output_line << '\n';

        if (results[index].has_prob) {
            WriteProbabilityLine(
                prob_1_output,
                results[index].startLoc_in_raw_companion_sequence,
                results[index].endLoc_in_raw_companion_sequence,
                results[index].prob_1,
                false
            );
            WriteProbabilityLine(
                prob_0_output,
                results[index].startLoc_in_raw_companion_sequence,
                results[index].endLoc_in_raw_companion_sequence,
                results[index].prob_1,
                true
            );
        }
        else {
            prob_1_output << '\n';
            prob_0_output << '\n';
        }
    }

    cout << "Generated: " << output_path << '\n'
         << "Generated: " << prob_1_path << '\n'
         << "Generated: " << prob_0_path << '\n';
    return 0;
}
