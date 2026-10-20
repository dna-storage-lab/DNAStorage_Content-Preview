#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>

using namespace std;

struct ProbSegment {
    long long start_pos;
    long long end_pos;
    vector<double> probs;
};

vector<ProbSegment> read_prob_file(const string& filename) {
    ifstream fin(filename);

    if (!fin.is_open()) {
        cerr << "Cannot open file: " << filename << endl;
        exit(1);
    }

    vector<ProbSegment> data;
    string line;
    long long line_number = 0;

    while (getline(fin, line)) {
        line_number++;

        if (line.empty()) {
            continue;
        }

        stringstream ss(line);

        ProbSegment segment;

        if (!(ss >> segment.start_pos >> segment.end_pos)) {
            cerr << "Invalid line format in " << filename
                 << " at line " << line_number
                 << ". Each line should start with start_pos and end_pos."
                 << endl;
            exit(1);
        }

        double val;

        while (ss >> val) {
            segment.probs.push_back(val);
        }

        if (segment.probs.empty()) {
            cerr << "No probability values found in " << filename
                 << " at line " << line_number << endl;
            exit(1);
        }

        if (segment.start_pos < 0 || segment.end_pos < 0) {
            cerr << "Negative start/end position in " << filename
                 << " at line " << line_number << endl;
            exit(1);
        }

        if (segment.end_pos < segment.start_pos) {
            cerr << "end_pos < start_pos in " << filename
                 << " at line " << line_number << endl;
            exit(1);
        }

        /*
            Retain the length check here.
            The actual probability mapping depends only on:
            absolute_pos = start_pos + offset

            Therefore, whether end_pos denotes a closed or open interval does not
            affect merging as long as the number of probabilities is correct.
        */
        long long expected_len_closed = segment.end_pos - segment.start_pos + 1;

        if (expected_len_closed != static_cast<long long>(segment.probs.size())) {
            cerr << "Warning: probability length does not match closed interval start/end at line "
                 << line_number << " in " << filename << ". "
                 << "start=" << segment.start_pos
                 << ", end=" << segment.end_pos
                 << ", expected closed length=" << expected_len_closed
                 << ", actual length=" << segment.probs.size()
                 << ". The program will map probabilities by start_pos + offset."
                 << endl;
        }

        data.push_back(segment);
    }

    fin.close();

    return data;
}

long long parse_positive_length(const string& s) {
    char* end_ptr = nullptr;
    long long value = strtoll(s.c_str(), &end_ptr, 10);

    if (end_ptr == s.c_str() || *end_ptr != '\0' || value <= 0) {
        cerr << "Invalid ideal companion-sequence length: " << s << endl;
        exit(1);
    }

    return value;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        cerr << "Usage: ./fb_decode_MajorityVoting_v5 prob0.txt prob1.txt output.txt ideal_companion_sequence_length" << endl;
        cerr << "Example: ./fb_decode_MajorityVoting_v5 prob0.txt prob1.txt output.txt 1000" << endl;
        return 1;
    }

    string prob0_file = argv[1];
    string prob1_file = argv[2];
    string output_file = argv[3];

    long long ideal_companion_sequence_length_ll = parse_positive_length(argv[4]);
    size_t ideal_companion_sequence_length = static_cast<size_t>(ideal_companion_sequence_length_ll);

    vector<ProbSegment> prob0 = read_prob_file(prob0_file);
    vector<ProbSegment> prob1 = read_prob_file(prob1_file);

    if (prob0.empty() || prob1.empty()) {
        cerr << "Input probability file is empty." << endl;
        return 1;
    }

    if (prob0.size() != prob1.size()) {
        cerr << "Row number mismatch between prob0 and prob1 files." << endl;
        return 1;
    }

    int row_num = static_cast<int>(prob0.size());

    /*
        Key change:
        Initialize the arrays directly from the ideal companion-sequence length.
        Keep the output length fixed at ideal_companion_sequence_length.
    */
    vector<double> log_p0_all(ideal_companion_sequence_length, 0.0);
    vector<double> log_p1_all(ideal_companion_sequence_length, 0.0);

    vector<bool> zero_p0_all(ideal_companion_sequence_length, false);
    vector<bool> zero_p1_all(ideal_companion_sequence_length, false);

    vector<int> valid_read_count_all(ideal_companion_sequence_length, 0);

    long long ignored_out_of_range_count = 0;

    for (int r = 0; r < row_num; r++) {
        if (prob0[r].start_pos != prob1[r].start_pos ||
            prob0[r].end_pos != prob1[r].end_pos) {
            cerr << "Start/end position mismatch between prob0 and prob1 at row "
                 << r << endl;
            return 1;
        }

        if (prob0[r].probs.size() != prob1[r].probs.size()) {
            cerr << "Probability length mismatch between prob0 and prob1 at row "
                 << r << endl;
            return 1;
        }

        long long start_pos = prob0[r].start_pos;
        size_t prob_len = prob0[r].probs.size();

        for (size_t k = 0; k < prob_len; k++) {
            long long abs_pos_ll = start_pos + static_cast<long long>(k);

            /*
                Key change:
                Skip probability positions outside the ideal companion-sequence range.
                Merge only positions 0 through ideal_companion_sequence_length - 1.
            */
            if (abs_pos_ll < 0 || abs_pos_ll >= ideal_companion_sequence_length_ll) {
                ignored_out_of_range_count++;
                continue;
            }

            size_t abs_pos = static_cast<size_t>(abs_pos_ll);

            double p0 = prob0[r].probs[k];
            double p1 = prob1[r].probs[k];

            valid_read_count_all[abs_pos]++;

            if (p0 <= 0) {
                zero_p0_all[abs_pos] = true;
            } else {
                log_p0_all[abs_pos] += log(p0);
            }

            if (p1 <= 0) {
                zero_p1_all[abs_pos] = true;
            } else {
                log_p1_all[abs_pos] += log(p1);
            }
        }
    }

    ofstream fout(output_file);

    if (!fout.is_open()) {
        cerr << "Cannot create output file: " << output_file << endl;
        return 1;
    }

    string result;
    result.reserve(ideal_companion_sequence_length);

    int uncovered_count = 0;

    /*
        Key change:
        Fix the output range at 0 through ideal_companion_sequence_length - 1.
        Fill positions without probability coverage with 0.
    */
    for (size_t pos = 0; pos < ideal_companion_sequence_length; pos++) {
        if (valid_read_count_all[pos] == 0) {
            result += '0';
            uncovered_count++;
            continue;
        }

        double log_p0 = log_p0_all[pos];
        double log_p1 = log_p1_all[pos];

        if (zero_p0_all[pos]) {
            log_p0 = -numeric_limits<double>::infinity();
        }

        if (zero_p1_all[pos]) {
            log_p1 = -numeric_limits<double>::infinity();
        }

        if (log_p1 > log_p0) {
            result += '1';
        } else {
            result += '0';
        }
    }

    fout << result << endl;
    fout.close();

    cout << "Ideal companion-sequence length: " << ideal_companion_sequence_length << endl;
    cout << "Output length: " << result.size() << endl;
    cout << "Output range: 0-" << ideal_companion_sequence_length - 1 << endl;
    cout << "Uncovered positions filled with 0: " << uncovered_count << endl;
    cout << "Ignored out-of-range probability values: " << ignored_out_of_range_count << endl;
    cout << "Output file: " << output_file << endl;

    return 0;
}
