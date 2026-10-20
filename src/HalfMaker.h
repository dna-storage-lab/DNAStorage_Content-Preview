#ifndef HALF_MAKER_H
#define HALF_MAKER_H

#include <string>
#include <vector>

std::string Companion_Sequence_Half_Maker_FB(
    int T,
    std::string companion_sequence,
    std::string indel_base_sequence,
    std::vector<double>& llr_1
);

#endif
