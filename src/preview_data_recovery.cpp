#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::fprintf(
            stderr,
            "Usage: %s <information_bits.txt> <recovered_file> <bit_length>\n",
            argv[0]
        );
        return EXIT_FAILURE;
    }

    char* end = NULL;
    const long parsed_length = std::strtol(argv[3], &end, 10);
    if (argv[3][0] == '\0' || *end != '\0' || parsed_length <= 0 ||
        parsed_length % 8 != 0) {
        std::fprintf(stderr, "bit_length must be a positive multiple of 8.\n");
        return EXIT_FAILURE;
    }

    const size_t bit_count = static_cast<size_t>(parsed_length);
    const size_t byte_count = bit_count / 8;
    std::vector<char> bits(bit_count);
    std::vector<unsigned char> recovered(byte_count, 0);

    FILE* input = std::fopen(argv[1], "rb");
    if (input == NULL) {
        std::fprintf(stderr, "Failed to open information-bit input: %s\n", argv[1]);
        return EXIT_FAILURE;
    }

    const size_t bits_read = std::fread(&bits[0], sizeof(char), bit_count, input);
    std::fclose(input);
    if (bits_read != bit_count) {
        std::fprintf(stderr, "Information-bit input is shorter than %zu bits.\n", bit_count);
        return EXIT_FAILURE;
    }

    for (size_t byte_index = 0; byte_index < byte_count; ++byte_index) {
        for (size_t bit_index = 0; bit_index < 8; ++bit_index) {
            const char bit = bits[byte_index * 8 + bit_index];
            if (bit != '0' && bit != '1') {
                std::fprintf(stderr, "Invalid bit at position %zu.\n", byte_index * 8 + bit_index);
                return EXIT_FAILURE;
            }
            recovered[byte_index] |= static_cast<unsigned char>(
                (bit - '0') << (7 - bit_index)
            );
        }
    }

    FILE* output = std::fopen(argv[2], "wb");
    if (output == NULL) {
        std::fprintf(stderr, "Failed to create recovered output: %s\n", argv[2]);
        return EXIT_FAILURE;
    }

    const size_t bytes_written =
        std::fwrite(&recovered[0], sizeof(unsigned char), byte_count, output);
    std::fclose(output);
    if (bytes_written != byte_count) {
        std::fprintf(stderr, "Failed to write the complete recovered output.\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
