# Network-Analysis-Tool
book_of_documents_t_s.txt is the written vernacular version in simplified Chinese.
book_of_documents.txt is the original version in Traditional Chinese.

network_generation.cpp works for the "character" network.
network_generation_w.cpp works for the "word" network.

Enable or Disable corresponding blocks of code to perform the desired computations.

To compile: g++ [the cpp file] -std=c++17 -Wall -O3 -o [output file name]
For example: g++ network_generation.cpp -std=c++17 -Wall -O3 -o output_runnable

To run: [output file name] [the text]
For example: ./output_runnable book_of_documents
